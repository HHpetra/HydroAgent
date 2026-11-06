"""SCS Curve Number excess rainfall, using the HEC-HMS cumulative method.

Depths are millimeters. Timestamp labels are RFC3339 strings. Nothing in this
module reads the wall clock; a cut such as ``as_of`` is an argument.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from enum import StrEnum
from typing import TypeAlias

_S_NUMERATOR_MM = 25400.0
_S_OFFSET_MM = 254.0
_COMPOSITE_IMPERVIOUS_CN = 98.0
_BALANCE_TOLERANCE_MM = 1e-9

_RFC3339 = re.compile(
    r"^\d{4}-\d{2}-\d{2}[Tt]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:[Zz]|[+-]\d{2}:\d{2})$"
)

class ScsCnErrorCode(StrEnum):
    """Stable identifier for an SCS-CN domain failure."""

    INVALID_CURVE_NUMBER = "invalid_curve_number"
    INVALID_LAMBDA = "invalid_lambda"
    INVALID_IMPERVIOUS_FRACTION = "invalid_impervious_fraction"
    INVALID_AMC_CONDITION = "invalid_amc_condition"
    INVALID_IMPERVIOUS_TREATMENT = "invalid_impervious_treatment"
    EMPTY_RAINFALL = "empty_rainfall"
    INVALID_RAINFALL = "invalid_rainfall"
    INVALID_TIMESTAMP = "invalid_timestamp"
    UNORDERED_TIMESTAMPS = "unordered_timestamps"
    WATER_BALANCE = "water_balance"


class ScsCnError(Exception):
    """Illegal SCS-CN input, or a result that breaks water balance.

    ``code`` identifies the failure. ``message`` describes the offending value.
    """

    def __init__(self, code: ScsCnErrorCode, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code.value}: {message}")


class AmcCondition(StrEnum):
    """Antecedent moisture condition for a published curve number CN II."""

    AMC_I = "AMC_I"
    AMC_II = "AMC_II"
    AMC_III = "AMC_III"


class ImperviousTreatment(StrEnum):
    """How directly connected impervious area enters the excess calculation.

    ``SPLIT`` is the HEC-HMS percent-impervious method: rainfall on that
    fraction becomes excess immediately, and SCS-CN applies to the remainder.

    ``COMPOSITE_CN`` folds the fraction into one curve number,
    ``CN + f * (98 - CN)``, then applies SCS-CN to the whole area.
    """

    SPLIT = "split"
    COMPOSITE_CN = "composite_cn"


def _require_real(value: object, code: ScsCnErrorCode, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ScsCnError(code, f"{name} must be a real number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise ScsCnError(code, f"{name} must be finite, got {value!r}")
    return number


def _require_curve_number(value: object) -> float:
    number = _require_real(value, ScsCnErrorCode.INVALID_CURVE_NUMBER, "cn_ii")
    if not 0.0 < number <= 100.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_CURVE_NUMBER,
            f"cn_ii must be in (0, 100], got {value!r}",
        )
    return number


def _require_lambda(value: object) -> float:
    number = _require_real(value, ScsCnErrorCode.INVALID_LAMBDA, "lambda_coeff")
    if not 0.0 <= number <= 1.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_LAMBDA,
            f"lambda_coeff must be in [0, 1], got {value!r}",
        )
    return number


def _require_impervious_fraction(value: object) -> float:
    number = _require_real(
        value,
        ScsCnErrorCode.INVALID_IMPERVIOUS_FRACTION,
        "impervious_fraction",
    )
    if not 0.0 <= number <= 1.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_IMPERVIOUS_FRACTION,
            f"impervious_fraction must be in [0, 1], got {value!r}",
        )
    return number


def _require_amc(value: object) -> AmcCondition:
    if isinstance(value, AmcCondition):
        return value
    if isinstance(value, str):
        try:
            return AmcCondition(value)
        except ValueError:
            pass
    raise ScsCnError(
        ScsCnErrorCode.INVALID_AMC_CONDITION,
        f"amc_condition must be AMC_I, AMC_II, or AMC_III, got {value!r}",
    )


def _require_treatment(value: object) -> ImperviousTreatment:
    if isinstance(value, ImperviousTreatment):
        return value
    if isinstance(value, str):
        try:
            return ImperviousTreatment(value)
        except ValueError:
            pass
    raise ScsCnError(
        ScsCnErrorCode.INVALID_IMPERVIOUS_TREATMENT,
        "impervious_treatment must be 'split' or 'composite_cn', "
        f"got {value!r}",
    )


def _require_depth(value: object, index: int | None) -> float:
    where = "rainfall depth" if index is None else f"rainfall depth at index {index}"
    number = _require_real(value, ScsCnErrorCode.INVALID_RAINFALL, where)
    if number < 0.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_RAINFALL,
            f"{where} must be >= 0 mm, got {value!r}",
        )
    return number


def _parse_rfc3339(value: object, *, role: str) -> datetime:
    if not isinstance(value, str) or _RFC3339.fullmatch(value) is None:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_TIMESTAMP,
            f"{role} must be an RFC3339 timestamp, got {value!r}",
        )
    normalized = f"{value[:-1]}+00:00" if value[-1] in "Zz" else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_TIMESTAMP,
            f"{role} is not a valid RFC3339 timestamp: {value}",
        ) from exc
    if parsed.tzinfo is None:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_TIMESTAMP,
            f"{role} must include a timezone offset, got {value!r}",
        )
    return parsed


@dataclass(frozen=True, slots=True)
class ScsCnParameters:
    """SCS-CN loss parameters for one storm.

    ``cn_ii`` is the curve number for average antecedent moisture. Pass
    ``AMC_II`` when ``cn_ii`` is already the curve number to use. ``lambda_coeff``
    is the initial-abstraction ratio Ia / S (0.20 is the SCS standard).
    """

    cn_ii: float
    amc_condition: AmcCondition | str = AmcCondition.AMC_II
    lambda_coeff: float = 0.2
    impervious_fraction: float = 0.0
    impervious_treatment: ImperviousTreatment | str = ImperviousTreatment.SPLIT

    def __post_init__(self) -> None:
        object.__setattr__(self, "cn_ii", _require_curve_number(self.cn_ii))
        object.__setattr__(self, "amc_condition", _require_amc(self.amc_condition))
        object.__setattr__(self, "lambda_coeff", _require_lambda(self.lambda_coeff))
        object.__setattr__(
            self,
            "impervious_fraction",
            _require_impervious_fraction(self.impervious_fraction),
        )
        object.__setattr__(
            self,
            "impervious_treatment",
            _require_treatment(self.impervious_treatment),
        )


@dataclass(frozen=True, slots=True)
class TimedRainfall:
    """Rainfall depth labelled with the original RFC3339 timestamp."""

    time: str
    depth_mm: float

    def __post_init__(self) -> None:
        _parse_rfc3339(self.time, role="rainfall time")
        object.__setattr__(self, "depth_mm", _require_depth(self.depth_mm, index=None))


RainfallSeries: TypeAlias = Sequence[int | float] | Sequence[TimedRainfall]


@dataclass(frozen=True, slots=True)
class ExcessStep:
    """One timestep of cumulative SCS-CN accounting, in millimetres."""

    time: str | None
    rainfall_mm: float
    cumulative_rainfall_mm: float
    excess_mm: float
    cumulative_excess_mm: float
    loss_mm: float
    cumulative_initial_abstraction_mm: float
    cumulative_infiltration_mm: float


@dataclass(frozen=True, slots=True)
class ScsCnResult:
    """Excess rainfall for one storm.

    ``curve_number`` is the value substituted into the SCS equation after
    antecedent-moisture and, when requested, composite-impervious adjustment.
    ``initial_abstraction_mm`` is the threshold Ia = λS, not the depth actually
    consumed. Consumed initial abstraction and infiltration are watershed-average
    depths, so they close with excess against rainfall.
    """

    curve_number: float
    potential_retention_mm: float
    initial_abstraction_mm: float
    total_rainfall_mm: float
    total_excess_mm: float
    initial_abstraction_consumed_mm: float
    cumulative_infiltration_mm: float
    runoff_coefficient: float | None
    steps: tuple[ExcessStep, ...]


@dataclass(frozen=True, slots=True)
class RemainderExcess:
    """Incremental excess strictly after ``t_cut``, taken from the full storm.

    Cumulative fields on ``steps`` still start at the beginning of the storm.
    ``storm`` is that full computation.
    """

    t_cut: str
    remainder_rainfall_mm: float
    remainder_excess_mm: float
    steps: tuple[ExcessStep, ...]
    storm: ScsCnResult


@dataclass(frozen=True, slots=True)
class ScsCnCase:
    """One rainfall series paired with one parameter set."""

    rainfall: RainfallSeries
    parameters: ScsCnParameters


def _clamp_curve_number(curve_number: float) -> float:
    # CN II = 100 maps to exactly 100; binary division can land one ulp above.
    if 100.0 < curve_number <= 100.0 + 1e-9:
        return 100.0
    if not math.isfinite(curve_number) or not 0.0 < curve_number <= 100.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_CURVE_NUMBER,
            f"effective curve number must be in (0, 100], got {curve_number!r}",
        )
    return curve_number


def _curve_number_for_amc(cn_ii: float, amc: AmcCondition) -> float:
    if cn_ii == 100.0:
        return 100.0
    if amc is AmcCondition.AMC_II:
        return cn_ii
    if amc is AmcCondition.AMC_I:
        converted = (4.2 * cn_ii) / (10.0 - 0.058 * cn_ii)
    else:
        converted = (23.0 * cn_ii) / (10.0 + 0.13 * cn_ii)
    return _clamp_curve_number(converted)


def _composite_curve_number(curve_number: float, impervious_fraction: float) -> float:
    converted = curve_number + impervious_fraction * (
        _COMPOSITE_IMPERVIOUS_CN - curve_number
    )
    return _clamp_curve_number(converted)


def _potential_retention_mm(curve_number: float) -> float:
    if curve_number == 100.0:
        return 0.0
    retention = _S_NUMERATOR_MM / curve_number - _S_OFFSET_MM
    if not math.isfinite(retention) or retention < 0.0:
        raise ScsCnError(
            ScsCnErrorCode.INVALID_CURVE_NUMBER,
            f"curve number {curve_number} does not yield a finite potential retention",
        )
    return retention


def _resolve(parameters: ScsCnParameters) -> tuple[float, float, float, float]:
    """Return the SCS curve number, S, Ia, and the split impervious weight.

    Composite adjustment is applied after the antecedent-moisture conversion
    and replaces the split weight: the impervious fraction is already inside CN.
    """

    curve_number = _curve_number_for_amc(
        parameters.cn_ii,
        _require_amc(parameters.amc_condition),
    )
    split_fraction = parameters.impervious_fraction
    if _require_treatment(parameters.impervious_treatment) is ImperviousTreatment.COMPOSITE_CN:
        curve_number = _composite_curve_number(
            curve_number,
            parameters.impervious_fraction,
        )
        split_fraction = 0.0
    retention = _potential_retention_mm(curve_number)
    return (
        curve_number,
        retention,
        parameters.lambda_coeff * retention,
        split_fraction,
    )


def _cumulative_excess(
    rainfall_mm: float,
    initial_abstraction_mm: float,
    retention_mm: float,
) -> float:
    if rainfall_mm <= initial_abstraction_mm:
        return 0.0
    if retention_mm == 0.0:
        return rainfall_mm
    above = rainfall_mm - initial_abstraction_mm
    return above * above / (above + retention_mm)


def _cumulative_infiltration(
    rainfall_mm: float,
    initial_abstraction_mm: float,
    retention_mm: float,
) -> float:
    if rainfall_mm <= initial_abstraction_mm or retention_mm == 0.0:
        return 0.0
    above = rainfall_mm - initial_abstraction_mm
    return retention_mm * above / (above + retention_mm)


def _watershed_totals(
    rainfall_mm: float,
    initial_abstraction_mm: float,
    retention_mm: float,
    impervious_fraction: float,
) -> tuple[float, float, float]:
    """Watershed-average cumulative excess, consumed Ia, and infiltration."""

    pervious = 1.0 - impervious_fraction
    excess = impervious_fraction * rainfall_mm + pervious * _cumulative_excess(
        rainfall_mm,
        initial_abstraction_mm,
        retention_mm,
    )
    consumed = pervious * min(rainfall_mm, initial_abstraction_mm)
    infiltration = pervious * _cumulative_infiltration(
        rainfall_mm,
        initial_abstraction_mm,
        retention_mm,
    )
    return excess, consumed, infiltration


def _incremental_excess(current: float, previous: float, rainfall_mm: float) -> float:
    increment = current - previous
    if increment < 0.0:
        if increment < -_BALANCE_TOLERANCE_MM:
            raise ScsCnError(
                ScsCnErrorCode.WATER_BALANCE,
                "cumulative excess decreased between timesteps",
            )
        return 0.0
    if increment > rainfall_mm:
        if increment - rainfall_mm > _BALANCE_TOLERANCE_MM:
            raise ScsCnError(
                ScsCnErrorCode.WATER_BALANCE,
                "incremental excess exceeds incremental rainfall",
            )
        return rainfall_mm
    return increment


def _normalize_rainfall(rainfall: RainfallSeries) -> list[tuple[str | None, float]]:
    if isinstance(rainfall, str | bytes) or not isinstance(rainfall, Sequence):
        raise ScsCnError(
            ScsCnErrorCode.INVALID_RAINFALL,
            "rainfall must be a sequence of depths or timed rainfall points",
        )
    if len(rainfall) == 0:
        raise ScsCnError(ScsCnErrorCode.EMPTY_RAINFALL, "rainfall series is empty")

    if isinstance(rainfall[0], TimedRainfall):
        return _normalize_timed(rainfall)
    return _normalize_depths(rainfall)


def _normalize_timed(rainfall: Sequence[object]) -> list[tuple[str | None, float]]:
    steps: list[tuple[str | None, float]] = []
    previous: datetime | None = None
    for index, point in enumerate(rainfall):
        if not isinstance(point, TimedRainfall):
            raise ScsCnError(
                ScsCnErrorCode.INVALID_RAINFALL,
                f"rainfall at index {index} must be TimedRainfall, matching the first point",
            )
        instant = _parse_rfc3339(point.time, role=f"rainfall time at index {index}")
        if previous is not None and instant <= previous:
            raise ScsCnError(
                ScsCnErrorCode.UNORDERED_TIMESTAMPS,
                f"rainfall time at index {index} must be strictly after the previous time",
            )
        previous = instant
        steps.append((point.time, point.depth_mm))
    return steps


def _normalize_depths(rainfall: Sequence[object]) -> list[tuple[str | None, float]]:
    steps: list[tuple[str | None, float]] = []
    for index, point in enumerate(rainfall):
        if isinstance(point, TimedRainfall):
            raise ScsCnError(
                ScsCnErrorCode.INVALID_RAINFALL,
                f"rainfall at index {index} is timed, but the series started with a bare depth",
            )
        steps.append((None, _require_depth(point, index)))
    return steps


def _check_balance(result: ScsCnResult) -> None:
    summed_excess = math.fsum(step.excess_mm for step in result.steps)
    summed_loss = math.fsum(step.loss_mm for step in result.steps)
    summed_rain = math.fsum(step.rainfall_mm for step in result.steps)
    residual_increments = abs(summed_excess - result.total_excess_mm)
    residual_partition = abs(summed_excess + summed_loss - summed_rain)
    residual_components = abs(
        result.total_excess_mm
        + result.initial_abstraction_consumed_mm
        + result.cumulative_infiltration_mm
        - result.total_rainfall_mm
    )
    if (
        residual_increments > _BALANCE_TOLERANCE_MM
        or residual_partition > _BALANCE_TOLERANCE_MM
        or residual_components > _BALANCE_TOLERANCE_MM
        or abs(summed_rain - result.total_rainfall_mm) > _BALANCE_TOLERANCE_MM
    ):
        raise ScsCnError(
            ScsCnErrorCode.WATER_BALANCE,
            "SCS-CN water balance residual exceeds "
            f"{_BALANCE_TOLERANCE_MM} mm "
            f"(increments={residual_increments}, partition={residual_partition}, "
            f"components={residual_components})",
        )


def compute_excess(rainfall: RainfallSeries, parameters: ScsCnParameters) -> ScsCnResult:
    """Convert one rainfall series into excess rainfall.

    Cumulative excess is evaluated at each timestep from the storm total so
    far, then differenced. The curve-number equation is not applied to the
    individual timestep depths.
    """

    points = _normalize_rainfall(rainfall)
    curve_number, retention, abstraction, split = _resolve(parameters)

    steps: list[ExcessStep] = []
    cumulative_rainfall = 0.0
    previous_excess = 0.0
    for time, depth in points:
        cumulative_rainfall += depth
        excess, consumed, infiltration = _watershed_totals(
            cumulative_rainfall,
            abstraction,
            retention,
            split,
        )
        incremental = _incremental_excess(excess, previous_excess, depth)
        steps.append(
            ExcessStep(
                time=time,
                rainfall_mm=depth,
                cumulative_rainfall_mm=cumulative_rainfall,
                excess_mm=incremental,
                cumulative_excess_mm=excess,
                loss_mm=depth - incremental,
                cumulative_initial_abstraction_mm=consumed,
                cumulative_infiltration_mm=infiltration,
            )
        )
        previous_excess = excess

    result_steps = tuple(steps)
    last = result_steps[-1]
    total_rainfall = last.cumulative_rainfall_mm
    total_excess = last.cumulative_excess_mm
    result = ScsCnResult(
        curve_number=curve_number,
        potential_retention_mm=retention,
        initial_abstraction_mm=abstraction,
        total_rainfall_mm=total_rainfall,
        total_excess_mm=total_excess,
        initial_abstraction_consumed_mm=last.cumulative_initial_abstraction_mm,
        cumulative_infiltration_mm=last.cumulative_infiltration_mm,
        runoff_coefficient=(
            None if total_rainfall == 0.0 else total_excess / total_rainfall
        ),
        steps=result_steps,
    )
    _check_balance(result)
    return result


def excess_after(
    rainfall: Sequence[TimedRainfall],
    parameters: ScsCnParameters,
    t_cut: str,
) -> RemainderExcess:
    """Incremental excess for timesteps strictly after ``t_cut``.

    ``t_cut`` is an RFC3339 instant, usually the scenario ``as_of`` or the last
    arrival from an upstream station. The storm is computed in full first.
    Initial abstraction already taken by rainfall at or before ``t_cut`` stays
    taken; the tail is not restarted as a new storm.
    """

    cut = _parse_rfc3339(t_cut, role="t_cut")
    storm = compute_excess(rainfall, parameters)
    selected: list[ExcessStep] = []
    for step in storm.steps:
        if step.time is None:
            raise ScsCnError(
                ScsCnErrorCode.INVALID_TIMESTAMP,
                "excess_after requires an RFC3339 rainfall series",
            )
        if _parse_rfc3339(step.time, role="rainfall time") > cut:
            selected.append(step)
    tail = tuple(selected)
    return RemainderExcess(
        t_cut=t_cut,
        remainder_rainfall_mm=math.fsum(step.rainfall_mm for step in tail),
        remainder_excess_mm=math.fsum(step.excess_mm for step in tail),
        steps=tail,
        storm=storm,
    )


def compute_excess_batch(cases: Sequence[ScsCnCase]) -> tuple[ScsCnResult, ...]:
    """Run each rainfall / parameter pair independently, in order."""

    if isinstance(cases, str | bytes) or not isinstance(cases, Sequence):
        raise ScsCnError(
            ScsCnErrorCode.INVALID_RAINFALL,
            "batch cases must be a sequence of ScsCnCase",
        )
    results: list[ScsCnResult] = []
    for index, case in enumerate(cases):
        if not isinstance(case, ScsCnCase):
            raise ScsCnError(
                ScsCnErrorCode.INVALID_RAINFALL,
                f"batch case at index {index} must be ScsCnCase",
            )
        results.append(compute_excess(case.rainfall, case.parameters))
    return tuple(results)

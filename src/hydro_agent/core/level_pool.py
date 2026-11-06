"""Level-pool (Modified Puls) reservoir routing and reverse inflow.

Stage is in metres, storage in m³, discharge in m³/s, and clock times are
RFC3339 strings. Nothing in this module reads the wall clock. A scenario
``as_of`` is an argument: the pool is integrated only at that instant and
later, from the stage anchored there. Earlier instants keep the observed
stage, outflow, and gate opening, including points on the previous calendar
day.

Outflow is either an ungated discharge curve or a gated weir/orifice rating.
Passing both, or neither, is a domain error. With ``as_of`` set, the result
covers at least ``[as_of - 12 h, as_of + 12 h]``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from enum import StrEnum

_G = 9.81
_SQRT_2G = math.sqrt(2.0 * _G)
_US_PER_HOUR = 3_600_000_000
_WINDOW_HALF_HOURS = 12.0
_MAX_STEPS = 1_000_000
_STAGE_TOL_M = 1e-7
_BALANCE_TOL_M3 = 1e-9
_LEVEL_TOL_M = 1e-6
_OPENING_ATOL_M = 1e-9

_RFC3339 = re.compile(
    r"^\d{4}-\d{2}-\d{2}[Tt]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:[Zz]|[+-]\d{2}:\d{2})$"
)
_EPOCH = datetime(1970, 1, 1, tzinfo=UTC)


class LevelPoolErrorCode(StrEnum):
    """Stable identifier for a level-pool domain failure."""

    INVALID_STAGE = "invalid_stage"
    INVALID_STORAGE = "invalid_storage"
    INVALID_DISCHARGE = "invalid_discharge"
    INVALID_TIMESTAMP = "invalid_timestamp"
    UNORDERED_TIMESTAMPS = "unordered_timestamps"
    EMPTY_SERIES = "empty_series"
    INVALID_TIMESTEP = "invalid_timestep"
    MISALIGNED_TIMES = "misaligned_times"
    NON_MONOTONE_STORAGE = "non_monotone_storage"
    NON_MONOTONE_DISCHARGE = "non_monotone_discharge"
    OUTFLOW_CONTROL_CONFLICT = "outflow_control_conflict"
    STAGE_OUT_OF_RANGE = "stage_out_of_range"
    STORAGE_OUT_OF_RANGE = "storage_out_of_range"
    INVALID_GATE = "invalid_gate"
    INVALID_OPENING = "invalid_opening"
    INVALID_TARGET = "invalid_target"
    ROOT_NOT_BRACKETED = "root_not_bracketed"
    INVALID_PARAMETERS = "invalid_parameters"


class LevelPoolError(Exception):
    """Illegal level-pool input, or a step the storage curve cannot close.

    ``code`` identifies the failure. ``message`` describes the offending value.
    """

    def __init__(self, code: LevelPoolErrorCode, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code.value}: {message}")


class LevelPoolMode(StrEnum):
    """How outflow is chosen at each instant at or after ``as_of``."""

    UNGATED = "ungated"
    SCHEDULE = "schedule"
    CONSTANT_TARGET = "constant_target"
    PRE_RELEASE = "pre_release"
    Q_SAFE = "q_safe"
    Z_MAX = "z_max"
    UNGATED_SPILL = "ungated_spill"
    GATE_FAULT = "gate_fault"


def _require_real(value: object, code: LevelPoolErrorCode, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise LevelPoolError(code, f"{name} must be a real number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise LevelPoolError(code, f"{name} must be finite, got {value!r}")
    return number


def _require_int(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_GATE,
            f"{name} must be an integer, got {value!r}",
        )
    return value


def _parse_rfc3339(value: object, *, role: str) -> datetime:
    if not isinstance(value, str) or _RFC3339.fullmatch(value) is None:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTAMP,
            f"{role} must be an RFC3339 timestamp, got {value!r}",
        )
    normalized = f"{value[:-1]}+00:00" if value[-1] in "Zz" else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTAMP,
            f"{role} is not a valid RFC3339 timestamp: {value}",
        ) from exc
    if parsed.tzinfo is None:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTAMP,
            f"{role} must include a timezone offset, got {value!r}",
        )
    return parsed


def _to_us(instant: datetime) -> int:
    utc = instant.astimezone(UTC)
    delta = utc - _EPOCH
    return delta.days * 86_400 * 1_000_000 + delta.seconds * 1_000_000 + delta.microseconds


def _format_us(microseconds: int) -> str:
    instant = _EPOCH + timedelta(microseconds=microseconds)
    if instant.microsecond:
        text = instant.strftime("%Y-%m-%dT%H:%M:%S.%f").rstrip("0")
        return f"{text}Z"
    return instant.strftime("%Y-%m-%dT%H:%M:%SZ")


def _step_us(timestep_hours: float) -> int:
    number = _require_real(
        timestep_hours,
        LevelPoolErrorCode.INVALID_TIMESTEP,
        "timestep_hours",
    )
    if number <= 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be > 0, got {timestep_hours!r}",
        )
    step = round(number * _US_PER_HOUR)
    if step <= 0 or abs(number * _US_PER_HOUR - step) > 1e-3:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be an exact microsecond count, got {timestep_hours!r}",
        )
    window = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    if window % step != 0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            "timestep_hours must divide the 12 h scenario window",
        )
    return step


def _parse_series(
    series: object,
    name: str,
    *,
    code: LevelPoolErrorCode,
    allow_empty: bool = False,
) -> list[tuple[int, float]]:
    if series is None:
        if allow_empty:
            return []
        raise LevelPoolError(code, f"{name} is missing")
    if isinstance(series, str | bytes) or not isinstance(series, Sequence):
        raise LevelPoolError(code, f"{name} must be a sequence of (time, value) pairs")
    if len(series) == 0:
        if allow_empty:
            return []
        raise LevelPoolError(LevelPoolErrorCode.EMPTY_SERIES, f"{name} is empty")
    if len(series) > _MAX_STEPS:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            f"{name} has {len(series)} points; the limit is {_MAX_STEPS}",
        )
    parsed: list[tuple[int, float]] = []
    previous: int | None = None
    for index, item in enumerate(series):
        if not isinstance(item, tuple) or len(item) != 2:
            raise LevelPoolError(code, f"{name} at index {index} must be (time, value)")
        stamp, raw = item
        instant = _to_us(_parse_rfc3339(stamp, role=f"{name} time at index {index}"))
        value = _require_real(raw, code, f"{name} at index {index}")
        if previous is not None and instant <= previous:
            raise LevelPoolError(
                LevelPoolErrorCode.UNORDERED_TIMESTAMPS,
                f"{name} time at index {index} must be strictly after the previous time",
            )
        previous = instant
        parsed.append((instant, value))
    return parsed


def _bracket(xs: Sequence[float], x: float) -> tuple[int, int]:
    lo = 0
    hi = len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    return lo, hi


def _interp(xs: Sequence[float], ys: Sequence[float], x: float) -> float:
    lo, hi = _bracket(xs, x)
    span = xs[hi] - xs[lo]
    if span == 0.0:
        return ys[lo]
    weight = (x - xs[lo]) / span
    return ys[lo] + weight * (ys[hi] - ys[lo])


def _sample(
    table: Sequence[tuple[int, float]],
    instant: int,
    *,
    outside: float,
) -> float:
    """Linear value inside the table; ``outside`` beyond either end."""

    if not table or instant < table[0][0] or instant > table[-1][0]:
        return outside
    lo = 0
    hi = len(table) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if table[mid][0] <= instant:
            lo = mid
        else:
            hi = mid
    left_t, left_v = table[lo]
    right_t, right_v = table[hi]
    if instant == left_t:
        return left_v
    if instant == right_t or right_t == left_t:
        return right_v
    weight = (instant - left_t) / (right_t - left_t)
    return left_v + weight * (right_v - left_v)


def _root(func: Callable[[float], float], left: float, right: float) -> float:
    """Secant/bisection root of a strictly monotone function on ``[left, right]``."""

    a = left
    b = right
    fa = func(a)
    fb = func(b)
    if abs(fa) <= _BALANCE_TOL_M3:
        return a
    if abs(fb) <= _BALANCE_TOL_M3:
        return b
    if fa * fb > 0.0:
        end = "below" if fa > 0.0 else "above"
        raise LevelPoolError(
            LevelPoolErrorCode.ROOT_NOT_BRACKETED,
            "storage curve does not bracket the Modified Puls step "
            f"(Φ({a:.6f} m)={fa:.6g} m³, Φ({b:.6f} m)={fb:.6g} m³); "
            f"the pool would finish {end} the tabulated range",
        )
    span = abs(b - a)
    best = a if abs(fa) <= abs(fb) else b
    best_f = fa if abs(fa) <= abs(fb) else fb
    for _ in range(80):
        if fb != fa:
            guess = b - fb * (b - a) / (fb - fa)
        else:
            guess = 0.5 * (a + b)
        outside = guess <= min(a, b) or guess >= max(a, b)
        if outside or abs(guess - b) >= 0.5 * span:
            guess = 0.5 * (a + b)
        fg = func(guess)
        if abs(fg) < abs(best_f):
            best = guess
            best_f = fg
        converged = abs(b - a) <= _STAGE_TOL_M and abs(fg) <= 1e-6
        if abs(fg) <= _BALANCE_TOL_M3 or converged:
            return guess
        span = abs(b - a)
        if fa * fg < 0.0:
            b, fb = guess, fg
        else:
            a, fa = guess, fg
        if abs(fa) < abs(fb):
            a, b, fa, fb = b, a, fb, fa
    if abs(best_f) > 1e-6:
        raise LevelPoolError(
            LevelPoolErrorCode.ROOT_NOT_BRACKETED,
            f"Modified Puls step did not converge (residual {best_f:.6g} m³)",
        )
    return best


@dataclass(frozen=True, slots=True)
class StorageCurve:
    """Strictly increasing stage–storage table.

    Stage is metres and storage is cubic metres. Both columns increase
    together, so either one determines the other inside the table.
    """

    stages_m: tuple[float, ...]
    volumes_m3: tuple[float, ...]

    def __post_init__(self) -> None:
        stages, volumes = _validated_storage(self.stages_m, self.volumes_m3)
        object.__setattr__(self, "stages_m", stages)
        object.__setattr__(self, "volumes_m3", volumes)


def _validated_storage(
    stages: object,
    volumes: object,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if isinstance(stages, str | bytes) or not isinstance(stages, Sequence):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_STORAGE,
            "stages_m must be a sequence of stages",
        )
    if isinstance(volumes, str | bytes) or not isinstance(volumes, Sequence):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_STORAGE,
            "volumes_m3 must be a sequence of storages",
        )
    if len(stages) != len(volumes) or len(stages) < 2:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_STORAGE,
            "a storage curve needs at least two (stage, volume) points of equal length",
        )
    parsed_z: list[float] = []
    parsed_v: list[float] = []
    for index, (raw_z, raw_v) in enumerate(zip(stages, volumes, strict=True)):
        stage = _require_real(raw_z, LevelPoolErrorCode.INVALID_STAGE, f"stage at index {index}")
        volume = _require_real(
            raw_v,
            LevelPoolErrorCode.INVALID_STORAGE,
            f"volume at index {index}",
        )
        if volume < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_STORAGE,
                f"volume at index {index} must be >= 0 m³, got {raw_v!r}",
            )
        if parsed_z and (stage <= parsed_z[-1] or volume <= parsed_v[-1]):
            raise LevelPoolError(
                LevelPoolErrorCode.NON_MONOTONE_STORAGE,
                "storage curve stage and volume must both be strictly increasing; "
                f"point {index} is ({stage}, {volume}) after ({parsed_z[-1]}, {parsed_v[-1]})",
            )
        parsed_z.append(stage)
        parsed_v.append(volume)
    return tuple(parsed_z), tuple(parsed_v)


@dataclass(frozen=True, slots=True)
class DischargeCurve:
    """Ungated stage–discharge table.

    Stage increases strictly and discharge is non-decreasing and non-negative.
    Discharge is 0 at or below the crest, the highest tabulated stage whose
    discharge is 0. Below the first stage the discharge is also 0.
    """

    stages_m: tuple[float, ...]
    discharges_m3s: tuple[float, ...]

    def __post_init__(self) -> None:
        stages, discharges = _validated_discharge(self.stages_m, self.discharges_m3s)
        object.__setattr__(self, "stages_m", stages)
        object.__setattr__(self, "discharges_m3s", discharges)

    @property
    def crest_stage_m(self) -> float | None:
        """Highest tabulated stage with zero discharge, if the table has one."""

        crest: float | None = None
        for stage, discharge in zip(self.stages_m, self.discharges_m3s, strict=True):
            if discharge <= 0.0:
                crest = stage
            else:
                break
        return crest


def _validated_discharge(
    stages: object,
    discharges: object,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if isinstance(stages, str | bytes) or not isinstance(stages, Sequence):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_DISCHARGE,
            "stages_m must be a sequence of stages",
        )
    if isinstance(discharges, str | bytes) or not isinstance(discharges, Sequence):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_DISCHARGE,
            "discharges_m3s must be a sequence of discharges",
        )
    if len(stages) != len(discharges) or len(stages) < 2:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_DISCHARGE,
            "a discharge curve needs at least two (stage, discharge) points",
        )
    parsed_z: list[float] = []
    parsed_q: list[float] = []
    for index, (raw_z, raw_q) in enumerate(zip(stages, discharges, strict=True)):
        stage = _require_real(raw_z, LevelPoolErrorCode.INVALID_STAGE, f"stage at index {index}")
        discharge = _require_real(
            raw_q,
            LevelPoolErrorCode.INVALID_DISCHARGE,
            f"discharge at index {index}",
        )
        if discharge < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"discharge at index {index} must be >= 0 m³/s, got {raw_q!r}",
            )
        if parsed_z and stage <= parsed_z[-1]:
            raise LevelPoolError(
                LevelPoolErrorCode.NON_MONOTONE_DISCHARGE,
                "discharge-curve stage must be strictly increasing; "
                f"point {index} stage {stage} is not above {parsed_z[-1]}",
            )
        if parsed_q and discharge < parsed_q[-1]:
            raise LevelPoolError(
                LevelPoolErrorCode.NON_MONOTONE_DISCHARGE,
                "discharge-curve discharge must be non-decreasing; "
                f"point {index} discharge {discharge} is below {parsed_q[-1]}",
            )
        parsed_z.append(stage)
        parsed_q.append(discharge)
    return tuple(parsed_z), tuple(parsed_q)


@dataclass(frozen=True, slots=True)
class GateGeometry:
    """Spillway bays rated as a weir, an orifice, or a blend of the two.

    ``n_gates`` bays of clear width ``width_m`` sit on ``sill_elevation_m``.
    ``max_opening_m`` is the full lift. ``weir_coefficient`` is m and
    ``orifice_coefficient`` is μ. ``contraction_coefficient`` is the side
    contraction ε. ``xi_orifice`` and ``xi_weir`` bound the opening/head
    ratio: orifice at or below the first, weir at or above the second, and a
    linear blend between them. A fully open bay uses the weir formula.
    ``orifice_head_factor`` is ε_v in the orifice head ``H0 - ε_v e``.
    ``n_fault`` bays stay at a stuck opening and do not take the command.
    """

    n_gates: int
    width_m: float
    sill_elevation_m: float
    max_opening_m: float
    weir_coefficient: float
    orifice_coefficient: float
    contraction_coefficient: float = 1.0
    xi_orifice: float = 0.60
    xi_weir: float = 0.65
    orifice_head_factor: float = 0.0
    n_fault: int = 0

    def __post_init__(self) -> None:
        n_gates = _require_int(self.n_gates, "n_gates")
        if n_gates < 1:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"n_gates must be >= 1, got {self.n_gates!r}",
            )
        width = _require_real(self.width_m, LevelPoolErrorCode.INVALID_GATE, "width_m")
        if width <= 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"width_m must be > 0, got {self.width_m!r}",
            )
        sill = _require_real(
            self.sill_elevation_m,
            LevelPoolErrorCode.INVALID_GATE,
            "sill_elevation_m",
        )
        opening = _require_real(
            self.max_opening_m,
            LevelPoolErrorCode.INVALID_GATE,
            "max_opening_m",
        )
        if opening <= 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"max_opening_m must be > 0, got {self.max_opening_m!r}",
            )
        weir = _require_real(
            self.weir_coefficient,
            LevelPoolErrorCode.INVALID_GATE,
            "weir_coefficient",
        )
        orifice = _require_real(
            self.orifice_coefficient,
            LevelPoolErrorCode.INVALID_GATE,
            "orifice_coefficient",
        )
        if weir <= 0.0 or orifice <= 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                "weir_coefficient and orifice_coefficient must be > 0",
            )
        contraction = _require_real(
            self.contraction_coefficient,
            LevelPoolErrorCode.INVALID_GATE,
            "contraction_coefficient",
        )
        if contraction <= 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"contraction_coefficient must be > 0, got {self.contraction_coefficient!r}",
            )
        xi_orifice = _require_real(
            self.xi_orifice,
            LevelPoolErrorCode.INVALID_GATE,
            "xi_orifice",
        )
        xi_weir = _require_real(self.xi_weir, LevelPoolErrorCode.INVALID_GATE, "xi_weir")
        if xi_orifice < 0.0 or xi_weir < xi_orifice:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                "require 0 <= xi_orifice <= xi_weir, "
                f"got xi_orifice={self.xi_orifice!r}, xi_weir={self.xi_weir!r}",
            )
        head_factor = _require_real(
            self.orifice_head_factor,
            LevelPoolErrorCode.INVALID_GATE,
            "orifice_head_factor",
        )
        if head_factor < 0.0 or head_factor > 0.5:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"orifice_head_factor must be in [0, 0.5], got {self.orifice_head_factor!r}",
            )
        n_fault = _require_int(self.n_fault, "n_fault")
        if n_fault < 0 or n_fault > n_gates:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_GATE,
                f"n_fault must be in [0, n_gates], got {self.n_fault!r}",
            )
        object.__setattr__(self, "n_gates", n_gates)
        object.__setattr__(self, "width_m", width)
        object.__setattr__(self, "sill_elevation_m", sill)
        object.__setattr__(self, "max_opening_m", opening)
        object.__setattr__(self, "weir_coefficient", weir)
        object.__setattr__(self, "orifice_coefficient", orifice)
        object.__setattr__(self, "contraction_coefficient", contraction)
        object.__setattr__(self, "xi_orifice", xi_orifice)
        object.__setattr__(self, "xi_weir", xi_weir)
        object.__setattr__(self, "orifice_head_factor", head_factor)
        object.__setattr__(self, "n_fault", n_fault)


@dataclass(frozen=True, slots=True)
class GateCommand:
    """Opening that realises a target discharge, and the discharge it can pass."""

    opening_m: float
    discharge_m3s: float
    limited_by_capacity: bool


@dataclass(frozen=True, slots=True)
class LevelPoolPoint:
    """One instant of a level-pool trace.

    ``gate_opening_m`` is ``None`` for an ungated spillway. ``hours_from_as_of``
    is measured from the scenario instant, negative before it.
    """

    time: str
    hours_from_as_of: float
    stage_m: float
    storage_m3: float
    inflow_m3s: float
    outflow_m3s: float
    gate_opening_m: float | None


@dataclass(frozen=True, slots=True)
class LevelPoolResult:
    """Integrated pool trace and the flood-control summary at or after ``as_of``."""

    points: tuple[LevelPoolPoint, ...]
    as_of: str
    timestep_hours: float
    z_max_m: float
    t_z_max: str
    v_max_m3: float
    q_out_max_m3s: float
    t_q_out_max: str
    q_in_max_m3s: float
    t_q_in_max: str
    peak_reduction_m3s: float
    peak_reduction_pct: float
    storage_utilization_pct: float | None
    z_end_m: float
    time_reached_target: str | None


@dataclass(frozen=True, slots=True)
class ReservoirInflowPoint:
    """One inflow ordinate, with hours measured from the scenario ``as_of``."""

    time: str
    hours_from_as_of: float
    discharge_m3s: float


@dataclass(frozen=True, slots=True)
class ReservoirInflow:
    """Inflow spliced from observations, reverse routing, and a forecast."""

    points: tuple[ReservoirInflowPoint, ...]
    as_of: str
    t_fail: str
    timestep_hours: float

    def as_pairs(self) -> tuple[tuple[str, float], ...]:
        """``(time, discharge)`` pairs accepted by the routing functions."""

        return tuple((point.time, point.discharge_m3s) for point in self.points)


def storage_volume(curve: StorageCurve, stage_m: float) -> float:
    """Storage (m³) on a strictly monotone stage–storage curve."""

    _require_curve(curve)
    stage = _require_real(stage_m, LevelPoolErrorCode.INVALID_STAGE, "stage_m")
    _require_stage_inside(curve, stage)
    return _interp(curve.stages_m, curve.volumes_m3, stage)


def stage_for_storage(curve: StorageCurve, volume_m3: float) -> float:
    """Stage (m) on a strictly monotone stage–storage curve."""

    _require_curve(curve)
    volume = _require_real(volume_m3, LevelPoolErrorCode.INVALID_STORAGE, "volume_m3")
    low = curve.volumes_m3[0]
    high = curve.volumes_m3[-1]
    if volume < low or volume > high:
        raise LevelPoolError(
            LevelPoolErrorCode.STORAGE_OUT_OF_RANGE,
            f"volume {volume} m³ is outside [{low}, {high}]",
        )
    return _interp(curve.volumes_m3, curve.stages_m, volume)


def ungated_discharge(curve: DischargeCurve, stage_m: float) -> float:
    """Free spillway discharge (m³/s). Zero at or below the crest."""

    _require_discharge_curve(curve)
    stage = _require_real(stage_m, LevelPoolErrorCode.INVALID_STAGE, "stage_m")
    crest = curve.crest_stage_m
    if stage < curve.stages_m[0] or (crest is not None and stage <= crest):
        return 0.0
    if stage >= curve.stages_m[-1]:
        z0 = curve.stages_m[-2]
        z1 = curve.stages_m[-1]
        q0 = curve.discharges_m3s[-2]
        q1 = curve.discharges_m3s[-1]
        slope = (q1 - q0) / (z1 - z0)
        return max(0.0, q1 + slope * (stage - z1))
    return max(0.0, _interp(curve.stages_m, curve.discharges_m3s, stage))


def _bay_discharge(gates: GateGeometry, stage_m: float, opening_m: float, n_bays: int) -> float:
    if n_bays <= 0 or opening_m <= _OPENING_ATOL_M:
        return 0.0
    head = stage_m - gates.sill_elevation_m
    if head <= 0.0:
        return 0.0
    opening = min(opening_m, gates.max_opening_m)
    q_weir = (
        gates.weir_coefficient
        * gates.contraction_coefficient
        * n_bays
        * gates.width_m
        * _SQRT_2G
        * head
        * math.sqrt(head)
    )
    orifice_head = head - gates.orifice_head_factor * opening
    if orifice_head <= 0.0:
        q_orifice = 0.0
    else:
        q_orifice = (
            gates.orifice_coefficient
            * n_bays
            * gates.width_m
            * opening
            * math.sqrt(2.0 * _G * orifice_head)
        )
    q_orifice = min(q_orifice, q_weir)
    if opening >= gates.max_opening_m - _OPENING_ATOL_M:
        return q_weir
    ratio = opening / head
    if ratio >= gates.xi_weir:
        return q_weir
    if ratio <= gates.xi_orifice or gates.xi_weir <= gates.xi_orifice:
        return q_orifice
    theta = (ratio - gates.xi_orifice) / (gates.xi_weir - gates.xi_orifice)
    return (1.0 - theta) * q_orifice + theta * q_weir


def gate_discharge(
    gates: GateGeometry,
    stage_m: float,
    opening_m: float,
    *,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
) -> float:
    """Discharge (m³/s) through healthy bays at ``opening_m`` plus stuck bays.

    Healthy bays are ``n_gates - n_fault``. Stuck bays pass ``stuck_opening_m``
    (closed when that opening is 0). A fully open healthy bay uses the weir.
    """

    _require_gates(gates)
    stage = _require_real(stage_m, LevelPoolErrorCode.INVALID_STAGE, "stage_m")
    opening = _require_real(opening_m, LevelPoolErrorCode.INVALID_OPENING, "opening_m")
    if opening < 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_OPENING,
            f"opening_m must be >= 0, got {opening_m!r}",
        )
    fault = _fault_count(gates, n_fault)
    stuck = _require_real(
        stuck_opening_m,
        LevelPoolErrorCode.INVALID_OPENING,
        "stuck_opening_m",
    )
    if stuck < 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_OPENING,
            f"stuck_opening_m must be >= 0, got {stuck_opening_m!r}",
        )
    healthy = gates.n_gates - fault
    return _bay_discharge(gates, stage, opening, healthy) + _bay_discharge(
        gates,
        stage,
        min(stuck, gates.max_opening_m),
        fault,
    )


def gate_opening(
    gates: GateGeometry,
    stage_m: float,
    target_m3s: float,
    *,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
) -> GateCommand:
    """Opening whose discharge equals ``target_m3s``, or full lift at capacity.

    When the target is at least the fully open discharge, the command is
    ``max_opening_m`` and the discharge is that capacity. The search is exact
    where discharge increases with opening.
    """

    _require_gates(gates)
    stage = _require_real(stage_m, LevelPoolErrorCode.INVALID_STAGE, "stage_m")
    target = _require_real(target_m3s, LevelPoolErrorCode.INVALID_TARGET, "target_m3s")
    if target < 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            f"target_m3s must be >= 0, got {target_m3s!r}",
        )
    fault = _fault_count(gates, n_fault)
    stuck = _require_real(
        stuck_opening_m,
        LevelPoolErrorCode.INVALID_OPENING,
        "stuck_opening_m",
    )
    if stuck < 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_OPENING,
            f"stuck_opening_m must be >= 0, got {stuck_opening_m!r}",
        )
    healthy = gates.n_gates - fault
    stuck_q = _bay_discharge(gates, stage, min(stuck, gates.max_opening_m), fault)

    def total(opening: float) -> float:
        return _bay_discharge(gates, stage, opening, healthy) + stuck_q

    q_max = total(gates.max_opening_m)
    if target <= stuck_q + 1e-9:
        return GateCommand(
            opening_m=0.0,
            discharge_m3s=stuck_q,
            limited_by_capacity=target > stuck_q + 1e-9,
        )
    if target >= q_max - 1e-8:
        return GateCommand(
            opening_m=gates.max_opening_m,
            discharge_m3s=q_max,
            limited_by_capacity=True,
        )
    lo = 0.0
    hi = gates.max_opening_m
    best = hi
    best_q = q_max
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        q_mid = total(mid)
        if abs(q_mid - target) < abs(best_q - target):
            best = mid
            best_q = q_mid
        if abs(q_mid - target) <= 1e-8 or hi - lo <= 1e-10:
            return GateCommand(opening_m=mid, discharge_m3s=q_mid, limited_by_capacity=False)
        if q_mid < target:
            lo = mid
        else:
            hi = mid
    return GateCommand(opening_m=best, discharge_m3s=best_q, limited_by_capacity=False)


def _require_curve(curve: object) -> None:
    if not isinstance(curve, StorageCurve):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_STORAGE,
            f"storage must be StorageCurve, got {curve!r}",
        )


def _require_discharge_curve(curve: object) -> None:
    if not isinstance(curve, DischargeCurve):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_DISCHARGE,
            f"discharge must be DischargeCurve, got {curve!r}",
        )


def _require_gates(gates: object) -> None:
    if not isinstance(gates, GateGeometry):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_GATE,
            f"gates must be GateGeometry, got {gates!r}",
        )


def _require_stage_inside(curve: StorageCurve, stage: float, *, name: str = "stage") -> None:
    low = curve.stages_m[0]
    high = curve.stages_m[-1]
    if stage < low or stage > high:
        raise LevelPoolError(
            LevelPoolErrorCode.STAGE_OUT_OF_RANGE,
            f"{name} {stage} m is outside [{low}, {high}]",
        )


def _fault_count(gates: GateGeometry, n_fault: int | None) -> int:
    fault = gates.n_fault if n_fault is None else _require_int(n_fault, "n_fault")
    if fault < 0 or fault > gates.n_gates:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_GATE,
            f"n_fault must be in [0, n_gates={gates.n_gates}], got {fault}",
        )
    return fault


def _require_mode_controls(
    *,
    mode: LevelPoolMode,
    discharge: DischargeCurve | None,
    gates: GateGeometry | None,
) -> None:
    has_discharge = discharge is not None
    has_gates = gates is not None
    if has_discharge == has_gates:
        raise LevelPoolError(
            LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT,
            "provide exactly one outflow control: a discharge curve or gate geometry",
        )
    if mode is LevelPoolMode.UNGATED and not has_discharge:
        raise LevelPoolError(
            LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT,
            "ungated routing uses a discharge curve, not gate geometry",
        )
    if mode is not LevelPoolMode.UNGATED and not has_gates:
        raise LevelPoolError(
            LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT,
            f"{mode.value} routing uses gate geometry, not a discharge curve",
        )


@dataclass(frozen=True, slots=True)
class _Job:
    storage: StorageCurve
    inflow: tuple[tuple[int, float], ...]
    mode: LevelPoolMode
    z0_m: float
    as_of: str
    as_of_us: int
    step_us: int
    timestep_hours: float
    discharge: DischargeCurve | None
    gates: GateGeometry | None
    n_fault: int
    stuck_opening_m: float
    q_safe_m3s: float | None
    q_target_m3s: float | None
    z_target_m: float | None
    z_hold_m: float | None
    z_lim_m: float | None
    wait_z_m: float | None
    tail_inflow_m3s: float
    schedule: tuple[tuple[int, float], ...] | None
    observed_stage: tuple[tuple[int, float], ...]
    observed_outflow: tuple[tuple[int, float], ...]
    observed_opening: tuple[tuple[int, float], ...]
    z_flood_limit_m: float | None
    z_design_m: float | None


def _positive_discharge(value: object, name: str) -> float:
    number = _require_real(value, LevelPoolErrorCode.INVALID_DISCHARGE, name)
    if number < 0.0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_DISCHARGE,
            f"{name} must be >= 0 m³/s, got {value!r}",
        )
    return number


def _optional_stage(value: object, name: str) -> float | None:
    if value is None:
        return None
    return _require_real(value, LevelPoolErrorCode.INVALID_STAGE, name)


def _build_job(
    storage: object,
    inflow: object,
    *,
    mode: LevelPoolMode,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    discharge: DischargeCurve | None,
    gates: GateGeometry | None,
    n_fault: int | None,
    stuck_opening_m: float,
    q_safe_m3s: float | None,
    q_target_m3s: float | None,
    z_target_m: float | None,
    z_hold_m: float | None,
    z_lim_m: float | None,
    wait_closed_until_z: float | None,
    tail_inflow_m3s: float,
    schedule: object,
    observed_stage: object,
    observed_outflow: object,
    observed_gate_opening: object,
    z_flood_limit_m: float | None,
    z_design_m: float | None,
) -> _Job:
    _require_curve(storage)
    assert isinstance(storage, StorageCurve)
    if not isinstance(mode, LevelPoolMode):
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_PARAMETERS,
            f"mode must be LevelPoolMode, got {mode!r}",
        )
    _require_mode_controls(mode=mode, discharge=discharge, gates=gates)
    if gates is not None:
        _require_gates(gates)
    if discharge is not None:
        _require_discharge_curve(discharge)
    stage0 = _require_real(z0_m, LevelPoolErrorCode.INVALID_STAGE, "z0_m")
    _require_stage_inside(storage, stage0)
    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of"))
    step = _step_us(timestep_hours)
    parsed_inflow = _parse_series(inflow, "inflow", code=LevelPoolErrorCode.INVALID_DISCHARGE)
    for index, (_instant, discharge_value) in enumerate(parsed_inflow):
        if discharge_value < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"inflow at index {index} must be >= 0 m³/s, got {discharge_value!r}",
            )
    stages = tuple(
        _parse_series(
            observed_stage,
            "observed_stage",
            code=LevelPoolErrorCode.INVALID_STAGE,
            allow_empty=True,
        )
    )
    outflows = tuple(
        _parse_series(
            observed_outflow,
            "observed_outflow",
            code=LevelPoolErrorCode.INVALID_DISCHARGE,
            allow_empty=True,
        )
    )
    openings = tuple(
        _parse_series(
            observed_gate_opening,
            "observed_gate_opening",
            code=LevelPoolErrorCode.INVALID_OPENING,
            allow_empty=True,
        )
    )
    for index, (_instant, outflow) in enumerate(outflows):
        if outflow < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"observed_outflow at index {index} must be >= 0 m³/s",
            )
    for instant, observed in stages:
        if instant > as_of_us:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TIMESTAMP,
                "observed_stage contains an instant after as_of",
            )
        _require_stage_inside(storage, observed)
    for instant, _outflow in outflows:
        if instant > as_of_us:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TIMESTAMP,
                "observed_outflow contains an instant after as_of",
            )
    for instant, opening in openings:
        if instant > as_of_us:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TIMESTAMP,
                "observed_gate_opening contains an instant after as_of",
            )
        if opening < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_OPENING,
                f"observed_gate_opening at {_format_us(instant)} must be >= 0",
            )
    fault = 0
    stuck = 0.0
    if gates is not None:
        fault = _fault_count(gates, n_fault)
        stuck = _require_real(
            stuck_opening_m,
            LevelPoolErrorCode.INVALID_OPENING,
            "stuck_opening_m",
        )
        if stuck < 0.0 or stuck > gates.max_opening_m:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_OPENING,
                "stuck_opening_m must be in [0, max_opening_m], "
                f"got {stuck_opening_m!r}",
            )
    q_safe = None if q_safe_m3s is None else _positive_discharge(q_safe_m3s, "q_safe_m3s")
    q_target = (
        None if q_target_m3s is None else _positive_discharge(q_target_m3s, "q_target_m3s")
    )
    z_target = _optional_stage(z_target_m, "z_target_m")
    z_hold = _optional_stage(z_hold_m, "z_hold_m")
    z_lim = _optional_stage(z_lim_m, "z_lim_m")
    wait_z = _optional_stage(wait_closed_until_z, "wait_closed_until_z")
    tail = _positive_discharge(tail_inflow_m3s, "tail_inflow_m3s")
    parsed_schedule: tuple[tuple[int, float], ...] | None = None
    if schedule is not None:
        rows = _parse_series(schedule, "schedule", code=LevelPoolErrorCode.INVALID_OPENING)
        if gates is None:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_PARAMETERS,
                "a gate schedule requires gate geometry",
            )
        checked: list[tuple[int, float]] = []
        for index, (instant, opening) in enumerate(rows):
            if opening < 0.0 or opening > gates.max_opening_m:
                raise LevelPoolError(
                    LevelPoolErrorCode.INVALID_OPENING,
                    f"schedule opening at index {index} must be in [0, max_opening_m]",
                )
            checked.append((instant, opening))
        parsed_schedule = tuple(checked)
    for name, stage in (
        ("z_target_m", z_target),
        ("z_hold_m", z_hold),
        ("z_lim_m", z_lim),
        ("wait_closed_until_z", wait_z),
        ("z_flood_limit_m", z_flood_limit_m),
        ("z_design_m", z_design_m),
    ):
        if stage is not None:
            _require_stage_inside(storage, stage, name=name)
    _check_mode_arguments(
        mode=mode,
        q_safe=q_safe,
        q_target=q_target,
        z_lim=z_lim,
        schedule=parsed_schedule,
        n_fault=fault,
        gates=gates,
        z0=stage0,
    )
    return _Job(
        storage=storage,
        inflow=tuple(parsed_inflow),
        mode=mode,
        z0_m=stage0,
        as_of=as_of,
        as_of_us=as_of_us,
        step_us=step,
        timestep_hours=float(timestep_hours),
        discharge=discharge,
        gates=gates,
        n_fault=fault,
        stuck_opening_m=stuck,
        q_safe_m3s=q_safe,
        q_target_m3s=q_target,
        z_target_m=z_target,
        z_hold_m=z_hold,
        z_lim_m=z_lim,
        wait_z_m=wait_z,
        tail_inflow_m3s=tail,
        schedule=parsed_schedule,
        observed_stage=stages,
        observed_outflow=outflows,
        observed_opening=openings,
        z_flood_limit_m=_optional_stage(z_flood_limit_m, "z_flood_limit_m"),
        z_design_m=_optional_stage(z_design_m, "z_design_m"),
    )


def _check_mode_arguments(
    *,
    mode: LevelPoolMode,
    q_safe: float | None,
    q_target: float | None,
    z_lim: float | None,
    schedule: tuple[tuple[int, float], ...] | None,
    n_fault: int,
    gates: GateGeometry | None,
    z0: float,
) -> None:
    if mode is LevelPoolMode.SCHEDULE and not schedule:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_PARAMETERS,
            "schedule routing requires a non-empty opening schedule",
        )
    if mode is LevelPoolMode.CONSTANT_TARGET and q_target is None:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            "constant-target routing requires q_target_m3s",
        )
    if mode in (LevelPoolMode.PRE_RELEASE, LevelPoolMode.Q_SAFE) and q_safe is None:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            f"{mode.value} routing requires q_safe_m3s",
        )
    if mode is LevelPoolMode.Z_MAX and z_lim is None:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            "route_to_z_max requires z_lim_m",
        )
    if mode is LevelPoolMode.Z_MAX and z_lim is not None and z_lim < z0:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            f"z_lim_m ({z_lim}) must be >= z0_m ({z0})",
        )
    if mode is LevelPoolMode.GATE_FAULT and n_fault < 1:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_GATE,
            "route_gate_fault requires n_fault >= 1",
        )
    if mode is not LevelPoolMode.UNGATED and gates is None:
        raise LevelPoolError(
            LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT,
            f"{mode.value} routing requires gate geometry",
        )


@dataclass
class _Pool:
    job: _Job
    dt_s: float
    z_target: float
    z_hold: float
    phase: str
    reached_us: int | None

    def volume(self, stage: float) -> float:
        return _interp(self.job.storage.stages_m, self.job.storage.volumes_m3, stage)

    def stage(self, volume: float) -> float:
        return stage_for_storage(self.job.storage, volume)

    def q_max(self, stage: float) -> float:
        gates = self.job.gates
        if gates is None:
            return 0.0
        return gate_discharge(
            gates,
            stage,
            gates.max_opening_m,
            n_fault=self.job.n_fault,
            stuck_opening_m=self.job.stuck_opening_m,
        )

    def end_volume(self, v1: float, i1: float, i2: float, o1: float, o2: float) -> float:
        return v1 + 0.5 * (i1 + i2) * self.dt_s - 0.5 * (o1 + o2) * self.dt_s

    def outflow_for(self, v1: float, v2: float, i1: float, i2: float, o1: float) -> float:
        return (i1 + i2) - o1 - 2.0 * (v2 - v1) / self.dt_s

    def solve_outflow(
        self,
        z1: float,
        o1: float,
        i1: float,
        i2: float,
        outflow_of: Callable[[float], float],
    ) -> tuple[float, float]:
        v1 = self.volume(z1)
        low = self.job.storage.stages_m[0]
        high = self.job.storage.stages_m[-1]

        def phi(stage: float) -> float:
            outflow = outflow_of(stage)
            return (
                self.volume(stage)
                - v1
                - 0.5 * (i1 + i2) * self.dt_s
                + 0.5 * (o1 + outflow) * self.dt_s
            )

        stage = _root(phi, low, high)
        _require_stage_inside(self.job.storage, stage)
        outflow = float(outflow_of(stage))
        return stage, outflow

    def release(
        self,
        z1: float,
        o1: float,
        i1: float,
        i2: float,
        i3: float,
        q_cap: float,
        z_floor: float,
    ) -> tuple[float, float, str]:
        """Release up to ``q_cap`` without finishing below ``z_floor``."""

        def capped(stage: float) -> float:
            return min(q_cap, self.q_max(stage))

        z_free, o_free = self.solve_outflow(z1, o1, i1, i2, capped)
        if z_free >= z_floor - 1e-8 and self._next_stays_up(z_free, o_free, i2, i3, z_floor):
            return z_free, o_free, "free"
        landed = self._land(z1, o1, i1, i2, i3, q_cap, z_floor)
        if landed is not None:
            return landed[0], landed[1], "land"
        magic = self._magic(z1, o1, i1, i2, i3, q_cap, z_floor)
        if magic is not None:
            return magic[0], magic[1], "magic"
        if z_free >= z_floor - 1e-8:
            return z_free, max(0.0, o_free), "free"
        return z_free, max(0.0, o_free), "below"

    def _next_stays_up(self, z2: float, o2: float, i2: float, i3: float, z_floor: float) -> bool:
        v3 = self.end_volume(self.volume(z2), i2, i3, o2, 0.0)
        return v3 >= self.volume(z_floor) - 1e-4

    def _land(
        self,
        z1: float,
        o1: float,
        i1: float,
        i2: float,
        i3: float,
        q_cap: float,
        z_floor: float,
    ) -> tuple[float, float] | None:
        o_land = self.outflow_for(self.volume(z1), self.volume(z_floor), i1, i2, o1)
        cap = min(q_cap, self.q_max(z_floor))
        if o_land > cap + 1e-8 or o_land < -1e-8:
            return None
        if o_land > i2 + i3 + 1e-6:
            return None
        if o_land < 0.0:
            v2 = self.end_volume(self.volume(z1), i1, i2, o1, 0.0)
            return self.stage(v2), 0.0
        return z_floor, o_land

    def _magic(
        self,
        z1: float,
        o1: float,
        i1: float,
        i2: float,
        i3: float,
        q_cap: float,
        z_floor: float,
    ) -> tuple[float, float] | None:
        # End this step so the next one can reach the floor with a closed gate.
        o_magic = (
            (self.volume(z1) - self.volume(z_floor)) / self.dt_s
            + 0.5 * (i1 + i2)
            - 0.5 * o1
            + 0.5 * (i2 + i3)
        )
        if o_magic < -1e-9 or o_magic > q_cap + 1e-9:
            return None
        o_use = min(q_cap, max(0.0, o_magic))
        v2 = self.end_volume(self.volume(z1), i1, i2, o1, o_use)
        low = self.job.storage.volumes_m3[0]
        high = self.job.storage.volumes_m3[-1]
        if v2 < low or v2 > high:
            return None
        z2 = self.stage(v2)
        if z2 < z_floor - 1e-8 or o_use > self.q_max(z2) + 1e-6:
            return None
        return z2, o_use

    def hold_level(
        self,
        z1: float,
        o1: float,
        i1: float,
        i2: float,
        i3: float,
        z_level: float,
        q_cap: float,
    ) -> tuple[float, float, str]:
        o_land = self.outflow_for(self.volume(z1), self.volume(z_level), i1, i2, o1)
        cap = min(q_cap, self.q_max(z_level))
        if -1e-8 <= o_land <= cap + 1e-8 and o_land <= i2 + i3 + 1e-6:
            if o_land < 0.0:
                v2 = self.end_volume(self.volume(z1), i1, i2, o1, 0.0)
                return self.stage(v2), 0.0, "hold"
            return z_level, o_land, "hold"
        if o_land > cap:
            return (*self.release(z1, o1, i1, i2, i3, q_cap, z_level)[:2], "rise")
        return (*self.release(z1, o1, i1, i2, i3, q_cap, z_level)[:2], "ease")


def _initial_phase(job: _Job) -> str:
    if job.mode is LevelPoolMode.UNGATED:
        return "ungated"
    if job.mode is LevelPoolMode.SCHEDULE:
        return "schedule"
    if job.mode is LevelPoolMode.CONSTANT_TARGET:
        return "constant"
    if job.mode is LevelPoolMode.PRE_RELEASE:
        return "draw"
    if job.mode is LevelPoolMode.Q_SAFE:
        return "wait" if job.wait_z_m is not None else "regulate"
    if job.mode is LevelPoolMode.Z_MAX:
        return "store"
    if job.mode is LevelPoolMode.GATE_FAULT and job.q_safe_m3s is not None:
        return "constant"
    return "spill"


def _scheduled_opening(job: _Job, instant: int) -> float:
    if not job.schedule:
        return 0.0
    opening = 0.0
    for time, value in job.schedule:
        if time <= instant:
            opening = value
        else:
            break
    return opening


def _targets(pool: _Pool) -> None:
    job = pool.job
    gates = job.gates
    if job.mode is not LevelPoolMode.PRE_RELEASE or gates is None:
        return
    sill = gates.sill_elevation_m
    target = sill if job.z_target_m is None else job.z_target_m
    if target < sill:
        target = sill
    _require_stage_inside(job.storage, target)
    if target >= job.z0_m:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            f"pre-release z_target_m ({target}) must be below z0_m ({job.z0_m}); "
            "a target below the sill is raised to the sill",
        )
    hold = job.z0_m if job.z_hold_m is None else job.z_hold_m
    _require_stage_inside(job.storage, hold)
    if hold < target:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TARGET,
            f"z_hold_m ({hold}) must be >= the drawdown target ({target})",
        )
    pool.z_target = target
    pool.z_hold = hold


def _step_pool(
    pool: _Pool,
    *,
    instant: int,
    z1: float,
    o1: float,
    i1: float,
    i2: float,
    i3: float,
) -> tuple[float, float]:
    job = pool.job
    q_safe = 0.0 if job.q_safe_m3s is None else job.q_safe_m3s
    q_target = 0.0 if job.q_target_m3s is None else job.q_target_m3s
    guard = 0
    while True:
        guard += 1
        if guard > 8:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_PARAMETERS,
                f"routing phase did not settle at {_format_us(instant)} ({pool.phase})",
            )
        phase = pool.phase
        if phase == "ungated":
            curve = job.discharge
            assert curve is not None

            def ungated(stage: float, rating: DischargeCurve = curve) -> float:
                return ungated_discharge(rating, stage)

            return pool.solve_outflow(z1, o1, i1, i2, ungated)
        if phase == "schedule":
            opening = _scheduled_opening(job, instant)
            gates = job.gates
            assert gates is not None

            def scheduled(
                stage: float,
                lift: float = opening,
                geometry: GateGeometry = gates,
            ) -> float:
                return gate_discharge(
                    geometry,
                    stage,
                    lift,
                    n_fault=job.n_fault,
                    stuck_opening_m=job.stuck_opening_m,
                )

            return pool.solve_outflow(z1, o1, i1, i2, scheduled)
        if phase in ("constant", "spill"):
            gates = job.gates
            assert gates is not None
            if phase == "spill":
                return pool.solve_outflow(z1, o1, i1, i2, pool.q_max)
            cap = q_safe if job.mode is LevelPoolMode.GATE_FAULT else q_target

            def capped_target(stage: float, limit: float = cap) -> float:
                return min(limit, pool.q_max(stage))

            return pool.solve_outflow(z1, o1, i1, i2, capped_target)
        if phase == "draw":
            if i2 > q_safe:
                pool.phase = "flood"
                continue
            z2, o2, how = pool.release(z1, o1, i1, i2, i3, q_safe, pool.z_target)
            if how == "land" or z2 <= pool.z_target + _LEVEL_TOL_M:
                pool.phase = "hold_target"
                if pool.reached_us is None:
                    pool.reached_us = instant
            return z2, o2
        if phase == "hold_target":
            o_land = pool.outflow_for(pool.volume(z1), pool.volume(pool.z_target), i1, i2, o1)
            if i2 > q_safe or o_land > max(pool.q_max(pool.z_target), 0.0) + 1e-6:
                pool.phase = "flood"
                continue
            z2, o2, _how = pool.hold_level(z1, o1, i1, i2, i3, pool.z_target, q_safe)
            if pool.reached_us is None and z2 <= pool.z_target + _LEVEL_TOL_M:
                pool.reached_us = instant
            return z2, o2
        if phase == "flood":
            if i2 <= q_safe and z1 >= pool.z_hold - _LEVEL_TOL_M:
                pool.phase = "recede"
                continue
            if i2 <= q_safe and z1 < pool.z_hold - _LEVEL_TOL_M:
                pool.phase = "refill"
                continue
            return pool.release(z1, o1, i1, i2, i3, q_safe, pool.z_target)[:2]
        if phase == "refill":
            cap = min(q_safe, pool.q_max(pool.z_hold))
            o_land = pool.outflow_for(pool.volume(z1), pool.volume(pool.z_hold), i1, i2, o1)
            v_shut = pool.end_volume(pool.volume(z1), i1, i2, o1, 0.0)
            z_shut = pool.stage(v_shut)
            if 0.0 <= o_land <= cap and o_land <= i2 + i3 + 1e-6 and z_shut >= pool.z_hold - 1e-8:
                pool.phase = "hold_flood"
                return pool.z_hold, o_land
            if z_shut <= pool.z_hold + 1e-8:
                return z_shut, 0.0
            pool.phase = "recede"
            continue
        if phase == "recede":
            z2, o2, how = pool.release(z1, o1, i1, i2, i3, q_safe, pool.z_hold)
            if how in ("land", "magic") or z2 <= pool.z_hold + _LEVEL_TOL_M:
                pool.phase = "hold_flood"
            return z2, o2
        if phase == "hold_flood":
            z2, o2, how = pool.hold_level(z1, o1, i1, i2, i3, pool.z_hold, q_safe)
            if how == "rise":
                pool.phase = "flood"
                continue
            return z2, o2
        if phase == "wait":
            assert job.wait_z_m is not None
            v_shut = pool.end_volume(pool.volume(z1), i1, i2, o1, 0.0)
            if pool.stage(v_shut) <= job.wait_z_m + 1e-8:
                return pool.stage(v_shut), 0.0
            pool.phase = "release"
            continue
        if phase == "release":
            assert job.wait_z_m is not None
            z2, o2, _how = pool.release(z1, o1, i1, i2, i3, q_safe, job.wait_z_m)
            if z2 <= job.wait_z_m + _LEVEL_TOL_M and i2 <= q_safe:
                pool.phase = "hold_wait"
            return z2, o2
        if phase == "hold_wait":
            assert job.wait_z_m is not None
            z2, o2, how = pool.hold_level(z1, o1, i1, i2, i3, job.wait_z_m, q_safe)
            if how == "rise":
                pool.phase = "release"
                continue
            return z2, o2
        if phase == "regulate":
            # Cap the release at q_safe and do not draw the pool below z_hold.
            z2, o2, how = pool.hold_level(z1, o1, i1, i2, i3, pool.z_hold, q_safe)
            if how == "rise":
                pool.phase = "regulate_release"
                continue
            return z2, o2
        if phase == "regulate_release":
            z2, o2, how = pool.release(z1, o1, i1, i2, i3, q_safe, pool.z_hold)
            if how in ("land", "magic") or (z2 <= pool.z_hold + _LEVEL_TOL_M and i2 <= q_safe):
                pool.phase = "regulate"
            return z2, o2
        if phase == "store":
            assert job.z_lim_m is not None
            v_shut = pool.end_volume(pool.volume(z1), i1, i2, o1, 0.0)
            if pool.stage(v_shut) <= job.z_lim_m + 1e-8:
                return pool.stage(v_shut), 0.0
            pool.phase = "pin"
            continue
        if phase == "pin":
            assert job.z_lim_m is not None
            if (
                job.q_safe_m3s is not None
                and job.z_hold_m is not None
                and i2 <= job.tail_inflow_m3s
            ):
                pool.phase = "draw_lim"
                continue
            cap = pool.q_max(job.z_lim_m)
            z2, o2, how = pool.hold_level(z1, o1, i1, i2, i3, job.z_lim_m, cap)
            if how == "rise":
                return pool.solve_outflow(z1, o1, i1, i2, pool.q_max)
            return z2, o2
        if phase == "draw_lim":
            assert job.z_hold_m is not None
            z2, o2, how = pool.release(z1, o1, i1, i2, i3, q_safe, job.z_hold_m)
            if how in ("land", "magic") or z2 <= job.z_hold_m + _LEVEL_TOL_M:
                pool.phase = "hold_lim"
            return z2, o2
        if phase == "hold_lim":
            assert job.z_hold_m is not None
            z2, o2, how = pool.hold_level(z1, o1, i1, i2, i3, job.z_hold_m, q_safe)
            if how == "rise":
                pool.phase = "pin"
                continue
            return z2, o2
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_PARAMETERS,
            f"unknown routing phase {phase!r}",
        )


def _opening_after(
    pool: _Pool,
    instant: int,
    stage: float,
    outflow: float,
) -> float | None:
    job = pool.job
    gates = job.gates
    if gates is None:
        return None
    if pool.phase in ("spill",) or (
        job.mode is LevelPoolMode.GATE_FAULT and job.q_safe_m3s is None
    ):
        return gates.max_opening_m
    if job.mode is LevelPoolMode.SCHEDULE:
        return _scheduled_opening(job, instant)
    if job.mode is LevelPoolMode.UNGATED_SPILL:
        return gates.max_opening_m
    command = gate_opening(
        gates,
        stage,
        max(0.0, outflow),
        n_fault=job.n_fault,
        stuck_opening_m=job.stuck_opening_m,
    )
    return command.opening_m


def _grid(job: _Job) -> list[int]:
    half = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    start = job.as_of_us - half
    end = job.as_of_us + half
    instants = [item[0] for item in job.inflow]
    instants.extend(item[0] for item in job.observed_stage)
    instants.extend(item[0] for item in job.observed_outflow)
    instants.extend(item[0] for item in job.observed_opening)
    for instant in instants:
        start = min(start, instant)
        end = max(end, instant)
    step = job.step_us
    # Snap outward onto the grid phased at as_of. An earlier off-grid instant
    # moves the start to the previous grid point, never forward past that data.
    start_remainder = (job.as_of_us - start) % step
    if start_remainder:
        start -= step - start_remainder
    end_mod = (end - job.as_of_us) % step
    if end_mod:
        end += step - end_mod
    count = (end - start) // step
    if count > _MAX_STEPS:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            f"routing window has {count} steps; the limit is {_MAX_STEPS}",
        )
    return [start + index * step for index in range(count + 1)]


def _route_job(job: _Job) -> LevelPoolResult:
    pool = _Pool(
        job=job,
        dt_s=job.step_us / 1_000_000.0,
        z_target=job.z0_m,
        z_hold=job.z0_m if job.z_hold_m is None else job.z_hold_m,
        phase=_initial_phase(job),
        reached_us=None,
    )
    if job.mode is LevelPoolMode.Q_SAFE and job.wait_z_m is None and job.z_hold_m is not None:
        pool.z_hold = job.z_hold_m
    _targets(pool)
    if job.mode is LevelPoolMode.Q_SAFE and job.wait_z_m is not None:
        pool.z_hold = job.wait_z_m
    grid = _grid(job)
    try:
        as_of_index = grid.index(job.as_of_us)
    except ValueError as exc:
        raise LevelPoolError(
            LevelPoolErrorCode.MISALIGNED_TIMES,
            "as_of is not on the routing grid",
        ) from exc
    inflows = [_sample(job.inflow, instant, outside=0.0) for instant in grid]
    stages = [0.0] * len(grid)
    outflows = [0.0] * len(grid)
    openings: list[float | None] = [None] * len(grid)
    gated = job.gates is not None
    for index, instant in enumerate(grid):
        if instant > job.as_of_us:
            break
        if instant == job.as_of_us:
            stages[index] = job.z0_m
            outflows[index] = _sample(job.observed_outflow, instant, outside=0.0)
            openings[index] = (
                None
                if not gated
                else _sample(job.observed_opening, instant, outside=0.0)
            )
            continue
        if job.observed_stage:
            stages[index] = _sample(job.observed_stage, instant, outside=job.z0_m)
        else:
            stages[index] = job.z0_m
        outflows[index] = _sample(job.observed_outflow, instant, outside=0.0)
        openings[index] = (
            None if not gated else _sample(job.observed_opening, instant, outside=0.0)
        )
        _require_stage_inside(job.storage, stages[index])
    for index in range(as_of_index + 1, len(grid)):
        i3 = inflows[index + 1] if index + 1 < len(grid) else 0.0
        z2, o2 = _step_pool(
            pool,
            instant=grid[index],
            z1=stages[index - 1],
            o1=outflows[index - 1],
            i1=inflows[index - 1],
            i2=inflows[index],
            i3=i3,
        )
        if o2 < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.ROOT_NOT_BRACKETED,
                f"outflow became negative at {_format_us(grid[index])}",
            )
        _require_stage_inside(job.storage, z2)
        stages[index] = z2
        outflows[index] = o2
        openings[index] = _opening_after(pool, grid[index], z2, o2)
    points: list[LevelPoolPoint] = []
    for index, instant in enumerate(grid):
        stage = stages[index]
        points.append(
            LevelPoolPoint(
                time=_format_us(instant),
                hours_from_as_of=(instant - job.as_of_us) / _US_PER_HOUR,
                stage_m=stage,
                storage_m3=storage_volume(job.storage, stage),
                inflow_m3s=inflows[index],
                outflow_m3s=outflows[index],
                gate_opening_m=openings[index],
            )
        )
    return _summarize(job, tuple(points), pool.reached_us)


def _summarize(
    job: _Job,
    points: tuple[LevelPoolPoint, ...],
    reached_us: int | None,
) -> LevelPoolResult:
    future = [point for point in points if point.hours_from_as_of >= 0.0]
    if not future:
        raise LevelPoolError(LevelPoolErrorCode.EMPTY_SERIES, "routing produced no instants")
    z_peak = future[0]
    o_peak = future[0]
    for point in future[1:]:
        if point.stage_m > z_peak.stage_m:
            z_peak = point
        if point.outflow_m3s > o_peak.outflow_m3s:
            o_peak = point
    i_peak = points[0]
    for point in points[1:]:
        if point.inflow_m3s > i_peak.inflow_m3s:
            i_peak = point
    reduction = i_peak.inflow_m3s - o_peak.outflow_m3s
    if i_peak.inflow_m3s <= 0.0:
        reduction_pct = 0.0
    else:
        reduction_pct = 100.0 * reduction / i_peak.inflow_m3s
    utilization: float | None = None
    if job.z_flood_limit_m is not None or job.z_design_m is not None:
        if job.z_flood_limit_m is None or job.z_design_m is None:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TARGET,
                "storage utilization needs both z_flood_limit_m and z_design_m",
            )
        v_flood = storage_volume(job.storage, job.z_flood_limit_m)
        v_design = storage_volume(job.storage, job.z_design_m)
        span = v_design - v_flood
        if span <= 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TARGET,
                "z_design_m must store more than z_flood_limit_m",
            )
        utilization = 100.0 * max(0.0, z_peak.storage_m3 - v_flood) / span
    reached = None if reached_us is None else _format_us(reached_us)
    return LevelPoolResult(
        points=points,
        as_of=job.as_of,
        timestep_hours=job.timestep_hours,
        z_max_m=z_peak.stage_m,
        t_z_max=z_peak.time,
        v_max_m3=z_peak.storage_m3,
        q_out_max_m3s=o_peak.outflow_m3s,
        t_q_out_max=o_peak.time,
        q_in_max_m3s=i_peak.inflow_m3s,
        t_q_in_max=i_peak.time,
        peak_reduction_m3s=reduction,
        peak_reduction_pct=reduction_pct,
        storage_utilization_pct=utilization,
        z_end_m=points[-1].stage_m,
        time_reached_target=reached,
    )


def _route(
    storage: object,
    inflow: object,
    *,
    mode: LevelPoolMode,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    discharge: DischargeCurve | None = None,
    gates: GateGeometry | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    q_safe_m3s: float | None = None,
    q_target_m3s: float | None = None,
    z_target_m: float | None = None,
    z_hold_m: float | None = None,
    z_lim_m: float | None = None,
    wait_closed_until_z: float | None = None,
    tail_inflow_m3s: float = 0.0,
    schedule: object = None,
    observed_stage: object = None,
    observed_outflow: object = None,
    observed_gate_opening: object = None,
    z_flood_limit_m: float | None = None,
    z_design_m: float | None = None,
) -> LevelPoolResult:
    job = _build_job(
        storage,
        inflow,
        mode=mode,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        discharge=discharge,
        gates=gates,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        q_safe_m3s=q_safe_m3s,
        q_target_m3s=q_target_m3s,
        z_target_m=z_target_m,
        z_hold_m=z_hold_m,
        z_lim_m=z_lim_m,
        wait_closed_until_z=wait_closed_until_z,
        tail_inflow_m3s=tail_inflow_m3s,
        schedule=schedule,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
        z_flood_limit_m=z_flood_limit_m,
        z_design_m=z_design_m,
    )
    return _route_job(job)


def route_level_pool(
    storage: StorageCurve,
    inflow: Sequence[tuple[str, float]],
    *,
    mode: LevelPoolMode,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    discharge: DischargeCurve | None = None,
    gates: GateGeometry | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    q_safe_m3s: float | None = None,
    q_target_m3s: float | None = None,
    z_target_m: float | None = None,
    z_hold_m: float | None = None,
    z_lim_m: float | None = None,
    wait_closed_until_z: float | None = None,
    tail_inflow_m3s: float = 0.0,
    schedule: Sequence[tuple[str, float]] | None = None,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
    z_flood_limit_m: float | None = None,
    z_design_m: float | None = None,
) -> LevelPoolResult:
    """Route one inflow through one outflow control.

    Exactly one of ``discharge`` and ``gates`` is accepted. ``mode`` selects
    ungated spill, a gate schedule, a constant target, pre-release, a safe
    discharge, a stage cap, full gate lift, or a stuck-gate case.
    """

    return _route(
        storage,
        inflow,
        mode=mode,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        discharge=discharge,
        gates=gates,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        q_safe_m3s=q_safe_m3s,
        q_target_m3s=q_target_m3s,
        z_target_m=z_target_m,
        z_hold_m=z_hold_m,
        z_lim_m=z_lim_m,
        wait_closed_until_z=wait_closed_until_z,
        tail_inflow_m3s=tail_inflow_m3s,
        schedule=schedule,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
        z_flood_limit_m=z_flood_limit_m,
        z_design_m=z_design_m,
    )


def _reject_extra_control(
    *,
    discharge: DischargeCurve | None,
    gates: GateGeometry | None,
) -> None:
    if discharge is not None and gates is not None:
        raise LevelPoolError(
            LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT,
            "provide exactly one outflow control: a discharge curve or gate geometry",
        )


def route_ungated(
    storage: StorageCurve,
    discharge: DischargeCurve,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    gates: GateGeometry | None = None,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    z_flood_limit_m: float | None = None,
    z_design_m: float | None = None,
) -> LevelPoolResult:
    """Modified Puls routing through an ungated discharge curve."""

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.UNGATED,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        discharge=discharge,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        z_flood_limit_m=z_flood_limit_m,
        z_design_m=z_design_m,
    )


def route_schedule(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    schedule: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    discharge: DischargeCurve | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Route with a step-held gate-opening schedule after ``as_of``."""

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.SCHEDULE,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        schedule=schedule,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def route_constant_target(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    q_target_m3s: float,
    discharge: DischargeCurve | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Release ``q_target_m3s`` whenever the gates can pass it."""

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.CONSTANT_TARGET,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        q_target_m3s=q_target_m3s,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def route_pre_release(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    q_safe_m3s: float,
    discharge: DischargeCurve | None = None,
    z_target_m: float | None = None,
    z_hold_m: float | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
    z_flood_limit_m: float | None = None,
    z_design_m: float | None = None,
) -> LevelPoolResult:
    """Draw the pool toward the sill, then hold the flood-limit stage.

    ``z_target_m`` defaults to the sill and is raised to the sill when a lower
    target is requested. After the rise, the pool holds ``z_hold_m`` (default
    ``z0_m``) and is not drawn back toward the sill, even when the sill was
    never reached.
    """

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.PRE_RELEASE,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        q_safe_m3s=q_safe_m3s,
        z_target_m=z_target_m,
        z_hold_m=z_hold_m,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
        z_flood_limit_m=z_flood_limit_m,
        z_design_m=z_design_m,
    )


def route_to_q_safe(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    q_safe_m3s: float,
    discharge: DischargeCurve | None = None,
    wait_closed_until_z: float | None = None,
    z_hold_m: float | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Cap outflow at the safe discharge.

    With ``wait_closed_until_z``, gates stay shut until that stage is forced
    upward, then release up to ``q_safe_m3s`` and hold the same stage on the
    recession. Without it, the pool is held at ``z_hold_m`` (default ``z0_m``)
    by passing inflow up to the cap, and is not drawn down in advance.
    """

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.Q_SAFE,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        q_safe_m3s=q_safe_m3s,
        wait_closed_until_z=wait_closed_until_z,
        z_hold_m=z_hold_m,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def route_to_z_max(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    z_lim_m: float,
    discharge: DischargeCurve | None = None,
    q_safe_m3s: float | None = None,
    z_hold_m: float | None = None,
    tail_inflow_m3s: float = 0.0,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Store with the gates closed until ``z_lim_m``, then match inflow.

    When ``q_safe_m3s`` and ``z_hold_m`` are both set, a recession that has
    fallen to ``tail_inflow_m3s`` is released at the safe discharge down to
    ``z_hold_m`` and held there.
    """

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.Z_MAX,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        z_lim_m=z_lim_m,
        q_safe_m3s=q_safe_m3s,
        z_hold_m=z_hold_m,
        tail_inflow_m3s=tail_inflow_m3s,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def route_ungated_spill(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    discharge: DischargeCurve | None = None,
    n_fault: int | None = None,
    stuck_opening_m: float = 0.0,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Hold every available gate fully open after ``as_of``."""

    _reject_extra_control(discharge=discharge, gates=gates)
    return _route(
        storage,
        inflow,
        mode=LevelPoolMode.UNGATED_SPILL,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def route_gate_fault(
    storage: StorageCurve,
    gates: GateGeometry,
    inflow: Sequence[tuple[str, float]],
    *,
    z0_m: float,
    as_of: str,
    timestep_hours: float,
    n_fault: int,
    discharge: DischargeCurve | None = None,
    stuck_opening_m: float = 0.0,
    q_safe_m3s: float | None = None,
    observed_stage: Sequence[tuple[str, float]] | None = None,
    observed_outflow: Sequence[tuple[str, float]] | None = None,
    observed_gate_opening: Sequence[tuple[str, float]] | None = None,
) -> LevelPoolResult:
    """Route with ``n_fault`` bays stuck and the rest spilling or capped.

    Stuck bays stay at ``stuck_opening_m`` (closed by default). The remaining
    bays spill fully open, or release ``q_safe_m3s`` when that cap is set.
    """

    _reject_extra_control(discharge=discharge, gates=gates)
    mode = LevelPoolMode.GATE_FAULT
    return _route(
        storage,
        inflow,
        mode=mode,
        z0_m=z0_m,
        as_of=as_of,
        timestep_hours=timestep_hours,
        gates=gates,
        n_fault=n_fault,
        stuck_opening_m=stuck_opening_m,
        q_safe_m3s=q_safe_m3s,
        q_target_m3s=q_safe_m3s,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_gate_opening,
    )


def reverse_inflow(
    storage: StorageCurve,
    observed_stage: Sequence[tuple[str, float]],
    observed_outflow: Sequence[tuple[str, float]],
    *,
    t_fail: str,
    as_of: str,
    timestep_hours: float,
    observed_inflow: Sequence[tuple[str, float]] | None = None,
    forecast_inflow: Sequence[tuple[str, float]] | None = None,
) -> ReservoirInflow:
    """Recover inflow on ``(t_fail, as_of]`` from stage and outflow.

    The interval average comes from level-pool continuity. Nodal inflow is the
    trapezoidal inverse of that average, starting from the observed inflow at
    ``t_fail``, then clipped at 0. Observed inflow at and before ``t_fail`` and
    forecast inflow after ``as_of`` are spliced on either side. The series
    covers at least ``[as_of - 12 h, as_of + 12 h]``.
    """

    _require_curve(storage)
    stages = _parse_series(
        observed_stage,
        "observed_stage",
        code=LevelPoolErrorCode.INVALID_STAGE,
    )
    outflows = _parse_series(
        observed_outflow,
        "observed_outflow",
        code=LevelPoolErrorCode.INVALID_DISCHARGE,
    )
    if len(stages) != len(outflows) or any(
        left[0] != right[0] for left, right in zip(stages, outflows, strict=True)
    ):
        raise LevelPoolError(
            LevelPoolErrorCode.MISALIGNED_TIMES,
            "observed_stage and observed_outflow must share the same instants",
        )
    for index, (_instant, outflow) in enumerate(outflows):
        if outflow < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"observed_outflow at index {index} must be >= 0 m³/s",
            )
    fail_us = _to_us(_parse_rfc3339(t_fail, role="t_fail"))
    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of"))
    if fail_us >= as_of_us:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTAMP,
            "t_fail must be earlier than as_of",
        )
    step = _step_us(timestep_hours)
    by_time = {
        instant: (stage, outflow)
        for (instant, stage), (_outflow_time, outflow) in zip(stages, outflows, strict=True)
    }
    if fail_us not in by_time or as_of_us not in by_time:
        raise LevelPoolError(
            LevelPoolErrorCode.MISALIGNED_TIMES,
            "observed stage and outflow must include both t_fail and as_of",
        )
    if (as_of_us - fail_us) % step != 0:
        raise LevelPoolError(
            LevelPoolErrorCode.MISALIGNED_TIMES,
            "t_fail and as_of must be an integer number of timesteps apart",
        )
    cursor = fail_us
    while cursor <= as_of_us:
        if cursor not in by_time:
            raise LevelPoolError(
                LevelPoolErrorCode.MISALIGNED_TIMES,
                f"observed stage and outflow are missing {_format_us(cursor)}",
            )
        _require_stage_inside(storage, by_time[cursor][0])
        cursor += step
    known = _parse_series(
        observed_inflow,
        "observed_inflow",
        code=LevelPoolErrorCode.INVALID_DISCHARGE,
        allow_empty=True,
    )
    forecast = _parse_series(
        forecast_inflow,
        "forecast_inflow",
        code=LevelPoolErrorCode.INVALID_DISCHARGE,
        allow_empty=True,
    )
    for index, (instant, discharge) in enumerate(known):
        if discharge < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"observed_inflow at index {index} must be >= 0 m³/s",
            )
        if instant > fail_us:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TIMESTAMP,
                "observed_inflow must not extend past t_fail; that gap is reconstructed",
            )
    known_at_fail = [discharge for instant, discharge in known if instant == fail_us]
    if not known_at_fail:
        raise LevelPoolError(
            LevelPoolErrorCode.MISALIGNED_TIMES,
            "observed_inflow must include t_fail so the gap can be reconstructed",
        )
    for index, (instant, discharge) in enumerate(forecast):
        if discharge < 0.0:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_DISCHARGE,
                f"forecast_inflow at index {index} must be >= 0 m³/s",
            )
        if instant <= as_of_us:
            raise LevelPoolError(
                LevelPoolErrorCode.INVALID_TIMESTAMP,
                "forecast_inflow must be later than as_of",
            )
    recovered: dict[int, float] = {fail_us: known_at_fail[0]}
    previous = known_at_fail[0]
    instant = fail_us + step
    dt_s = step / 1_000_000.0
    while instant <= as_of_us:
        stage_prev, outflow_prev = by_time[instant - step]
        stage, outflow = by_time[instant]
        average = 0.5 * (outflow + outflow_prev) + (
            storage_volume(storage, stage) - storage_volume(storage, stage_prev)
        ) / dt_s
        nodal = 2.0 * average - previous
        if nodal < 0.0:
            nodal = 0.0
        recovered[instant] = nodal
        previous = nodal
        instant += step
    half = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    start = min(as_of_us - half, fail_us, known[0][0])
    end = as_of_us + half
    if forecast:
        end = max(end, forecast[-1][0])
    start_remainder = (as_of_us - start) % step
    if start_remainder:
        start -= step - start_remainder
    if (end - as_of_us) % step != 0:
        end += step - (end - as_of_us) % step
    count = (end - start) // step
    if count > _MAX_STEPS:
        raise LevelPoolError(
            LevelPoolErrorCode.INVALID_TIMESTEP,
            f"reverse inflow window has {count} steps; the limit is {_MAX_STEPS}",
        )
    points: list[ReservoirInflowPoint] = []
    for index in range(count + 1):
        instant = start + index * step
        if instant in recovered and instant > fail_us:
            discharge = recovered[instant]
        elif instant <= fail_us:
            discharge = _sample(known, instant, outside=0.0)
        else:
            discharge = _sample(forecast, instant, outside=0.0)
        points.append(
            ReservoirInflowPoint(
                time=_format_us(instant),
                hours_from_as_of=(instant - as_of_us) / _US_PER_HOUR,
                discharge_m3s=discharge,
            )
        )
    return ReservoirInflow(
        points=tuple(points),
        as_of=as_of,
        t_fail=t_fail,
        timestep_hours=float(timestep_hours),
    )

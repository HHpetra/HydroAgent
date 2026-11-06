"""SCS triangular unit hydrograph and degraded triangular inflow.

Excess is in millimetres, discharge in m³/s, area in km², and clock times are
RFC3339 strings. Nothing in this module reads the wall clock. A scenario
``as_of`` or a join instant is an argument.

One millimetre of excess over ``F`` km² yields ``1000 F`` m³ after the discrete
unit hydrograph is normalized. ``t_join`` selects remainder excess and keeps
the inflow continuous at the measured join discharge. ``as_of`` anchors the
degraded triangle and pads the output so that ``[as_of - 12 h, as_of + 12 h]``
is present; slots with no hydrograph value are 0 m³/s.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from enum import StrEnum
from typing import TypeAlias

from hydro_agent.core.scs_cn import (
    ExcessStep,
    RainfallSeries,
    ScsCnParameters,
    TimedRainfall,
    compute_excess,
)

_US_PER_HOUR = 3_600_000_000
_WINDOW_HALF_HOURS = 12.0
_SCS_RECESSION_RATIO = 1.67
_TRIANGLE_RECESSION_TO_RISE = 2.0
_LAG_FRACTION_OF_TC = 0.6
_MAX_STEPS = 1_000_000

_RFC3339 = re.compile(r"^\d{4}-\d{2}-\d{2}[Tt]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:[Zz]|[+-]\d{2}:\d{2})$")

_EPOCH = datetime(1970, 1, 1, tzinfo=UTC)


class ScsUhErrorCode(StrEnum):
    """Stable identifier for an SCS unit-hydrograph domain failure."""

    INVALID_AREA = "invalid_area"
    INVALID_TIMESTEP = "invalid_timestep"
    INVALID_TIMING = "invalid_timing"
    INVALID_LAG = "invalid_lag"
    INVALID_TIME_OF_CONCENTRATION = "invalid_time_of_concentration"
    INVALID_LENGTH = "invalid_length"
    INVALID_SLOPE = "invalid_slope"
    INVALID_CURVE_NUMBER = "invalid_curve_number"
    INVALID_RECESSION_RATIO = "invalid_recession_ratio"
    INVALID_EXCESS = "invalid_excess"
    EMPTY_EXCESS = "empty_excess"
    INVALID_TIMESTAMP = "invalid_timestamp"
    UNORDERED_TIMESTAMPS = "unordered_timestamps"
    MISALIGNED_TIMES = "misaligned_times"
    INVALID_PEAK = "invalid_peak"
    INVALID_HOURS_TO_PEAK = "invalid_hours_to_peak"
    INVALID_FALL_HOURS = "invalid_fall_hours"
    INVALID_DISCHARGE = "invalid_discharge"
    INVALID_BASEFLOW = "invalid_baseflow"
    INVALID_JOIN = "invalid_join"
    INVALID_REGIONAL = "invalid_regional"


class ScsUhError(Exception):
    """Illegal unit-hydrograph input.

    ``code`` identifies the failure. ``message`` describes the offending value.
    """

    def __init__(self, code: ScsUhErrorCode, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code.value}: {message}")


class TimingMethod(StrEnum):
    """How watershed geometry or a given duration becomes the SCS lag."""

    LAG = "lag"
    TIME_OF_CONCENTRATION = "time_of_concentration"
    NRCS_LAG = "nrcs_lag"
    KIRPICH = "kirpich"
    REGIONAL = "regional"


class RegionalQuantity(StrEnum):
    """What a regional power formula returns."""

    LAG = "lag"
    TIME_OF_CONCENTRATION = "time_of_concentration"
    RISING_TIME = "rising_time"


def _require_real(value: object, code: ScsUhErrorCode, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ScsUhError(code, f"{name} must be a real number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise ScsUhError(code, f"{name} must be finite, got {value!r}")
    return number


def _parse_rfc3339(value: object, *, role: str) -> datetime:
    if not isinstance(value, str) or _RFC3339.fullmatch(value) is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTAMP,
            f"{role} must be an RFC3339 timestamp, got {value!r}",
        )
    normalized = f"{value[:-1]}+00:00" if value[-1] in "Zz" else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTAMP,
            f"{role} is not a valid RFC3339 timestamp: {value}",
        ) from exc
    if parsed.tzinfo is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTAMP,
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
    number = _require_real(timestep_hours, ScsUhErrorCode.INVALID_TIMESTEP, "timestep_hours")
    if number <= 0.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be > 0, got {timestep_hours!r}",
        )
    step = round(number * _US_PER_HOUR)
    if step <= 0 or abs(number * _US_PER_HOUR - step) > 1e-3:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be an exact microsecond count, got {timestep_hours!r}",
        )
    return step


def _require_divides_window(step: int) -> None:
    window = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    if window % step != 0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTEP,
            "timestep_hours must divide the 12 h scenario window",
        )


def _require_multiple(hours: float, step: int, *, code: ScsUhErrorCode, name: str) -> int:
    span = round(hours * _US_PER_HOUR)
    if span <= 0 or abs(hours * _US_PER_HOUR - span) > 1e-3 or span % step != 0:
        raise ScsUhError(
            code,
            f"{name} must be a positive multiple of timestep_hours, got {hours!r}",
        )
    return span


def _enum_value(
    enum_type: type[StrEnum],
    value: object,
    code: ScsUhErrorCode,
    name: str,
) -> StrEnum:
    if isinstance(value, enum_type):
        return value
    if isinstance(value, str):
        try:
            return enum_type(value)
        except ValueError:
            pass
    allowed = ", ".join(item.value for item in enum_type)
    raise ScsUhError(code, f"{name} must be one of {allowed}, got {value!r}")


@dataclass(frozen=True, slots=True)
class RegionalTiming:
    """Regional power formula ``c · F^α · L^β · J^γ`` in hours.

    ``quantity`` says whether the result is lag, time of concentration, or the
    unit-hydrograph rising time τ. Rising time is the SCS time to peak.
    """

    coefficient: float
    area_exponent: float = 0.0
    length_exponent: float = 0.0
    slope_exponent: float = 0.0
    quantity: RegionalQuantity | str = RegionalQuantity.RISING_TIME

    def __post_init__(self) -> None:
        coefficient = _require_real(
            self.coefficient,
            ScsUhErrorCode.INVALID_REGIONAL,
            "coefficient",
        )
        if coefficient <= 0.0:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_REGIONAL,
                f"coefficient must be > 0, got {self.coefficient!r}",
            )
        object.__setattr__(self, "coefficient", coefficient)
        for name in ("area_exponent", "length_exponent", "slope_exponent"):
            object.__setattr__(
                self,
                name,
                _require_real(getattr(self, name), ScsUhErrorCode.INVALID_REGIONAL, name),
            )
        object.__setattr__(
            self,
            "quantity",
            _enum_value(
                RegionalQuantity,
                self.quantity,
                ScsUhErrorCode.INVALID_REGIONAL,
                "quantity",
            ),
        )


@dataclass(frozen=True, slots=True)
class Watershed:
    """Area and the inputs that determine SCS lag.

    When ``timing_method`` is omitted, a given lag is used first, then a given
    time of concentration, then a regional formula, then NRCS lag when a curve
    number is present, then Kirpich.
    """

    area_km2: float
    lag_hours: float | None = None
    time_of_concentration_hours: float | None = None
    channel_length_km: float | None = None
    slope: float | None = None
    curve_number: float | None = None
    timing_method: TimingMethod | str | None = None
    regional: RegionalTiming | None = None

    def __post_init__(self) -> None:
        area = _require_real(self.area_km2, ScsUhErrorCode.INVALID_AREA, "area_km2")
        if area <= 0.0:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_AREA,
                f"area_km2 must be > 0, got {self.area_km2!r}",
            )
        object.__setattr__(self, "area_km2", area)
        object.__setattr__(
            self,
            "lag_hours",
            _optional_duration(self.lag_hours, "lag_hours", allow_zero=True),
        )
        object.__setattr__(
            self,
            "time_of_concentration_hours",
            _optional_duration(
                self.time_of_concentration_hours,
                "time_of_concentration_hours",
                allow_zero=False,
            ),
        )
        object.__setattr__(
            self,
            "channel_length_km",
            _optional_positive(
                self.channel_length_km,
                ScsUhErrorCode.INVALID_LENGTH,
                "channel_length_km",
            ),
        )
        object.__setattr__(
            self,
            "slope",
            _optional_positive(self.slope, ScsUhErrorCode.INVALID_SLOPE, "slope"),
        )
        object.__setattr__(self, "curve_number", _optional_curve_number(self.curve_number))
        if self.timing_method is not None:
            object.__setattr__(
                self,
                "timing_method",
                _enum_value(
                    TimingMethod,
                    self.timing_method,
                    ScsUhErrorCode.INVALID_TIMING,
                    "timing_method",
                ),
            )
        if self.regional is not None and not isinstance(self.regional, RegionalTiming):
            raise ScsUhError(
                ScsUhErrorCode.INVALID_REGIONAL,
                f"regional must be RegionalTiming, got {self.regional!r}",
            )


def _optional_duration(value: object, name: str, *, allow_zero: bool) -> float | None:
    if value is None:
        return None
    code = (
        ScsUhErrorCode.INVALID_LAG
        if name == "lag_hours"
        else ScsUhErrorCode.INVALID_TIME_OF_CONCENTRATION
    )
    number = _require_real(value, code, name)
    if number < 0.0 or (number == 0.0 and not allow_zero):
        raise ScsUhError(code, f"{name} must be {'>= 0' if allow_zero else '> 0'}, got {value!r}")
    return number


def _optional_positive(value: object, code: ScsUhErrorCode, name: str) -> float | None:
    if value is None:
        return None
    number = _require_real(value, code, name)
    if number <= 0.0:
        raise ScsUhError(code, f"{name} must be > 0, got {value!r}")
    return number


def _optional_curve_number(value: object) -> float | None:
    if value is None:
        return None
    number = _require_real(value, ScsUhErrorCode.INVALID_CURVE_NUMBER, "curve_number")
    if not 0.0 < number <= 100.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_CURVE_NUMBER,
            f"curve_number must be in (0, 100], got {value!r}",
        )
    return number


@dataclass(frozen=True, slots=True)
class UnitHydrographOrdinate:
    """One discrete SCS ordinate for 1 mm of excess."""

    hours: float
    discharge_m3s: float
    dimensionless: float


@dataclass(frozen=True, slots=True)
class ScsUnitHydrograph:
    """SCS triangular unit hydrograph for 1 mm of excess.

    ``peak_discharge_m3s`` is the analytical triangle peak. ``ordinates`` are
    the discrete samples after rescaling so their rectangular volume is
    ``1000 * area_km2`` cubic metres.
    """

    area_km2: float
    timestep_hours: float
    lag_hours: float
    time_of_concentration_hours: float
    time_to_peak_hours: float
    recession_hours: float
    base_hours: float
    peak_discharge_m3s: float
    peak_rate_coefficient: float
    recession_ratio: float
    ordinates: tuple[UnitHydrographOrdinate, ...]


@dataclass(frozen=True, slots=True)
class InflowPoint:
    """One inflow ordinate.

    ``time`` is omitted when the series was built from bare excess depths.
    ``hours_from_origin`` is measured from ``as_of`` when that instant is set,
    otherwise from the join, otherwise from the start of the first excess
    interval.
    """

    time: str | None
    hours_from_origin: float
    discharge_m3s: float


@dataclass(frozen=True, slots=True)
class InflowHydrograph:
    """Inflow at the dam, with the peak and the rectangular volume of the ordinates.

    ``volume_m3`` is ``Σ Q · 3600 · Δt`` over the returned ordinates. Empty
    scenario-window slots are included as 0 m³/s and do not change the volume.
    """

    points: tuple[InflowPoint, ...]
    peak_discharge_m3s: float
    peak_time: str | None
    peak_hours_from_origin: float
    volume_m3: float
    timestep_hours: float
    as_of: str | None

    def discharge_at(self, time: str) -> float:
        """Discharge at an RFC3339 instant that the series contains."""

        target = _to_us(_parse_rfc3339(time, role="discharge time"))
        for point in self.points:
            if point.time is None:
                continue
            point_us = _to_us(_parse_rfc3339(point.time, role="inflow time"))
            if point_us == target:
                return point.discharge_m3s
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTAMP,
            f"inflow series has no ordinate at {time}",
        )


ExcessInput: TypeAlias = (
    Sequence[float] | Sequence[TimedRainfall] | Sequence[tuple[str, float]] | Sequence[ExcessStep]
)


def _require_length(watershed: Watershed) -> float:
    if watershed.channel_length_km is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_LENGTH,
            "channel_length_km is required for this timing method",
        )
    return watershed.channel_length_km


def _require_slope(watershed: Watershed) -> float:
    if watershed.slope is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_SLOPE,
            "slope is required for this timing method",
        )
    return watershed.slope


def _require_cn(watershed: Watershed) -> float:
    if watershed.curve_number is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_CURVE_NUMBER,
            "curve_number is required for NRCS lag",
        )
    return watershed.curve_number


def nrcs_lag_hours(channel_length_km: float, slope: float, curve_number: float) -> float:
    """NRCS lag in hours. Length is km, slope is m/m, retention is derived from CN."""

    retention_mm = 25400.0 / curve_number - 254.0
    length_ft = 3280.84 * channel_length_km
    slope_percent = 100.0 * slope
    retention_term = retention_mm / 25.4 + 1.0
    return (math.pow(length_ft, 0.8) * math.pow(retention_term, 0.7)) / (
        1900.0 * math.sqrt(slope_percent)
    )


def kirpich_tc_hours(channel_length_km: float, slope: float) -> float:
    """Kirpich time of concentration in hours. Length is km and slope is m/m."""

    return 0.0663 * math.pow(channel_length_km, 0.77) * math.pow(slope, -0.385)


def _regional_hours(watershed: Watershed, regional: RegionalTiming) -> float:
    value = regional.coefficient
    if regional.area_exponent != 0.0:
        value *= math.pow(watershed.area_km2, regional.area_exponent)
    if regional.length_exponent != 0.0:
        value *= math.pow(_require_length(watershed), regional.length_exponent)
    elif watershed.channel_length_km is not None:
        value *= math.pow(watershed.channel_length_km, regional.length_exponent)
    if regional.slope_exponent != 0.0:
        value *= math.pow(_require_slope(watershed), regional.slope_exponent)
    elif watershed.slope is not None:
        value *= math.pow(watershed.slope, regional.slope_exponent)
    if not math.isfinite(value) or value <= 0.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_REGIONAL,
            f"regional timing formula returned {value!r}",
        )
    return value


def _infer_method(watershed: Watershed) -> TimingMethod:
    if isinstance(watershed.timing_method, TimingMethod):
        return watershed.timing_method
    if watershed.lag_hours is not None:
        return TimingMethod.LAG
    if watershed.time_of_concentration_hours is not None:
        return TimingMethod.TIME_OF_CONCENTRATION
    if watershed.regional is not None:
        return TimingMethod.REGIONAL
    if (
        watershed.channel_length_km is not None
        and watershed.slope is not None
        and watershed.curve_number is not None
    ):
        return TimingMethod.NRCS_LAG
    if watershed.channel_length_km is not None and watershed.slope is not None:
        return TimingMethod.KIRPICH
    raise ScsUhError(
        ScsUhErrorCode.INVALID_TIMING,
        "provide lag_hours, time_of_concentration_hours, or watershed geometry",
    )


def _resolve_lag_hours(watershed: Watershed, timestep_hours: float) -> tuple[float, float]:
    """Return SCS lag and the reported time of concentration, both in hours."""

    method = _infer_method(watershed)
    if method is TimingMethod.LAG:
        if watershed.lag_hours is None:
            raise ScsUhError(ScsUhErrorCode.INVALID_LAG, "lag_hours is required")
        lag = watershed.lag_hours
    elif method is TimingMethod.TIME_OF_CONCENTRATION:
        if watershed.time_of_concentration_hours is None:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_TIME_OF_CONCENTRATION,
                "time_of_concentration_hours is required",
            )
        lag = _LAG_FRACTION_OF_TC * watershed.time_of_concentration_hours
    elif method is TimingMethod.NRCS_LAG:
        lag = nrcs_lag_hours(
            _require_length(watershed),
            _require_slope(watershed),
            _require_cn(watershed),
        )
    elif method is TimingMethod.KIRPICH:
        lag = _LAG_FRACTION_OF_TC * kirpich_tc_hours(
            _require_length(watershed),
            _require_slope(watershed),
        )
    else:
        if watershed.regional is None:
            raise ScsUhError(ScsUhErrorCode.INVALID_REGIONAL, "regional timing is required")
        regional = watershed.regional
        hours = _regional_hours(watershed, regional)
        quantity = regional.quantity
        if quantity is RegionalQuantity.LAG:
            lag = hours
        elif quantity is RegionalQuantity.TIME_OF_CONCENTRATION:
            lag = _LAG_FRACTION_OF_TC * hours
        else:
            lag = hours - timestep_hours / 2.0
            if lag < 0.0:
                raise ScsUhError(
                    ScsUhErrorCode.INVALID_REGIONAL,
                    "regional rising time must be at least half a timestep",
                )
    if not math.isfinite(lag) or lag < 0.0:
        raise ScsUhError(ScsUhErrorCode.INVALID_LAG, f"lag must be >= 0 h, got {lag!r}")
    return lag, lag / _LAG_FRACTION_OF_TC


def _triangle(time_hours: float, time_to_peak: float, base_hours: float, peak: float) -> float:
    if time_hours <= 0.0 or time_hours >= base_hours:
        return 0.0
    if time_hours <= time_to_peak:
        return peak * time_hours / time_to_peak
    return peak * (base_hours - time_hours) / (base_hours - time_to_peak)


def scs_unit_hydrograph(
    watershed: Watershed,
    *,
    timestep_hours: float,
    recession_ratio: float = _SCS_RECESSION_RATIO,
) -> ScsUnitHydrograph:
    """Build the discrete SCS triangle for 1 mm, normalized to ``1000 F`` m³."""

    _step_us(timestep_hours)
    ratio = _require_real(
        recession_ratio,
        ScsUhErrorCode.INVALID_RECESSION_RATIO,
        "recession_ratio",
    )
    if ratio <= 0.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_RECESSION_RATIO,
            f"recession_ratio must be > 0, got {recession_ratio!r}",
        )
    lag, time_of_concentration = _resolve_lag_hours(watershed, timestep_hours)
    time_to_peak = timestep_hours / 2.0 + lag
    if time_to_peak <= 0.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_LAG,
            "time to peak must be > 0; increase lag or the timestep",
        )
    recession = ratio * time_to_peak
    base = time_to_peak + recession
    peak_rate = (1.0 / 3.6) * 2.0 / (1.0 + ratio)
    peak = peak_rate * watershed.area_km2 / time_to_peak

    steps = 0
    while steps * timestep_hours < base - 1e-9:
        steps += 1
        if steps > _MAX_STEPS:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_TIMESTEP,
                "unit hydrograph has more than 1000000 ordinates",
            )
    raw = [
        _triangle(index * timestep_hours, time_to_peak, base, peak) for index in range(steps + 1)
    ]
    weight = math.fsum(raw)
    if weight <= 0.0:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTEP,
            "discrete triangle samples are all zero; reduce the timestep",
        )
    target_sum = (1000.0 * watershed.area_km2) / (3600.0 * timestep_hours)
    scale = target_sum / weight
    scaled = [value * scale for value in raw]
    scaled[max(range(len(scaled)), key=scaled.__getitem__)] += target_sum - math.fsum(scaled)
    ordinates = tuple(
        UnitHydrographOrdinate(
            hours=index * timestep_hours,
            discharge_m3s=discharge,
            dimensionless=0.0 if peak == 0.0 else raw[index] / peak,
        )
        for index, discharge in enumerate(scaled)
    )
    return ScsUnitHydrograph(
        area_km2=watershed.area_km2,
        timestep_hours=timestep_hours,
        lag_hours=lag,
        time_of_concentration_hours=time_of_concentration,
        time_to_peak_hours=time_to_peak,
        recession_hours=recession,
        base_hours=base,
        peak_discharge_m3s=peak,
        peak_rate_coefficient=peak_rate,
        recession_ratio=ratio,
        ordinates=ordinates,
    )


def _require_discharge(value: object, code: ScsUhErrorCode, name: str) -> float:
    number = _require_real(value, code, name)
    if number < 0.0:
        raise ScsUhError(code, f"{name} must be >= 0, got {value!r}")
    return number


def _observed_map(
    observed: Sequence[tuple[str, float]] | None,
    step: int,
    phase: int | None,
) -> dict[int, float]:
    if observed is None:
        return {}
    if isinstance(observed, str | bytes) or not isinstance(observed, Sequence):
        raise ScsUhError(
            ScsUhErrorCode.INVALID_DISCHARGE,
            "observed_inflow must be a sequence of (time, discharge) pairs",
        )
    mapped: dict[int, float] = {}
    previous: int | None = None
    for index, item in enumerate(observed):
        if not isinstance(item, tuple) or len(item) != 2:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_DISCHARGE,
                f"observed_inflow at index {index} must be (time, discharge)",
            )
        instant = _to_us(_parse_rfc3339(item[0], role=f"observed inflow time at index {index}"))
        discharge = _require_discharge(
            item[1],
            ScsUhErrorCode.INVALID_DISCHARGE,
            f"observed inflow at index {index}",
        )
        if previous is not None and instant <= previous:
            raise ScsUhError(
                ScsUhErrorCode.UNORDERED_TIMESTAMPS,
                f"observed inflow time at index {index} must be strictly after the previous time",
            )
        if phase is not None and (instant - phase) % step != 0:
            raise ScsUhError(
                ScsUhErrorCode.MISALIGNED_TIMES,
                f"observed inflow at index {index} is not on the timestep grid",
            )
        previous = instant
        mapped[instant] = discharge
    return mapped


def _pulses_from_excess(excess: ExcessInput) -> list[tuple[str | None, float]]:
    if isinstance(excess, str | bytes) or not isinstance(excess, Sequence):
        raise ScsUhError(
            ScsUhErrorCode.INVALID_EXCESS,
            "excess must be a sequence of depths or timed excess pulses",
        )
    if len(excess) == 0:
        raise ScsUhError(ScsUhErrorCode.EMPTY_EXCESS, "excess series is empty")
    first = excess[0]
    if isinstance(first, TimedRainfall):
        return _timed_pulses(excess, TimedRainfall)
    if isinstance(first, ExcessStep):
        return _excess_steps(excess)
    if isinstance(first, tuple):
        return _pair_pulses(excess)
    return _bare_pulses(excess)


def _timed_pulses(
    excess: Sequence[object],
    kind: type[TimedRainfall],
) -> list[tuple[str | None, float]]:
    pulses: list[tuple[str | None, float]] = []
    previous: datetime | None = None
    for index, item in enumerate(excess):
        if not isinstance(item, kind):
            raise ScsUhError(
                ScsUhErrorCode.INVALID_EXCESS,
                f"excess at index {index} must be TimedRainfall, matching the first pulse",
            )
        instant = _parse_rfc3339(item.time, role=f"excess time at index {index}")
        if previous is not None and instant <= previous:
            raise ScsUhError(
                ScsUhErrorCode.UNORDERED_TIMESTAMPS,
                f"excess time at index {index} must be strictly after the previous time",
            )
        previous = instant
        pulses.append((item.time, item.depth_mm))
    return pulses


def _excess_steps(excess: Sequence[object]) -> list[tuple[str | None, float]]:
    pulses: list[tuple[str | None, float]] = []
    seen_time = False
    seen_bare = False
    previous: datetime | None = None
    for index, item in enumerate(excess):
        if not isinstance(item, ExcessStep):
            raise ScsUhError(
                ScsUhErrorCode.INVALID_EXCESS,
                f"excess at index {index} must be ExcessStep, matching the first pulse",
            )
        depth = _require_discharge(
            item.excess_mm,
            ScsUhErrorCode.INVALID_EXCESS,
            f"excess at index {index}",
        )
        if item.time is None:
            seen_bare = True
            pulses.append((None, depth))
            continue
        seen_time = True
        instant = _parse_rfc3339(item.time, role=f"excess time at index {index}")
        if previous is not None and instant <= previous:
            raise ScsUhError(
                ScsUhErrorCode.UNORDERED_TIMESTAMPS,
                f"excess time at index {index} must be strictly after the previous time",
            )
        previous = instant
        pulses.append((item.time, depth))
    if seen_time and seen_bare:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_EXCESS,
            "excess steps must be all timed or all bare",
        )
    return pulses


def _pair_pulses(excess: Sequence[object]) -> list[tuple[str | None, float]]:
    pulses: list[tuple[str | None, float]] = []
    previous: datetime | None = None
    for index, item in enumerate(excess):
        if not isinstance(item, tuple) or len(item) != 2 or not isinstance(item[0], str):
            raise ScsUhError(
                ScsUhErrorCode.INVALID_EXCESS,
                f"excess at index {index} must be (time, excess_mm)",
            )
        instant = _parse_rfc3339(item[0], role=f"excess time at index {index}")
        if previous is not None and instant <= previous:
            raise ScsUhError(
                ScsUhErrorCode.UNORDERED_TIMESTAMPS,
                f"excess time at index {index} must be strictly after the previous time",
            )
        previous = instant
        depth = _require_discharge(
            item[1],
            ScsUhErrorCode.INVALID_EXCESS,
            f"excess at index {index}",
        )
        pulses.append((item[0], depth))
    return pulses


def _bare_pulses(excess: Sequence[object]) -> list[tuple[str | None, float]]:
    pulses: list[tuple[str | None, float]] = []
    for index, item in enumerate(excess):
        if isinstance(item, TimedRainfall | ExcessStep | tuple):
            raise ScsUhError(
                ScsUhErrorCode.INVALID_EXCESS,
                f"excess at index {index} is timed, but the series started with a bare depth",
            )
        pulses.append(
            (
                None,
                _require_discharge(
                    item,
                    ScsUhErrorCode.INVALID_EXCESS,
                    f"excess at index {index}",
                ),
            )
        )
    return pulses


def _direct_runoff(
    pulse_ends: Sequence[tuple[int, float]],
    ordinates: Sequence[float],
    step: int,
) -> dict[int, float]:
    """Convolution ``Q(n) = Σ ΔR_k U_{n-k+1}`` keyed by absolute microseconds."""

    response: dict[int, float] = {}
    for end_us, excess_mm in pulse_ends:
        if excess_mm == 0.0:
            continue
        start_us = end_us - step
        for index, ordinate in enumerate(ordinates):
            if ordinate == 0.0:
                continue
            key = start_us + index * step
            response[key] = response.get(key, 0.0) + excess_mm * ordinate
    return response


def _dense(
    start: int,
    end: int,
    step: int,
    values: dict[int, float],
    baseline: float,
) -> dict[int, float]:
    filled: dict[int, float] = {}
    count = (end - start) // step
    for index in range(count + 1):
        key = start + index * step
        filled[key] = values.get(key, 0.0) + baseline
    return filled


def _window_times(as_of_us: int, step: int) -> list[int]:
    half = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    start = as_of_us - half
    count = (2 * half) // step
    return [start + index * step for index in range(count + 1)]


def _value_at(
    instant: int,
    *,
    series: dict[int, float],
    observed: dict[int, float],
    splice: int | None,
    in_window: bool,
) -> float | None:
    if splice is not None and instant < splice:
        if instant in observed:
            return observed[instant]
        return 0.0 if in_window else None
    if splice is not None and instant == splice:
        if instant in observed:
            return observed[instant]
        if instant in series:
            return series[instant]
        return 0.0 if in_window else None
    if instant in series:
        return series[instant]
    return 0.0 if in_window else None


def _finish(
    points: list[InflowPoint],
    timestep_hours: float,
    as_of: str | None,
) -> InflowHydrograph:
    if not points:
        raise ScsUhError(ScsUhErrorCode.EMPTY_EXCESS, "inflow series is empty")
    peak = points[0]
    for point in points[1:]:
        if point.discharge_m3s > peak.discharge_m3s:
            peak = point
    volume = math.fsum(point.discharge_m3s for point in points) * 3600.0 * timestep_hours
    return InflowHydrograph(
        points=tuple(points),
        peak_discharge_m3s=peak.discharge_m3s,
        peak_time=peak.time,
        peak_hours_from_origin=peak.hours_from_origin,
        volume_m3=volume,
        timestep_hours=timestep_hours,
        as_of=as_of,
    )


def _points_from_samples(
    samples: Sequence[tuple[int, float]],
    *,
    origin_us: int,
    timestep_hours: float,
    as_of: str | None,
    calendar: bool,
) -> InflowHydrograph:
    points = [
        InflowPoint(
            time=_format_us(instant) if calendar else None,
            hours_from_origin=(instant - origin_us) / _US_PER_HOUR,
            discharge_m3s=discharge,
        )
        for instant, discharge in samples
    ]
    return _finish(points, timestep_hours, as_of)


def _collect_samples(
    series: dict[int, float],
    observed: dict[int, float],
    *,
    splice: int | None,
    as_of_us: int | None,
    step: int,
) -> list[tuple[int, float]]:
    if as_of_us is None:
        keys = set(series)
        for instant in observed:
            if splice is None or instant < splice:
                keys.add(instant)
        if splice is not None:
            keys.add(splice)
        ordered = sorted(keys)
        in_window = {instant: False for instant in ordered}
    else:
        window = _window_times(as_of_us, step)
        window_set = set(window)
        keys = set(window)
        keys.update(series)
        if splice is not None:
            keys.update(instant for instant in observed if instant <= splice)
        ordered = sorted(keys)
        in_window = {instant: instant in window_set for instant in ordered}

    samples: list[tuple[int, float]] = []
    for instant in ordered:
        discharge = _value_at(
            instant,
            series=series,
            observed=observed,
            splice=splice,
            in_window=in_window[instant],
        )
        if discharge is not None:
            samples.append((instant, discharge))
    return samples


def route_excess(
    excess: ExcessInput,
    watershed: Watershed,
    *,
    timestep_hours: float,
    recession_ratio: float = _SCS_RECESSION_RATIO,
    baseflow_m3s: float = 0.0,
    t_join: str | None = None,
    as_of: str | None = None,
    observed_inflow: Sequence[tuple[str, float]] | None = None,
) -> InflowHydrograph:
    """Convolve excess rainfall with the SCS triangular unit hydrograph.

    A bare depth sequence is one pulse per timestep, starting at hour 0.
    Timed pulses are interval-end excess in millimetres. When ``t_join`` is
    set, only pulses that end after that instant are convolved, and the result
    is continuous with the observed discharge at the join. ``baseflow_m3s`` is
    added on the convolved support when there is no join; it is not copied
    into empty scenario-window slots.
    """

    pulses = _pulses_from_excess(excess)
    return _route_pulses(
        pulses,
        watershed,
        timestep_hours=timestep_hours,
        recession_ratio=recession_ratio,
        baseflow_m3s=baseflow_m3s,
        t_join=t_join,
        as_of=as_of,
        observed_inflow=observed_inflow,
    )


def route_rainfall(
    rainfall: RainfallSeries,
    loss: ScsCnParameters,
    watershed: Watershed,
    *,
    timestep_hours: float,
    recession_ratio: float = _SCS_RECESSION_RATIO,
    baseflow_m3s: float = 0.0,
    t_join: str | None = None,
    as_of: str | None = None,
    observed_inflow: Sequence[tuple[str, float]] | None = None,
) -> InflowHydrograph:
    """Run SCS-CN on the given hyetograph, then SCS-TUH convolution.

    The inflow is the convolution of that storm's incremental excess. A storm
    total is never substituted for the hyetograph. When ``t_join`` is set, the
    full storm is still used to consume initial abstraction, and only excess
    increments after ``t_join`` are convolved.
    """

    storm = compute_excess(rainfall, loss)
    pulses = [(step.time, step.excess_mm) for step in storm.steps]
    return _route_pulses(
        pulses,
        watershed,
        timestep_hours=timestep_hours,
        recession_ratio=recession_ratio,
        baseflow_m3s=baseflow_m3s,
        t_join=t_join,
        as_of=as_of,
        observed_inflow=observed_inflow,
    )


def _route_pulses(
    pulses: Sequence[tuple[str | None, float]],
    watershed: Watershed,
    *,
    timestep_hours: float,
    recession_ratio: float,
    baseflow_m3s: float,
    t_join: str | None,
    as_of: str | None,
    observed_inflow: Sequence[tuple[str, float]] | None,
) -> InflowHydrograph:
    if not pulses:
        raise ScsUhError(ScsUhErrorCode.EMPTY_EXCESS, "excess series is empty")
    step = _step_us(timestep_hours)
    baseflow = _require_discharge(baseflow_m3s, ScsUhErrorCode.INVALID_BASEFLOW, "baseflow_m3s")
    unit = scs_unit_hydrograph(
        watershed,
        timestep_hours=timestep_hours,
        recession_ratio=recession_ratio,
    )
    ordinates = tuple(item.discharge_m3s for item in unit.ordinates)
    calendar = pulses[0][0] is not None
    if any((pulse[0] is None) != (not calendar) for pulse in pulses):
        raise ScsUhError(
            ScsUhErrorCode.INVALID_EXCESS,
            "excess pulses must be all timed or all bare",
        )

    join_us = _to_us(_parse_rfc3339(t_join, role="t_join")) if t_join is not None else None
    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of")) if as_of is not None else None
    if as_of_us is not None:
        _require_divides_window(step)
    if observed_inflow is not None and join_us is None and as_of_us is None:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_JOIN,
            "observed_inflow requires t_join or as_of",
        )
    if (join_us is not None or as_of_us is not None) and not calendar:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_TIMESTAMP,
            "t_join and as_of require timed excess pulses",
        )

    phase = as_of_us if as_of_us is not None else join_us
    ends = _pulse_ends(pulses, step, phase, calendar)
    if join_us is not None:
        ends = [item for item in ends if item[0] > join_us]

    observed = _observed_map(observed_inflow, step, phase if calendar else None)
    response = _direct_runoff(ends, ordinates, step)
    if join_us is not None:
        datum = observed.get(join_us, baseflow)
        support = [instant for instant in response if instant >= join_us]
        if support:
            series = _dense(min(support), max(support), step, response, datum)
        else:
            series = {}
        series[join_us] = datum
        splice: int | None = join_us
        origin = join_us
    else:
        if response:
            series = _dense(min(response), max(response), step, response, baseflow)
        elif ends:
            series = {end: baseflow for end, _excess in ends}
        else:
            series = {}
        splice = as_of_us if observed_inflow is not None and as_of_us is not None else None
        origin = ends[0][0] - step if ends else 0

    if as_of_us is not None:
        origin = as_of_us
    samples = _collect_samples(
        series,
        observed,
        splice=splice,
        as_of_us=as_of_us,
        step=step,
    )
    return _points_from_samples(
        samples,
        origin_us=origin,
        timestep_hours=timestep_hours,
        as_of=as_of,
        calendar=calendar or as_of_us is not None,
    )


def _pulse_ends(
    pulses: Sequence[tuple[str | None, float]],
    step: int,
    phase: int | None,
    calendar: bool,
) -> list[tuple[int, float]]:
    if not calendar:
        return [((index + 1) * step, depth) for index, (_time, depth) in enumerate(pulses)]

    ends: list[tuple[int, float]] = []
    previous: int | None = None
    for index, (stamp, depth) in enumerate(pulses):
        if stamp is None:
            raise ScsUhError(
                ScsUhErrorCode.INVALID_TIMESTAMP,
                "timed excess is missing a timestamp",
            )
        instant = _to_us(_parse_rfc3339(stamp, role=f"excess time at index {index}"))
        if previous is not None:
            delta = instant - previous
            if delta <= 0 or delta % step != 0:
                raise ScsUhError(
                    ScsUhErrorCode.MISALIGNED_TIMES,
                    "excess time at index "
                    f"{index} is not an integer number of timesteps after the previous pulse",
                )
        if phase is not None and (instant - phase) % step != 0:
            raise ScsUhError(
                ScsUhErrorCode.MISALIGNED_TIMES,
                f"excess time at index {index} is not on the timestep grid",
            )
        previous = instant
        ends.append((instant, depth))
    return ends


def generalized_triangle_inflow(
    *,
    peak_discharge_m3s: float,
    as_of: str,
    inflow_at_as_of_m3s: float,
    timestep_hours: float,
    hours_to_peak: float | None = None,
    tau_hours: float | None = None,
    fall_hours: float | None = None,
    observed_inflow: Sequence[tuple[str, float]] | None = None,
) -> InflowHydrograph:
    """Degraded triangular inflow from a commanded peak and a time to peak.

    ``hours_to_peak`` is the rise from ``as_of`` and wins over ``tau_hours``.
    ``tau_hours`` is used only when ``hours_to_peak`` is omitted, and then the
    rise also starts at ``as_of``. Omitted ``fall_hours`` uses a rising to
    recession ratio of 1:2. The series covers at least
    ``[as_of - 12 h, as_of + 12 h]``; unmeasured slots are 0 m³/s.
    """

    peak = _require_discharge(peak_discharge_m3s, ScsUhErrorCode.INVALID_PEAK, "peak_discharge_m3s")
    inflow_now = _require_discharge(
        inflow_at_as_of_m3s,
        ScsUhErrorCode.INVALID_DISCHARGE,
        "inflow_at_as_of_m3s",
    )
    step = _step_us(timestep_hours)
    _require_divides_window(step)
    if hours_to_peak is not None:
        rise = _require_real(hours_to_peak, ScsUhErrorCode.INVALID_HOURS_TO_PEAK, "hours_to_peak")
        rise_name = "hours_to_peak"
    elif tau_hours is not None:
        rise = _require_real(tau_hours, ScsUhErrorCode.INVALID_HOURS_TO_PEAK, "tau_hours")
        rise_name = "tau_hours"
    else:
        raise ScsUhError(
            ScsUhErrorCode.INVALID_HOURS_TO_PEAK,
            "hours_to_peak or tau_hours is required",
        )
    rise_us = _require_multiple(
        rise,
        step,
        code=ScsUhErrorCode.INVALID_HOURS_TO_PEAK,
        name=rise_name,
    )
    if fall_hours is None:
        fall = _TRIANGLE_RECESSION_TO_RISE * rise
    else:
        fall = _require_real(fall_hours, ScsUhErrorCode.INVALID_FALL_HOURS, "fall_hours")
    fall_us = _require_multiple(
        fall,
        step,
        code=ScsUhErrorCode.INVALID_FALL_HOURS,
        name="fall_hours",
    )

    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of"))
    observed = _observed_map(observed_inflow, step, as_of_us)
    observed[as_of_us] = inflow_now
    end_us = as_of_us + rise_us + fall_us
    series: dict[int, float] = {}
    instant = as_of_us
    while instant <= end_us:
        series[instant] = _triangle_inflow(
            instant - as_of_us,
            rise_us,
            fall_us,
            inflow_now,
            peak,
        )
        instant += step

    samples = _collect_samples(
        series,
        observed,
        splice=as_of_us,
        as_of_us=as_of_us,
        step=step,
    )
    return _points_from_samples(
        samples,
        origin_us=as_of_us,
        timestep_hours=timestep_hours,
        as_of=as_of,
        calendar=True,
    )


def _triangle_inflow(
    elapsed_us: int,
    rise_us: int,
    fall_us: int,
    inflow_at_start: float,
    peak: float,
) -> float:
    if elapsed_us <= rise_us:
        return inflow_at_start + (peak - inflow_at_start) * (elapsed_us / rise_us)
    if elapsed_us <= rise_us + fall_us:
        return peak * ((rise_us + fall_us) - elapsed_us) / fall_us
    return 0.0

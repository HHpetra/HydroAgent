"""Muskingum and kinematic variable-lag channel routing.

Discharge is in m³/s, travel time in hours, and clock times are RFC3339
strings. Nothing in this module reads the wall clock. A scenario ``as_of`` is
an argument.

Muskingum routes a reach by the HEC-HMS storage equation. When the timestep is
shorter than ``2 K x``, the reach is split into ``ceil(K / Δt)`` subreaches and
routed in series. Variable lag maps each station discharge to a dam-site
arrival; a larger discharge has a shorter lag. ``blend_station_and_forecast``
keeps measured dam inflow through ``as_of``, uses the lagged station hydrograph
until the last arrival ``t_join``, and after that adds SCS-CN / SCS-TUH runoff
from rainfall whose interval ends after ``t_join``. It does not extend the
station hydrograph by slope. With ``as_of`` set, the blended series covers at
least ``[as_of - 12 h, as_of + 12 h]``; slots with no ordinate are 0 m³/s.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from enum import StrEnum

from hydro_agent.core.scs_cn import ScsCnParameters, TimedRainfall
from hydro_agent.core.scs_uh import InflowHydrograph, InflowPoint, Watershed, route_rainfall

_US_PER_HOUR = 3_600_000_000
_WINDOW_HALF_HOURS = 12.0
_MAX_STEPS = 1_000_000
_MAX_SUBREACHES = 100_000
_STABLE_TOLERANCE = 1e-8

_RFC3339 = re.compile(
    r"^\d{4}-\d{2}-\d{2}[Tt]\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:[Zz]|[+-]\d{2}:\d{2})$"
)
_EPOCH = datetime(1970, 1, 1, tzinfo=UTC)


class ChannelRoutingErrorCode(StrEnum):
    """Stable identifier for a channel-routing domain failure."""

    INVALID_DISCHARGE = "invalid_discharge"
    INVALID_TIMESTAMP = "invalid_timestamp"
    UNORDERED_TIMESTAMPS = "unordered_timestamps"
    EMPTY_SERIES = "empty_series"
    INVALID_TIMESTEP = "invalid_timestep"
    MISALIGNED_TIMES = "misaligned_times"
    INVALID_K = "invalid_k"
    INVALID_X = "invalid_x"
    UNSTABLE_TIMESTEP = "unstable_timestep"
    INVALID_LAG = "invalid_lag"
    INVALID_BETA = "invalid_beta"
    INVALID_REFERENCE_DISCHARGE = "invalid_reference_discharge"
    INVALID_AREA = "invalid_area"
    MISSING_FUTURE_RAINFALL = "missing_future_rainfall"
    STATION_NOT_AT_AS_OF = "station_not_at_as_of"
    MISALIGNED_JOIN = "misaligned_join"
    INVALID_PARAMETERS = "invalid_parameters"


class ChannelRoutingError(Exception):
    """Illegal channel-routing input.

    ``code`` identifies the failure. ``message`` describes the offending value.
    """

    def __init__(self, code: ChannelRoutingErrorCode, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(f"{code.value}: {message}")


def _require_real(value: object, code: ChannelRoutingErrorCode, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ChannelRoutingError(code, f"{name} must be a real number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise ChannelRoutingError(code, f"{name} must be finite, got {value!r}")
    return number


def _require_discharge(value: object, name: str) -> float:
    number = _require_real(value, ChannelRoutingErrorCode.INVALID_DISCHARGE, name)
    if number < 0.0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_DISCHARGE,
            f"{name} must be >= 0 m³/s, got {value!r}",
        )
    return number


def _parse_rfc3339(value: object, *, role: str) -> datetime:
    if not isinstance(value, str) or _RFC3339.fullmatch(value) is None:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTAMP,
            f"{role} must be an RFC3339 timestamp, got {value!r}",
        )
    normalized = f"{value[:-1]}+00:00" if value[-1] in "Zz" else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError as exc:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTAMP,
            f"{role} is not a valid RFC3339 timestamp: {value}",
        ) from exc
    if parsed.tzinfo is None:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTAMP,
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


def _hours_to_us(hours: float) -> int:
    return round(hours * _US_PER_HOUR)


def _step_us(timestep_hours: float) -> int:
    number = _require_real(
        timestep_hours,
        ChannelRoutingErrorCode.INVALID_TIMESTEP,
        "timestep_hours",
    )
    if number <= 0.0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be > 0, got {timestep_hours!r}",
        )
    step = round(number * _US_PER_HOUR)
    if step <= 0 or abs(number * _US_PER_HOUR - step) > 1e-3:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTEP,
            f"timestep_hours must be an exact microsecond count, got {timestep_hours!r}",
        )
    return step


def _require_divides_window(step: int) -> None:
    window = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    if window % step != 0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTEP,
            "timestep_hours must divide the 12 h scenario window",
        )


def _ceil_tol(value: float) -> int:
    nearest = round(value)
    if abs(value - nearest) <= _STABLE_TOLERANCE:
        return int(nearest)
    return math.ceil(value - 1e-12)


def _floor_tol(value: float) -> int:
    nearest = round(value)
    if abs(value - nearest) <= _STABLE_TOLERANCE:
        return int(nearest)
    return math.floor(value + 1e-12)


def _parse_series(series: object, name: str) -> list[tuple[int, float, str]]:
    if isinstance(series, str | bytes) or not isinstance(series, Sequence):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_DISCHARGE,
            f"{name} must be a sequence of (time, discharge) pairs",
        )
    if len(series) == 0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.EMPTY_SERIES,
            f"{name} is empty",
        )
    if len(series) > _MAX_STEPS:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTEP,
            f"{name} has {len(series)} points; the limit is {_MAX_STEPS}",
        )
    parsed: list[tuple[int, float, str]] = []
    previous: int | None = None
    for index, item in enumerate(series):
        if not isinstance(item, tuple) or len(item) != 2:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_DISCHARGE,
                f"{name} at index {index} must be (time, discharge)",
            )
        stamp, raw_discharge = item
        instant = _to_us(_parse_rfc3339(stamp, role=f"{name} time at index {index}"))
        discharge = _require_discharge(raw_discharge, f"{name} at index {index}")
        if previous is not None and instant <= previous:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.UNORDERED_TIMESTAMPS,
                f"{name} time at index {index} must be strictly after the previous time",
            )
        previous = instant
        parsed.append((instant, discharge, stamp))
    return parsed


def _hydrograph(
    samples: Sequence[tuple[int, float, str]],
    *,
    origin_us: int,
    timestep_hours: float,
    as_of: str | None,
) -> InflowHydrograph:
    if not samples:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.EMPTY_SERIES,
            "inflow series is empty",
        )
    points = tuple(
        InflowPoint(
            time=stamp,
            hours_from_origin=(instant - origin_us) / _US_PER_HOUR,
            discharge_m3s=discharge,
        )
        for instant, discharge, stamp in samples
    )
    peak = points[0]
    for point in points[1:]:
        if point.discharge_m3s > peak.discharge_m3s:
            peak = point
    volume = math.fsum(point.discharge_m3s for point in points) * 3600.0 * timestep_hours
    return InflowHydrograph(
        points=points,
        peak_discharge_m3s=peak.discharge_m3s,
        peak_time=peak.time,
        peak_hours_from_origin=peak.hours_from_origin,
        volume_m3=volume,
        timestep_hours=timestep_hours,
        as_of=as_of,
    )


@dataclass(frozen=True, slots=True)
class AreaScale:
    """Dam-site scale κ = (F_dam / F_sta) ** exponent, or an explicit ratio.

    An explicit ``ratio`` is κ itself. Otherwise both areas are required.
    ``exponent`` defaults to 1. Omitted scaling leaves κ = 1.
    """

    station_area_km2: float | None = None
    dam_area_km2: float | None = None
    exponent: float = 1.0
    ratio: float | None = None

    def __post_init__(self) -> None:
        exponent = _require_real(self.exponent, ChannelRoutingErrorCode.INVALID_AREA, "exponent")
        object.__setattr__(self, "exponent", exponent)
        if self.ratio is not None:
            ratio = _require_real(self.ratio, ChannelRoutingErrorCode.INVALID_AREA, "ratio")
            if ratio < 0.0:
                raise ChannelRoutingError(
                    ChannelRoutingErrorCode.INVALID_AREA,
                    f"ratio must be >= 0, got {self.ratio!r}",
                )
            object.__setattr__(self, "ratio", ratio)
            return
        if self.station_area_km2 is None and self.dam_area_km2 is None:
            return
        station = _require_real(
            self.station_area_km2,
            ChannelRoutingErrorCode.INVALID_AREA,
            "station_area_km2",
        )
        dam = _require_real(
            self.dam_area_km2,
            ChannelRoutingErrorCode.INVALID_AREA,
            "dam_area_km2",
        )
        if station <= 0.0 or dam <= 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_AREA,
                "station_area_km2 and dam_area_km2 must be > 0",
            )
        object.__setattr__(self, "station_area_km2", station)
        object.__setattr__(self, "dam_area_km2", dam)


def _area_ratio(scale: AreaScale | None) -> float:
    if scale is None or (
        scale.ratio is None and scale.station_area_km2 is None and scale.dam_area_km2 is None
    ):
        return 1.0
    if scale.ratio is not None:
        return scale.ratio
    assert scale.station_area_km2 is not None and scale.dam_area_km2 is not None
    return math.pow(scale.dam_area_km2 / scale.station_area_km2, scale.exponent)


@dataclass(frozen=True, slots=True)
class MuskingumCoefficients:
    """Subreach coefficients actually used for the cascade.

    ``c0 + c1 + c2`` is 1. ``k_sub_hours`` is ``K / n_subreaches``.
    """

    c0: float
    c1: float
    c2: float
    n_subreaches: int
    k_sub_hours: float
    timestep_hours: float

    def __post_init__(self) -> None:
        total = self.c0 + self.c1 + self.c2
        if abs(total - 1.0) > 1e-9:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.UNSTABLE_TIMESTEP,
                f"Muskingum coefficients sum to {total}, not 1",
            )


@dataclass(frozen=True, slots=True)
class MuskingumReach:
    """One channel reach: travel-time constant ``K`` (h) and weighting ``x``."""

    k_hours: float
    x: float
    area_scale: AreaScale | None = None

    def __post_init__(self) -> None:
        k_hours = _require_real(self.k_hours, ChannelRoutingErrorCode.INVALID_K, "k_hours")
        if k_hours <= 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_K,
                f"k_hours must be > 0, got {self.k_hours!r}",
            )
        x = _require_real(self.x, ChannelRoutingErrorCode.INVALID_X, "x")
        if x < 0.0 or x > 0.5:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_X,
                f"x must be in [0, 0.5], got {self.x!r}",
            )
        object.__setattr__(self, "k_hours", k_hours)
        object.__setattr__(self, "x", x)
        if self.area_scale is not None and not isinstance(self.area_scale, AreaScale):
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_AREA,
                f"area_scale must be AreaScale, got {self.area_scale!r}",
            )


@dataclass(frozen=True, slots=True)
class MuskingumResult:
    """Routed dam-site outflow, including the subreach coefficients."""

    hydrograph: InflowHydrograph
    coefficients: MuskingumCoefficients
    area_ratio: float


def _subreach_count(k_hours: float, x: float, timestep_hours: float) -> int:
    """HEC-HMS split: N = ceil(K / Δt) when Δt < 2 K x and that N is stable.

    N must satisfy ``2 K_sub x <= Δt <= 2 K_sub (1 - x)``. If the ceiling count
    falls outside that range, use the nearest integer that does. Splitting
    cannot repair a timestep that is already longer than ``2 K (1 - x)``.
    """

    if x == 0.0:
        n_min = 1
    else:
        n_min = max(1, _ceil_tol((2.0 * k_hours * x) / timestep_hours))
    n_max = _floor_tol((2.0 * k_hours * (1.0 - x)) / timestep_hours)
    if n_max < n_min:
        upper = 2.0 * k_hours * (1.0 - x)
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.UNSTABLE_TIMESTEP,
            "no Muskingum subreach count satisfies "
            f"2 K_sub x <= Δt <= 2 K_sub (1-x) for K={k_hours} h, x={x}, "
            f"Δt={timestep_hours} h (single-reach upper bound is {upper} h)",
        )
    if timestep_hours + _STABLE_TOLERANCE < 2.0 * k_hours * x:
        preferred = max(1, _ceil_tol(k_hours / timestep_hours))
        if preferred > _MAX_SUBREACHES:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.UNSTABLE_TIMESTEP,
                f"Muskingum subreach count {preferred} exceeds {_MAX_SUBREACHES}",
            )
        if n_min <= preferred <= n_max:
            return preferred
        return min(range(n_min, n_max + 1), key=lambda count: (abs(count - preferred), count))
    return 1


def _coefficients(
    k_sub_hours: float,
    x: float,
    timestep_hours: float,
) -> tuple[float, float, float]:
    half = 0.5 * timestep_hours
    denominator = k_sub_hours * (1.0 - x) + half
    if denominator <= 0.0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.UNSTABLE_TIMESTEP,
            "Muskingum denominator must be > 0",
        )
    c0 = (half - k_sub_hours * x) / denominator
    c1 = (half + k_sub_hours * x) / denominator
    c2 = (k_sub_hours * (1.0 - x) - half) / denominator
    return c0, c1, c2


def _one_reach(inflow: Sequence[float], c0: float, c1: float, c2: float) -> list[float]:
    """One Muskingum step with a steady start ``O0 = I0``."""

    outflow = [inflow[0]]
    previous_in = inflow[0]
    previous_out = inflow[0]
    for incoming in inflow[1:]:
        next_out = c0 * incoming + c1 * previous_in + c2 * previous_out
        outflow.append(next_out)
        previous_in = incoming
        previous_out = next_out
    return outflow


def route_muskingum(
    inflow: Sequence[tuple[str, float]],
    reach: MuskingumReach,
    *,
    timestep_hours: float,
) -> MuskingumResult:
    """Route ``inflow`` through a Muskingum reach and scale by the area ratio.

    The outflow is at the same instants as the inflow. ``O0 = I0``. When
    ``x = 0.5`` and ``Δt = K`` (or each subreach has ``K_sub = Δt``), the wave
    is translated by ``K`` with no peak attenuation. κ = 1 conserves volume
    once storage has returned to the initial steady state.
    """

    if not isinstance(reach, MuskingumReach):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_PARAMETERS,
            f"reach must be MuskingumReach, got {reach!r}",
        )
    parsed = _parse_series(inflow, "inflow")
    step = _step_us(timestep_hours)
    for index in range(1, len(parsed)):
        if parsed[index][0] - parsed[index - 1][0] != step:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.MISALIGNED_TIMES,
                f"inflow time at index {index} is not {timestep_hours} h after the previous time",
            )
    n_sub = _subreach_count(reach.k_hours, reach.x, timestep_hours)
    k_sub = reach.k_hours / n_sub
    c0, c1, c2 = _coefficients(k_sub, reach.x, timestep_hours)
    routed = [discharge for _instant, discharge, _stamp in parsed]
    for _ in range(n_sub):
        routed = _one_reach(routed, c0, c1, c2)
    ratio = _area_ratio(reach.area_scale)
    samples = [
        (instant, ratio * discharge, stamp)
        for (instant, _original, stamp), discharge in zip(parsed, routed, strict=True)
    ]
    hydrograph = _hydrograph(
        samples,
        origin_us=parsed[0][0],
        timestep_hours=timestep_hours,
        as_of=None,
    )
    return MuskingumResult(
        hydrograph=hydrograph,
        coefficients=MuskingumCoefficients(
            c0=c0,
            c1=c1,
            c2=c2,
            n_subreaches=n_sub,
            k_sub_hours=k_sub,
            timestep_hours=timestep_hours,
        ),
        area_ratio=ratio,
    )


@dataclass(frozen=True, slots=True)
class VariableLag:
    """Kinematic lag τ(Q) = τ_ref * (max(Q, Q_min) / Q_ref) ** (-β).

    ``beta = 0`` does not change the lag with discharge. The result is clamped to
    ``[tau_min_hours, tau_max_hours]`` when those bounds are set. A larger
    discharge has a shorter lag when ``beta > 0``.
    """

    tau_ref_hours: float
    q_ref_m3s: float = 1.0
    beta: float = 0.0
    q_min_m3s: float = 0.0
    tau_min_hours: float | None = None
    tau_max_hours: float | None = None
    area_scale: AreaScale | None = None

    def __post_init__(self) -> None:
        tau_ref = _require_real(
            self.tau_ref_hours,
            ChannelRoutingErrorCode.INVALID_LAG,
            "tau_ref_hours",
        )
        if tau_ref < 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_LAG,
                f"tau_ref_hours must be >= 0, got {self.tau_ref_hours!r}",
            )
        q_ref = _require_real(
            self.q_ref_m3s,
            ChannelRoutingErrorCode.INVALID_REFERENCE_DISCHARGE,
            "q_ref_m3s",
        )
        if q_ref <= 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_REFERENCE_DISCHARGE,
                f"q_ref_m3s must be > 0, got {self.q_ref_m3s!r}",
            )
        beta = _require_real(self.beta, ChannelRoutingErrorCode.INVALID_BETA, "beta")
        if beta < 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_BETA,
                f"beta must be >= 0 so higher discharge travels faster, got {self.beta!r}",
            )
        q_min = _require_real(
            self.q_min_m3s,
            ChannelRoutingErrorCode.INVALID_REFERENCE_DISCHARGE,
            "q_min_m3s",
        )
        if q_min < 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_REFERENCE_DISCHARGE,
                f"q_min_m3s must be >= 0, got {self.q_min_m3s!r}",
            )
        tau_min = _optional_lag(self.tau_min_hours, "tau_min_hours")
        tau_max = _optional_lag(self.tau_max_hours, "tau_max_hours")
        if tau_min is not None and tau_max is not None and tau_min > tau_max:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_LAG,
                f"tau_min_hours ({tau_min}) must be <= tau_max_hours ({tau_max})",
            )
        object.__setattr__(self, "tau_ref_hours", tau_ref)
        object.__setattr__(self, "q_ref_m3s", q_ref)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "q_min_m3s", q_min)
        object.__setattr__(self, "tau_min_hours", tau_min)
        object.__setattr__(self, "tau_max_hours", tau_max)
        if self.area_scale is not None and not isinstance(self.area_scale, AreaScale):
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_AREA,
                f"area_scale must be AreaScale, got {self.area_scale!r}",
            )


def _optional_lag(value: object, name: str) -> float | None:
    if value is None:
        return None
    number = _require_real(value, ChannelRoutingErrorCode.INVALID_LAG, name)
    if number < 0.0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_LAG,
            f"{name} must be >= 0, got {value!r}",
        )
    return number


def _lag_hours(lag: VariableLag, discharge: float) -> float:
    if lag.beta == 0.0:
        raw = lag.tau_ref_hours
    else:
        discharge_for_lag = max(discharge, lag.q_min_m3s)
        if discharge_for_lag <= 0.0:
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_LAG,
                "variable lag needs q_min_m3s > 0 when discharge is 0",
            )
        raw = lag.tau_ref_hours * math.pow(discharge_for_lag / lag.q_ref_m3s, -lag.beta)
    if lag.tau_min_hours is not None:
        raw = max(raw, lag.tau_min_hours)
    if lag.tau_max_hours is not None:
        raw = min(raw, lag.tau_max_hours)
    return raw


@dataclass(frozen=True, slots=True)
class LaggedArrival:
    """One station ordinate mapped to the dam site."""

    station_time: str
    arrival_time: str
    station_discharge_m3s: float
    discharge_m3s: float
    lag_hours: float


@dataclass(frozen=True, slots=True)
class VariableLagResult:
    """Dam-site arrivals ordered by arrival time, plus the scenario join."""

    arrivals: tuple[LaggedArrival, ...]
    as_of: str
    t_join: str
    area_ratio: float
    peak_discharge_m3s: float
    peak_arrival_time: str


def _monotone_envelope(mapped: Sequence[tuple[int, float, LaggedArrival]]) -> list[LaggedArrival]:
    """Keep the larger discharge when a fast wave arrives before a slow one.

    Summing the overlapped water would raise the dam-site peak above κ times
    the station peak. The envelope keeps that peak.
    """

    ordered = sorted(mapped, key=lambda item: (item[0], -item[1]))
    envelope: list[LaggedArrival] = []
    last_arrival: int | None = None
    for arrival_us, discharge, lagged in ordered:
        if last_arrival is not None and arrival_us <= last_arrival:
            if envelope and discharge > envelope[-1].discharge_m3s:
                envelope[-1] = lagged
            continue
        envelope.append(lagged)
        last_arrival = arrival_us
    return envelope


def _mapped_arrivals(
    station: Sequence[tuple[int, float, str]],
    lag: VariableLag,
    ratio: float,
) -> list[tuple[int, float, LaggedArrival]]:
    mapped: list[tuple[int, float, LaggedArrival]] = []
    for instant, discharge, stamp in station:
        lag_hours = _lag_hours(lag, discharge)
        arrival_us = instant + _hours_to_us(lag_hours)
        scaled = ratio * discharge
        mapped.append(
            (
                arrival_us,
                scaled,
                LaggedArrival(
                    station_time=stamp,
                    arrival_time=_format_us(arrival_us),
                    station_discharge_m3s=discharge,
                    discharge_m3s=scaled,
                    lag_hours=lag_hours,
                ),
            )
        )
    return mapped


def route_variable_lag(
    station_flow: Sequence[tuple[str, float]],
    lag: VariableLag,
    *,
    as_of: str,
) -> VariableLagResult:
    """Map station discharge to dam-site arrivals and report ``t_join``.

    ``t_join = as_of + τ(Q_sta(as_of))``. Each arrival discharge is κ times the
    station discharge, so the dam-site peak equals κ times the station peak.
    """

    if not isinstance(lag, VariableLag):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_PARAMETERS,
            f"lag must be VariableLag, got {lag!r}",
        )
    parsed = _parse_series(station_flow, "station_flow")
    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of"))
    if parsed[-1][0] > as_of_us:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_TIMESTAMP,
            "station_flow contains an observation after as_of "
            f"({as_of}); later station data is not part of this scenario",
        )
    at_as_of = [discharge for instant, discharge, _stamp in parsed if instant == as_of_us]
    if not at_as_of:
        latest = _format_us(parsed[-1][0])
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.STATION_NOT_AT_AS_OF,
            f"station_flow must include an observation at as_of ({as_of}); "
            f"the latest observation is {latest}",
        )
    ratio = _area_ratio(lag.area_scale)
    envelope = _monotone_envelope(_mapped_arrivals(parsed, lag, ratio))
    if not envelope:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.EMPTY_SERIES,
            "station_flow produced no arrivals",
        )
    peak = envelope[0]
    for arrival in envelope[1:]:
        if arrival.discharge_m3s > peak.discharge_m3s:
            peak = arrival
    join_lag = _lag_hours(lag, at_as_of[0])
    return VariableLagResult(
        arrivals=tuple(envelope),
        as_of=as_of,
        t_join=_format_us(as_of_us + _hours_to_us(join_lag)),
        area_ratio=ratio,
        peak_discharge_m3s=peak.discharge_m3s,
        peak_arrival_time=peak.arrival_time,
    )


@dataclass(frozen=True, slots=True)
class BlendedInflow:
    """Dam-site inflow spliced from observation, station lag, and remainder rain."""

    hydrograph: InflowHydrograph
    as_of: str
    t_join: str
    area_ratio: float
    join_discharge_m3s: float
    routed_peak_discharge_m3s: float
    routed_peak_time: str


def _interp(anchors: Sequence[tuple[int, float]], instant: int) -> float:
    if not anchors or instant < anchors[0][0]:
        return 0.0
    if instant >= anchors[-1][0]:
        return anchors[-1][1]
    lo = 0
    hi = len(anchors) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if anchors[mid][0] <= instant:
            lo = mid
        else:
            hi = mid
    left_t, left_q = anchors[lo]
    right_t, right_q = anchors[hi]
    if right_t == left_t:
        return right_q
    fraction = (instant - left_t) / (right_t - left_t)
    return left_q + (right_q - left_q) * fraction


def _grid(start: int, end: int, step: int) -> list[int]:
    if end < start:
        return []
    count = (end - start) // step
    return [start + index * step for index in range(count + 1)]


def _require_future_rainfall(rainfall: object, as_of_us: int, as_of: str) -> list[TimedRainfall]:
    if isinstance(rainfall, str | bytes) or not isinstance(rainfall, Sequence):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.MISSING_FUTURE_RAINFALL,
            "station_and_forecast blending needs a timed rainfall series",
        )
    if len(rainfall) == 0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.MISSING_FUTURE_RAINFALL,
            f"station_and_forecast blending needs rainfall after as_of ({as_of}); "
            "the hyetograph is empty",
        )
    points: list[TimedRainfall] = []
    seen_future = False
    for index, item in enumerate(rainfall):
        if not isinstance(item, TimedRainfall):
            raise ChannelRoutingError(
                ChannelRoutingErrorCode.INVALID_PARAMETERS,
                f"rainfall at index {index} must be TimedRainfall",
            )
        instant = _to_us(_parse_rfc3339(item.time, role=f"rainfall time at index {index}"))
        if instant > as_of_us:
            seen_future = True
        points.append(item)
    if not seen_future:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.MISSING_FUTURE_RAINFALL,
            "station_and_forecast blending needs rainfall after as_of "
            f"({as_of}); the hyetograph only contains rainfall at or before as_of",
        )
    return points


def _anchors_through_join(
    envelope: Sequence[LaggedArrival],
    t_join_us: int,
    q_join: float,
) -> list[tuple[int, float]]:
    anchors: list[tuple[int, float]] = []
    for arrival in envelope:
        arrival_us = _to_us(_parse_rfc3339(arrival.arrival_time, role="arrival time"))
        if arrival_us <= t_join_us:
            anchors.append((arrival_us, arrival.discharge_m3s))
    if anchors and anchors[-1][0] == t_join_us:
        anchors[-1] = (t_join_us, q_join)
    else:
        anchors.append((t_join_us, q_join))
    return anchors


def blend_station_and_forecast(
    station_flow: Sequence[tuple[str, float]],
    rainfall: Sequence[TimedRainfall],
    loss: ScsCnParameters,
    watershed: Watershed,
    lag: VariableLag,
    *,
    as_of: str,
    timestep_hours: float,
    observed_inflow: Sequence[tuple[str, float]] | None = None,
) -> BlendedInflow:
    """Splice observed inflow, lagged station flow, and remainder-rain runoff.

    For ``t <= as_of`` the series is the measured dam inflow, including a
    measured 0 when the flood has not arrived. For ``as_of < t <= t_join`` it
    is the variable-lag routing only. For ``t > t_join`` it is
    ``Q_arr(t_join) + max(0, Q_rem(t) - Q_rem(t_join))``, where ``Q_rem`` is the
    SCS-TUH convolution of excess pulses that end after ``t_join``. Initial
    abstraction comes from the whole storm. Hours after ``t_join`` and before
    the first remainder ordinate stay at the join discharge, because ``Q_rem``
    is still 0 there. After the remainder hydrograph ends, a scenario-window
    hour with no ordinate stays 0. ``t_join`` must fall on the timestep grid
    so the SCS-TUH splice can see it. Rainfall that ends at or before
    ``as_of`` is not a forecast: the call fails with
    ``missing_future_rainfall`` instead of inventing a slope.
    """

    if not isinstance(loss, ScsCnParameters):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_PARAMETERS,
            f"loss must be ScsCnParameters, got {loss!r}",
        )
    if not isinstance(watershed, Watershed):
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.INVALID_PARAMETERS,
            f"watershed must be Watershed, got {watershed!r}",
        )
    routed = route_variable_lag(station_flow, lag, as_of=as_of)
    as_of_us = _to_us(_parse_rfc3339(as_of, role="as_of"))
    t_join_us = _to_us(_parse_rfc3339(routed.t_join, role="t_join"))
    step = _step_us(timestep_hours)
    _require_divides_window(step)
    if (t_join_us - as_of_us) % step != 0:
        raise ChannelRoutingError(
            ChannelRoutingErrorCode.MISALIGNED_JOIN,
            f"t_join ({routed.t_join}) is not on the {timestep_hours} h grid from as_of",
        )
    hyetograph = _require_future_rainfall(rainfall, as_of_us, as_of)

    station = _parse_series(station_flow, "station_flow")
    q_station = next(discharge for instant, discharge, _stamp in station if instant == as_of_us)
    q_join = routed.area_ratio * q_station
    anchors = _anchors_through_join(routed.arrivals, t_join_us, q_join)

    # t <= as_of stays on the measured inflow. The join discharge is written
    # only when it falls strictly after as_of; a zero lag must not replace it.
    values: dict[int, float] = {}
    if t_join_us > as_of_us:
        for instant in _grid(as_of_us + step, t_join_us, step):
            values[instant] = _interp(anchors, instant)

    remainder = route_rainfall(
        hyetograph,
        loss,
        watershed,
        timestep_hours=timestep_hours,
        t_join=routed.t_join,
        observed_inflow=((routed.t_join, q_join),),
    )
    remainder_times: list[int] = []
    for point in remainder.points:
        if point.time is None:
            continue
        instant = _to_us(_parse_rfc3339(point.time, role="remainder inflow time"))
        if instant <= t_join_us:
            continue
        increment = point.discharge_m3s - q_join
        if increment < 0.0:
            increment = 0.0
        values[instant] = q_join + increment
        remainder_times.append(instant)
    if remainder_times:
        # Q_rem is 0 until the first convolved ordinate, so the splice holds
        # the join discharge instead of dropping to an empty-slot 0.
        first_remainder = min(remainder_times)
        for instant in _grid(t_join_us + step, first_remainder - step, step):
            values.setdefault(instant, q_join)

    if observed_inflow:
        observed = _parse_series(observed_inflow, "observed_inflow")
        for index, (instant, discharge, _stamp) in enumerate(observed):
            if instant > as_of_us:
                raise ChannelRoutingError(
                    ChannelRoutingErrorCode.INVALID_TIMESTAMP,
                    f"observed_inflow at index {index} is after as_of; "
                    "dam-site observations after the scenario time are not used",
                )
            if (instant - as_of_us) % step != 0:
                raise ChannelRoutingError(
                    ChannelRoutingErrorCode.MISALIGNED_TIMES,
                    f"observed_inflow at index {index} is not on the timestep grid",
                )
            values[instant] = discharge

    half = round(_WINDOW_HALF_HOURS * _US_PER_HOUR)
    for instant in _grid(as_of_us - half, as_of_us + half, step):
        values.setdefault(instant, 0.0)

    samples = [
        (instant, values[instant], _format_us(instant)) for instant in sorted(values)
    ]
    hydrograph = _hydrograph(
        samples,
        origin_us=as_of_us,
        timestep_hours=timestep_hours,
        as_of=as_of,
    )
    return BlendedInflow(
        hydrograph=hydrograph,
        as_of=as_of,
        t_join=routed.t_join,
        area_ratio=routed.area_ratio,
        join_discharge_m3s=q_join,
        routed_peak_discharge_m3s=routed.peak_discharge_m3s,
        routed_peak_time=routed.peak_arrival_time,
    )

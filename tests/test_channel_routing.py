"""Channel routing: Muskingum, kinematic variable lag, and station/forecast splicing."""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import pytest

from hydro_agent.core import (
    AreaScale,
    ChannelRoutingError,
    ChannelRoutingErrorCode,
    InflowHydrograph,
    MuskingumReach,
    ScsCnParameters,
    TimedRainfall,
    VariableLag,
    Watershed,
    blend_station_and_forecast,
    route_muskingum,
    route_rainfall,
    route_variable_lag,
)

_AS_OF = "2020-07-01T03:00:00Z"
_T_JOIN = "2020-07-01T06:00:00Z"
_DAY = datetime(2020, 7, 1, tzinfo=UTC)


def _stamp(instant: datetime) -> str:
    return instant.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hourly(
    discharges: tuple[float, ...],
    *,
    start: datetime | None = None,
) -> list[tuple[str, float]]:
    origin = _DAY if start is None else start
    return [
        (_stamp(origin + timedelta(hours=index)), discharge)
        for index, discharge in enumerate(discharges)
    ]


def _rain(depths_by_hour: dict[int, float]) -> list[TimedRainfall]:
    return [
        TimedRainfall(_stamp(_DAY + timedelta(hours=hour)), depths_by_hour[hour])
        for hour in sorted(depths_by_hour)
    ]


def _forecast_rain(*, tail: float = 8.0, early: float = 6.0) -> list[TimedRainfall]:
    depths = {
        1: 12.0,
        2: 12.0,
        3: 10.0,
        4: early,
        5: early,
        6: 6.0,
        7: 8.0,
        8: tail,
        9: tail,
        10: tail,
        11: tail,
        12: tail,
    }
    return _rain(depths)


def _station_flow() -> list[tuple[str, float]]:
    return _hourly((0.0, 30.0, 80.0, 120.0))


def _observed_still_dry() -> list[tuple[str, float]]:
    return _hourly((0.0, 0.0, 0.0, 0.0))


def _three_hour_lag() -> VariableLag:
    return VariableLag(tau_ref_hours=3.0, q_ref_m3s=100.0, beta=0.0)


def _volume(discharges: tuple[float, ...] | list[float], timestep_hours: float) -> float:
    return sum(discharges) * 3600.0 * timestep_hours


def _discharges(hydrograph: InflowHydrograph) -> tuple[float, ...]:
    return tuple(point.discharge_m3s for point in hydrograph.points)


def test_muskingum_with_x_half_and_timestep_equal_to_k_translates_without_attenuation() -> None:
    inflow = _hourly((0.0, 0.0, 5.0, 12.0, 7.0, 0.0, 0.0))
    routed = route_muskingum(
        inflow,
        MuskingumReach(k_hours=1.0, x=0.5),
        timestep_hours=1.0,
    )

    assert routed.coefficients.c0 == 0.0
    assert routed.coefficients.c1 == 1.0
    assert routed.coefficients.c2 == 0.0
    assert routed.coefficients.c0 + routed.coefficients.c1 + routed.coefficients.c2 == 1.0
    assert routed.coefficients.n_subreaches == 1
    assert _discharges(routed.hydrograph) == (0.0, 0.0, 0.0, 5.0, 12.0, 7.0, 0.0)
    assert routed.hydrograph.peak_discharge_m3s == 12.0
    assert routed.hydrograph.peak_time == "2020-07-01T04:00:00Z"
    inflow_volume = _volume((0, 0, 5, 12, 7, 0, 0), 1)
    assert routed.hydrograph.volume_m3 == pytest.approx(inflow_volume, rel=1e-6)


def test_muskingum_subreaches_translate_the_wave_by_the_full_travel_time() -> None:
    inflow = _hourly((0.0, 0.0, 4.0, 9.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    routed = route_muskingum(
        inflow,
        MuskingumReach(k_hours=4.0, x=0.5),
        timestep_hours=1.0,
    )

    assert routed.coefficients.n_subreaches == 4
    assert routed.coefficients.k_sub_hours == 1.0
    assert routed.coefficients.c0 == 0.0
    assert routed.coefficients.c1 == 1.0
    assert routed.coefficients.c2 == 0.0
    assert _discharges(routed.hydrograph) == (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 9.0, 4.0, 0.0)
    assert routed.hydrograph.peak_discharge_m3s == 9.0
    assert routed.hydrograph.peak_time == "2020-07-01T07:00:00Z"
    assert routed.hydrograph.volume_m3 == pytest.approx(
        _volume((0, 0, 4, 9, 4, 0, 0, 0, 0, 0), 1),
        rel=1e-6,
    )


def test_muskingum_attenuation_cuts_the_peak_delays_it_and_conserves_volume() -> None:
    rising = (0.0, 2.0, 8.0, 18.0, 12.0, 6.0, 3.0, 1.0)
    inflow = _hourly(rising + (0.0,) * 60)
    routed = route_muskingum(
        inflow,
        MuskingumReach(k_hours=2.0, x=0.2),
        timestep_hours=1.0,
    )
    inflow_peak_time = "2020-07-01T03:00:00Z"

    assert routed.coefficients.n_subreaches == 1
    coefficient_sum = routed.coefficients.c0 + routed.coefficients.c1 + routed.coefficients.c2
    assert coefficient_sum == pytest.approx(1.0)
    assert routed.hydrograph.peak_discharge_m3s < 18.0
    assert routed.hydrograph.peak_time is not None
    assert routed.hydrograph.peak_time > inflow_peak_time
    inflow_volume = _volume([discharge for _time, discharge in inflow], 1)
    assert routed.hydrograph.volume_m3 == pytest.approx(inflow_volume, rel=1e-6)
    assert routed.area_ratio == 1.0


def test_area_ratio_multiplies_muskingum_outflow() -> None:
    inflow = _hourly((0.0, 0.0, 5.0, 12.0, 7.0, 0.0, 0.0))
    scaled = route_muskingum(
        inflow,
        MuskingumReach(
            k_hours=1.0,
            x=0.5,
            area_scale=AreaScale(station_area_km2=40.0, dam_area_km2=80.0),
        ),
        timestep_hours=1.0,
    )
    squared = route_muskingum(
        inflow,
        MuskingumReach(
            k_hours=1.0,
            x=0.5,
            area_scale=AreaScale(station_area_km2=40.0, dam_area_km2=80.0, exponent=2.0),
        ),
        timestep_hours=1.0,
    )

    assert scaled.area_ratio == 2.0
    assert _discharges(scaled.hydrograph) == (0.0, 0.0, 0.0, 10.0, 24.0, 14.0, 0.0)
    assert squared.area_ratio == 4.0
    assert squared.hydrograph.peak_discharge_m3s == 48.0


def test_variable_lag_is_shorter_at_high_discharge_and_peak_scales_with_area() -> None:
    station = [
        ("2020-07-01T00:00:00Z", 50.0),
        ("2020-07-01T01:00:00Z", 200.0),
    ]
    routed = route_variable_lag(
        station,
        VariableLag(
            tau_ref_hours=4.0,
            q_ref_m3s=100.0,
            beta=1.0,
            area_scale=AreaScale(station_area_km2=25.0, dam_area_km2=100.0),
        ),
        as_of="2020-07-01T01:00:00Z",
    )
    squared = route_variable_lag(
        station,
        VariableLag(
            tau_ref_hours=4.0,
            q_ref_m3s=100.0,
            beta=1.0,
            area_scale=AreaScale(station_area_km2=25.0, dam_area_km2=100.0, exponent=2.0),
        ),
        as_of="2020-07-01T01:00:00Z",
    )
    by_station = {arrival.station_time: arrival for arrival in routed.arrivals}

    high = by_station["2020-07-01T01:00:00Z"]
    low = by_station["2020-07-01T00:00:00Z"]
    assert high.lag_hours == 2.0
    assert low.lag_hours == 8.0
    assert high.lag_hours < low.lag_hours
    assert routed.arrivals[0].arrival_time == "2020-07-01T03:00:00Z"
    assert routed.t_join == "2020-07-01T03:00:00Z"
    assert routed.peak_discharge_m3s == 800.0
    assert routed.peak_arrival_time == "2020-07-01T03:00:00Z"
    assert routed.peak_discharge_m3s == 4.0 * 200.0
    assert squared.area_ratio == 16.0
    assert squared.peak_discharge_m3s == 3200.0


def test_variable_lag_clamps_to_the_configured_bounds() -> None:
    routed = route_variable_lag(
        [
            ("2020-07-01T00:00:00Z", 1.0),
            ("2020-07-01T01:00:00Z", 10_000.0),
        ],
        VariableLag(
            tau_ref_hours=10.0,
            q_ref_m3s=100.0,
            beta=1.0,
            tau_min_hours=1.0,
            tau_max_hours=20.0,
        ),
        as_of="2020-07-01T01:00:00Z",
    )
    by_station = {arrival.station_time: arrival for arrival in routed.arrivals}

    assert by_station["2020-07-01T00:00:00Z"].lag_hours == 20.0
    assert by_station["2020-07-01T01:00:00Z"].lag_hours == 1.0
    assert by_station["2020-07-01T01:00:00Z"].arrival_time == "2020-07-01T02:00:00Z"


def test_station_and_forecast_splices_observed_lag_and_remainder_rain() -> None:
    loss = ScsCnParameters(cn_ii=100)
    watershed = Watershed(area_km2=30, lag_hours=1.0)
    rainfall = _forecast_rain()
    blended = blend_station_and_forecast(
        _station_flow(),
        rainfall,
        loss,
        watershed,
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=_observed_still_dry(),
    )
    heavier_tail = blend_station_and_forecast(
        _station_flow(),
        _forecast_rain(tail=40.0),
        loss,
        watershed,
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=_observed_still_dry(),
    )
    heavier_before_join = blend_station_and_forecast(
        _station_flow(),
        _forecast_rain(early=50.0),
        loss,
        watershed,
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=_observed_still_dry(),
    )
    oracle = route_rainfall(
        rainfall,
        loss,
        watershed,
        timestep_hours=1.0,
        t_join=_T_JOIN,
        observed_inflow=((_T_JOIN, 120.0),),
    )

    assert blended.as_of == _AS_OF
    assert blended.t_join == _T_JOIN
    assert blended.join_discharge_m3s == 120.0
    assert blended.routed_peak_discharge_m3s == 120.0
    assert blended.routed_peak_time == _T_JOIN
    assert blended.hydrograph.points[0].time == "2020-06-30T15:00:00Z"
    assert blended.hydrograph.points[-1].time is not None
    assert blended.hydrograph.points[-1].time >= "2020-07-01T15:00:00Z"
    assert blended.hydrograph.discharge_at("2020-06-30T15:00:00Z") == 0.0
    assert blended.hydrograph.discharge_at("2020-06-30T21:00:00Z") == 0.0
    for time, discharge in _observed_still_dry():
        assert blended.hydrograph.discharge_at(time) == discharge
    assert blended.hydrograph.discharge_at("2020-07-01T01:00:00Z") == 0.0
    assert blended.hydrograph.discharge_at("2020-07-01T04:00:00Z") == 30.0
    assert blended.hydrograph.discharge_at("2020-07-01T05:00:00Z") == 80.0
    assert blended.hydrograph.discharge_at(_T_JOIN) == 120.0
    for point in oracle.points:
        assert point.time is not None
        if point.time > _T_JOIN:
            assert blended.hydrograph.discharge_at(point.time) == pytest.approx(
                point.discharge_m3s
            )
    for point in blended.hydrograph.points:
        assert point.time is not None
        if point.time <= _T_JOIN:
            assert heavier_tail.hydrograph.discharge_at(point.time) == pytest.approx(
                point.discharge_m3s
            )
    assert heavier_tail.hydrograph.discharge_at("2020-07-01T10:00:00Z") != pytest.approx(
        blended.hydrograph.discharge_at("2020-07-01T10:00:00Z")
    )
    assert abs(
        heavier_tail.hydrograph.discharge_at("2020-07-01T10:00:00Z")
        - blended.hydrograph.discharge_at("2020-07-01T10:00:00Z")
    ) > 1.0
    assert _discharges(heavier_before_join.hydrograph) == pytest.approx(
        _discharges(blended.hydrograph)
    )


def test_remainder_after_join_keeps_initial_abstraction_from_the_full_storm() -> None:
    loss = ScsCnParameters(cn_ii=75)
    watershed = Watershed(area_km2=20, lag_hours=1.0)
    observed = _observed_still_dry()
    wet_before_join = blend_station_and_forecast(
        _station_flow(),
        _rain({2: 40.0, 9: 15.0}),
        loss,
        watershed,
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=observed,
    )
    dry_before_join = blend_station_and_forecast(
        _station_flow(),
        _rain({2: 0.0, 9: 15.0}),
        loss,
        watershed,
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=observed,
    )
    probe = "2020-07-01T10:00:00Z"

    assert wet_before_join.hydrograph.discharge_at(_AS_OF) == 0.0
    wet_at_join = wet_before_join.hydrograph.discharge_at(_T_JOIN)
    dry_at_join = dry_before_join.hydrograph.discharge_at(_T_JOIN)
    assert wet_at_join == dry_at_join == 120.0
    assert wet_before_join.hydrograph.discharge_at("2020-07-01T07:00:00Z") == 120.0
    assert wet_before_join.hydrograph.discharge_at("2020-07-01T08:00:00Z") == 120.0
    assert wet_before_join.hydrograph.discharge_at("2020-07-01T15:00:00Z") == 0.0
    assert wet_before_join.hydrograph.discharge_at(probe) > dry_before_join.hydrograph.discharge_at(
        probe
    )


def test_fractional_lag_interpolates_between_as_of_and_t_join() -> None:
    # τ(200) = 3 * (200/100)^(-1) = 1.5 h, so t_join is 04:30, not an hourly step.
    as_of = "2020-07-01T03:00:00Z"
    t_join = "2020-07-01T04:30:00Z"
    station = [
        ("2020-07-01T00:00:00Z", 100.0),
        ("2020-07-01T01:00:00Z", 50.0),
        ("2020-07-01T03:00:00Z", 200.0),
    ]
    observed = _hourly((0.0, 0.0, 0.0, 0.0))
    lag = VariableLag(tau_ref_hours=3.0, q_ref_m3s=100.0, beta=1.0)
    loss = ScsCnParameters(cn_ii=100)
    watershed = Watershed(area_km2=30, lag_hours=1.0)
    rainfall = _rain({2: 4.0, 3: 4.0, 4: 8.0, 5: 12.0, 6: 6.0, 7: 6.0})
    heavier_after_join = _rain({2: 4.0, 3: 4.0, 4: 8.0, 5: 12.0, 6: 40.0, 7: 40.0})
    blended = blend_station_and_forecast(
        station,
        rainfall,
        loss,
        watershed,
        lag,
        as_of=as_of,
        timestep_hours=1.0,
        observed_inflow=observed,
    )
    changed = blend_station_and_forecast(
        station,
        heavier_after_join,
        loss,
        watershed,
        lag,
        as_of=as_of,
        timestep_hours=1.0,
        observed_inflow=observed,
    )

    assert blended.t_join == t_join
    assert blended.join_discharge_m3s == 200.0
    assert blended.routed_peak_time == t_join
    assert all(point.time != t_join for point in blended.hydrograph.points)
    for time, discharge in observed:
        assert blended.hydrograph.discharge_at(time) == discharge
    assert blended.hydrograph.discharge_at(as_of) == 0.0
    # 04:00 is the only hourly instant in (03:00, 04:30]; anchors are 100 at 03:00 and 200 at 04:30.
    assert blended.hydrograph.discharge_at("2020-07-01T04:00:00Z") == pytest.approx(500.0 / 3.0)
    assert blended.hydrograph.discharge_at("2020-07-01T05:00:00Z") > 200.0
    assert changed.hydrograph.discharge_at("2020-07-01T04:00:00Z") == pytest.approx(
        blended.hydrograph.discharge_at("2020-07-01T04:00:00Z")
    )
    assert changed.hydrograph.discharge_at(as_of) == 0.0
    assert changed.hydrograph.discharge_at("2020-07-01T06:00:00Z") != pytest.approx(
        blended.hydrograph.discharge_at("2020-07-01T06:00:00Z")
    )


def test_remainder_tail_recedes_to_zero_after_the_unit_hydrograph() -> None:
    blended = blend_station_and_forecast(
        _station_flow(),
        _rain({7: 20.0}),
        ScsCnParameters(cn_ii=100),
        Watershed(area_km2=25, lag_hours=1.0),
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=_observed_still_dry(),
    )

    assert blended.hydrograph.discharge_at(_T_JOIN) == blended.join_discharge_m3s == 120.0
    assert blended.hydrograph.discharge_at("2020-07-01T08:00:00Z") > 120.0
    assert blended.hydrograph.discharge_at("2020-07-01T15:00:00Z") == 0.0


def test_no_positive_remainder_excess_leaves_later_window_slots_at_zero() -> None:
    blended = blend_station_and_forecast(
        _station_flow(),
        _rain({4: 4.0, 9: 5.0}),
        ScsCnParameters(cn_ii=75),
        Watershed(area_km2=20, lag_hours=1.0),
        _three_hour_lag(),
        as_of=_AS_OF,
        timestep_hours=1.0,
        observed_inflow=_observed_still_dry(),
    )

    assert blended.hydrograph.discharge_at(_AS_OF) == 0.0
    assert blended.hydrograph.discharge_at("2020-07-01T04:00:00Z") == 30.0
    assert blended.hydrograph.discharge_at(_T_JOIN) == 120.0
    assert blended.hydrograph.discharge_at("2020-07-01T07:00:00Z") == 0.0
    assert blended.hydrograph.discharge_at("2020-07-01T09:00:00Z") == 0.0
    assert blended.hydrograph.discharge_at("2020-07-01T12:00:00Z") == 0.0
    assert blended.hydrograph.discharge_at("2020-07-01T15:00:00Z") == 0.0


def test_station_and_forecast_rejects_a_hyetograph_that_ends_at_as_of() -> None:
    with pytest.raises(ChannelRoutingError) as error:
        blend_station_and_forecast(
            _station_flow(),
            _rain({1: 10.0, 2: 10.0, 3: 10.0}),
            ScsCnParameters(cn_ii=80),
            Watershed(area_km2=20, lag_hours=1.0),
            _three_hour_lag(),
            as_of=_AS_OF,
            timestep_hours=1.0,
            observed_inflow=_observed_still_dry(),
        )

    assert error.value.code is ChannelRoutingErrorCode.MISSING_FUTURE_RAINFALL
    assert "as_of" in error.value.message
    assert "rainfall" in error.value.message


def test_muskingum_uses_the_nearest_stable_subreach_count() -> None:
    routed = route_muskingum(
        _hourly((0.0, 10.0, 0.0, 0.0, 0.0, 0.0)),
        MuskingumReach(k_hours=5.1, x=0.49),
        timestep_hours=1.0,
    )

    assert routed.coefficients.n_subreaches == 5
    assert routed.coefficients.c0 >= 0.0
    assert routed.coefficients.c2 >= 0.0
    coefficient_sum = routed.coefficients.c0 + routed.coefficients.c1 + routed.coefficients.c2
    assert coefficient_sum == pytest.approx(1.0)


def test_muskingum_near_boundary_does_not_outflow_below_zero() -> None:
    routed = route_muskingum(
        _hourly((0.0, 10.0, 4.0, 0.0)),
        MuskingumReach(k_hours=1.000000001, x=0.5),
        timestep_hours=1.0,
    )

    assert routed.coefficients.c0 >= 0.0
    assert routed.coefficients.c2 >= 0.0
    coefficient_sum = routed.coefficients.c0 + routed.coefficients.c1 + routed.coefficients.c2
    assert coefficient_sum == pytest.approx(1.0)
    assert all(discharge >= 0.0 for discharge in _discharges(routed.hydrograph))


def test_muskingum_rejects_illegal_weighting_storage_and_an_unstable_timestep() -> None:
    with pytest.raises(ChannelRoutingError) as weighting:
        MuskingumReach(k_hours=1.0, x=0.8)
    assert weighting.value.code is ChannelRoutingErrorCode.INVALID_X

    with pytest.raises(ChannelRoutingError) as storage:
        MuskingumReach(k_hours=0.0, x=0.2)
    assert storage.value.code is ChannelRoutingErrorCode.INVALID_K

    start = datetime(2020, 7, 1, tzinfo=UTC)
    inflow = [
        (_stamp(start), 0.0),
        (_stamp(start + timedelta(hours=3)), 10.0),
        (_stamp(start + timedelta(hours=6)), 0.0),
    ]
    with pytest.raises(ChannelRoutingError) as timestep:
        route_muskingum(inflow, MuskingumReach(k_hours=1.0, x=0.0), timestep_hours=3.0)
    assert timestep.value.code is ChannelRoutingErrorCode.UNSTABLE_TIMESTEP

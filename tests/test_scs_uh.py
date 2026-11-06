"""Public behavior of SCS triangular unit-hydrograph routing and the degraded triangle."""

from __future__ import annotations

from datetime import UTC, datetime, timedelta

import pytest

from hydro_agent.core import (
    RegionalQuantity,
    RegionalTiming,
    ScsCnParameters,
    ScsUhError,
    ScsUhErrorCode,
    TimedRainfall,
    TimingMethod,
    Watershed,
    compute_excess,
    excess_after,
    generalized_triangle_inflow,
    route_excess,
    route_rainfall,
    scs_unit_hydrograph,
)

_AS_OF = "2020-01-01T06:00:00Z"
# 3000 ft, 4% slope, CN 75, evaluated with the US customary NRCS lag equation.
_NRCS_LAG_HOURS = 0.4443128518305718
_KIRPICH_TC_HOURS = 2.2988653185002845
_REGIONAL_RISING_HOURS = 15.811388300841896


def _utc(stamp: str) -> datetime:
    text = f"{stamp[:-1]}+00:00" if stamp.endswith(("Z", "z")) else stamp
    return datetime.fromisoformat(text).astimezone(UTC)


def _stamp(instant: datetime) -> str:
    return instant.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hourly_rain(start: datetime, depths: tuple[float, ...]) -> list[TimedRainfall]:
    return [
        TimedRainfall(_stamp(start + timedelta(hours=index + 1)), depth)
        for index, depth in enumerate(depths)
    ]


@pytest.mark.parametrize("area_km2", [10, 10.5, 250, 1234.5, 5000])
@pytest.mark.parametrize("timestep_hours", [0.25, 0.5, 1.0])
@pytest.mark.parametrize("lag_hours", [0.0, 0.5, 4.0])
def test_one_millimetre_of_excess_conserves_volume(
    area_km2: float,
    timestep_hours: float,
    lag_hours: float,
) -> None:
    hydrograph = route_excess(
        (1.0,),
        Watershed(area_km2=area_km2, lag_hours=lag_hours),
        timestep_hours=timestep_hours,
    )

    assert hydrograph.volume_m3 == pytest.approx(1000.0 * area_km2, rel=1e-6)
    assert hydrograph.peak_discharge_m3s > 0.0


def test_routing_a_rainfall_hyetograph_uses_scs_cn_excess_not_the_storm_total() -> None:
    rainfall = (5.0, 10.0, 20.0, 25.0, 15.0, 5.0)
    loss = ScsCnParameters(cn_ii=75)
    excess = compute_excess(rainfall, loss)
    hydrograph = route_rainfall(
        rainfall,
        loss,
        Watershed(area_km2=42, lag_hours=2),
        timestep_hours=1,
    )

    assert excess.total_excess_mm < excess.total_rainfall_mm
    assert hydrograph.volume_m3 == pytest.approx(1000.0 * 42 * excess.total_excess_mm, rel=1e-6)
    assert hydrograph.volume_m3 != pytest.approx(1000.0 * 42 * excess.total_rainfall_mm, rel=1e-3)


def test_delaying_the_rainfall_peak_delays_the_inflow_peak_by_the_same_amount() -> None:
    depths = (1.0, 4.0, 20.0, 4.0, 1.0)
    start = datetime(2020, 6, 1, tzinfo=UTC)
    watershed = Watershed(area_km2=50, lag_hours=1.5)
    loss = ScsCnParameters(cn_ii=100)
    early = route_rainfall(
        _hourly_rain(start, depths),
        loss,
        watershed,
        timestep_hours=1,
    )
    late = route_rainfall(
        _hourly_rain(start + timedelta(hours=3), depths),
        loss,
        watershed,
        timestep_hours=1,
    )

    assert early.peak_time is not None and late.peak_time is not None
    assert _utc(late.peak_time) - _utc(early.peak_time) == timedelta(hours=3)
    assert [point.discharge_m3s for point in late.points] == pytest.approx(
        [point.discharge_m3s for point in early.points]
    )
    assert all(
        _utc(late_point.time or "") - _utc(early_point.time or "") == timedelta(hours=3)
        for early_point, late_point in zip(early.points, late.points, strict=True)
    )


def test_hours_to_peak_triangle_hits_the_peak_and_recesses_one_to_two() -> None:
    hydrograph = generalized_triangle_inflow(
        peak_discharge_m3s=500,
        as_of=_AS_OF,
        inflow_at_as_of_m3s=50,
        hours_to_peak=4,
        timestep_hours=1,
    )

    assert hydrograph.peak_time == "2020-01-01T10:00:00Z"
    assert hydrograph.peak_hours_from_origin == pytest.approx(4.0)
    assert hydrograph.discharge_at(_AS_OF) == 50
    assert hydrograph.discharge_at("2020-01-01T08:00:00Z") == 275
    assert hydrograph.discharge_at("2020-01-01T10:00:00Z") == 500
    assert hydrograph.discharge_at("2020-01-01T17:00:00Z") == 62.5
    assert hydrograph.discharge_at("2020-01-01T18:00:00Z") == 0
    assert hydrograph.points[0].time == "2019-12-31T18:00:00Z"
    assert hydrograph.points[-1].time == "2020-01-01T18:00:00Z"


def test_hours_to_peak_overrides_tau_hours() -> None:
    hydrograph = generalized_triangle_inflow(
        peak_discharge_m3s=500,
        as_of=_AS_OF,
        inflow_at_as_of_m3s=50,
        hours_to_peak=4,
        tau_hours=9,
        timestep_hours=1,
    )

    assert hydrograph.peak_time == "2020-01-01T10:00:00Z"
    assert hydrograph.discharge_at("2020-01-01T10:00:00Z") == 500
    assert hydrograph.discharge_at("2020-01-01T15:00:00Z") == 187.5


def test_tau_hours_alone_measures_the_rise_from_as_of() -> None:
    hydrograph = generalized_triangle_inflow(
        peak_discharge_m3s=500,
        as_of=_AS_OF,
        inflow_at_as_of_m3s=50,
        tau_hours=9,
        timestep_hours=1,
    )

    assert hydrograph.peak_time == "2020-01-01T15:00:00Z"
    assert hydrograph.discharge_at("2020-01-01T15:00:00Z") == 500


def test_explicit_fall_hours_replaces_the_one_to_two_recession() -> None:
    hydrograph = generalized_triangle_inflow(
        peak_discharge_m3s=500,
        as_of=_AS_OF,
        inflow_at_as_of_m3s=50,
        hours_to_peak=4,
        fall_hours=5,
        timestep_hours=1,
    )

    assert hydrograph.discharge_at("2020-01-01T14:00:00Z") == 100
    assert hydrograph.discharge_at("2020-01-01T15:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T18:00:00Z") == 0


def test_remainder_after_join_is_continuous_and_ignores_earlier_rain_when_cn_is_100() -> None:
    join = "2020-01-01T02:00:00Z"
    start = datetime(2020, 1, 1, tzinfo=UTC)
    rainfall = _hourly_rain(start, (12.0, 12.0, 30.0, 6.0))
    heavier_before_join = _hourly_rain(start, (40.0, 12.0, 30.0, 6.0))
    watershed = Watershed(area_km2=20, lag_hours=2)
    loss = ScsCnParameters(cn_ii=100)
    observed = ((join, 80.0),)

    routed = route_rainfall(
        rainfall,
        loss,
        watershed,
        timestep_hours=1,
        t_join=join,
        observed_inflow=observed,
    )
    changed_history = route_rainfall(
        heavier_before_join,
        loss,
        watershed,
        timestep_hours=1,
        t_join=join,
        observed_inflow=observed,
    )
    tail = route_rainfall(rainfall[2:], loss, watershed, timestep_hours=1)
    full = route_rainfall(rainfall, loss, watershed, timestep_hours=1)

    assert routed.discharge_at(join) == 80
    for point in tail.points:
        assert point.time is not None
        assert routed.discharge_at(point.time) == pytest.approx(80.0 + point.discharge_m3s)
        assert changed_history.discharge_at(point.time) == pytest.approx(
            routed.discharge_at(point.time)
        )
    remainder_direct = routed.discharge_at("2020-01-01T03:00:00Z") - 80.0
    assert full.discharge_at("2020-01-01T03:00:00Z") != pytest.approx(remainder_direct)


def test_remainder_keeps_initial_abstraction_already_taken_before_the_join() -> None:
    join = "2020-01-01T02:00:00Z"
    start = datetime(2020, 1, 1, tzinfo=UTC)
    rainfall = _hourly_rain(start, (20.0, 20.0, 10.0))
    loss = ScsCnParameters(cn_ii=75)
    watershed = Watershed(area_km2=15, lag_hours=1.5)
    remainder = excess_after(rainfall, loss, join)

    routed = route_rainfall(
        rainfall,
        loss,
        watershed,
        timestep_hours=1,
        t_join=join,
        observed_inflow=((join, 0.0),),
    )
    restarted = route_rainfall(rainfall[2:], loss, watershed, timestep_hours=1)

    assert remainder.remainder_excess_mm > 0.0
    assert routed.volume_m3 == pytest.approx(
        1000.0 * 15 * remainder.remainder_excess_mm,
        rel=1e-6,
    )
    assert restarted.volume_m3 == 0.0
    assert routed.discharge_at(join) == 0.0


def test_as_of_window_fills_unmeasured_slots_with_zero_not_baseflow() -> None:
    hydrograph = generalized_triangle_inflow(
        peak_discharge_m3s=400,
        as_of=_AS_OF,
        inflow_at_as_of_m3s=50,
        hours_to_peak=2,
        timestep_hours=1,
        observed_inflow=(("2020-01-01T05:00:00Z", 50.0),),
    )

    assert hydrograph.discharge_at("2019-12-31T18:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T04:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T05:00:00Z") == 50
    assert hydrograph.discharge_at(_AS_OF) == 50
    assert hydrograph.discharge_at("2020-01-01T07:00:00Z") == 225
    assert hydrograph.discharge_at("2020-01-01T12:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T13:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T18:00:00Z") == 0
    assert hydrograph.points[0].time == "2019-12-31T18:00:00Z"
    assert hydrograph.points[-1].time == "2020-01-01T18:00:00Z"


def test_as_of_window_does_not_spread_convolution_baseflow_into_empty_slots() -> None:
    as_of = "2020-01-01T12:00:00Z"
    hydrograph = route_excess(
        (TimedRainfall("2020-01-01T08:00:00Z", 5.0),),
        Watershed(area_km2=50, lag_hours=1.5),
        timestep_hours=1,
        baseflow_m3s=30,
        as_of=as_of,
    )

    assert hydrograph.discharge_at("2020-01-01T00:00:00Z") == 0
    assert hydrograph.discharge_at("2020-01-01T08:00:00Z") > 30
    assert hydrograph.discharge_at("2020-01-02T00:00:00Z") == 0


def test_triangular_unit_hydrograph_uses_scs_time_to_peak_and_recession() -> None:
    hydrograph = scs_unit_hydrograph(
        Watershed(area_km2=100, lag_hours=2),
        timestep_hours=1,
    )

    assert hydrograph.time_to_peak_hours == pytest.approx(2.5)
    assert hydrograph.recession_ratio == pytest.approx(1.67)
    assert hydrograph.recession_hours == pytest.approx(1.67 * 2.5)
    assert hydrograph.base_hours == pytest.approx(2.67 * 2.5)
    assert hydrograph.peak_rate_coefficient == pytest.approx(2.0 / (3.6 * 2.67))
    assert hydrograph.peak_discharge_m3s == pytest.approx(
        2.0 * 1000.0 * 100.0 / (3600.0 * hydrograph.base_hours)
    )
    assert hydrograph.ordinates[0].discharge_m3s == 0.0
    assert hydrograph.ordinates[0].dimensionless == 0.0


def test_time_of_concentration_sets_lag_to_six_tenths() -> None:
    hydrograph = scs_unit_hydrograph(
        Watershed(area_km2=10, time_of_concentration_hours=5),
        timestep_hours=1,
    )

    assert hydrograph.lag_hours == pytest.approx(3.0)
    assert hydrograph.time_of_concentration_hours == pytest.approx(5.0)
    assert hydrograph.time_to_peak_hours == pytest.approx(3.5)


def test_nrcs_lag_matches_the_customary_equation_for_a_3000_foot_watershed() -> None:
    hydrograph = scs_unit_hydrograph(
        Watershed(
            area_km2=10,
            channel_length_km=0.9143999707392009,
            slope=0.04,
            curve_number=75,
            timing_method=TimingMethod.NRCS_LAG,
        ),
        timestep_hours=1,
    )

    assert hydrograph.lag_hours == pytest.approx(_NRCS_LAG_HOURS, rel=1e-12)
    assert hydrograph.time_of_concentration_hours == pytest.approx(_NRCS_LAG_HOURS / 0.6)


def test_kirpich_time_of_concentration_sets_lag() -> None:
    hydrograph = scs_unit_hydrograph(
        Watershed(
            area_km2=10,
            channel_length_km=10,
            slope=0.01,
            timing_method=TimingMethod.KIRPICH,
        ),
        timestep_hours=1,
    )

    assert hydrograph.time_of_concentration_hours == pytest.approx(_KIRPICH_TC_HOURS, rel=1e-12)
    assert hydrograph.lag_hours == pytest.approx(0.6 * _KIRPICH_TC_HOURS, rel=1e-12)


def test_regional_rising_time_is_the_unit_hydrograph_time_to_peak() -> None:
    hydrograph = scs_unit_hydrograph(
        Watershed(
            area_km2=100,
            channel_length_km=10,
            slope=0.01,
            timing_method=TimingMethod.REGIONAL,
            regional=RegionalTiming(
                coefficient=0.5,
                area_exponent=0.2,
                length_exponent=0.5,
                slope_exponent=-0.3,
                quantity=RegionalQuantity.RISING_TIME,
            ),
        ),
        timestep_hours=1,
    )

    assert hydrograph.time_to_peak_hours == pytest.approx(_REGIONAL_RISING_HOURS, rel=1e-12)


def test_missing_timing_and_nonpositive_area_are_rejected() -> None:
    with pytest.raises(ScsUhError) as area_error:
        scs_unit_hydrograph(Watershed(area_km2=0, lag_hours=1), timestep_hours=1)
    assert area_error.value.code is ScsUhErrorCode.INVALID_AREA

    with pytest.raises(ScsUhError) as timing_error:
        scs_unit_hydrograph(Watershed(area_km2=10), timestep_hours=1)
    assert timing_error.value.code is ScsUhErrorCode.INVALID_TIMING

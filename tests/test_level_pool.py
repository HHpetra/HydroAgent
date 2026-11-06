"""Level-pool routing: Modified Puls, gated control, and reverse inflow."""

from __future__ import annotations

import math
from datetime import UTC, datetime, timedelta

import pytest

from hydro_agent.core import (
    DischargeCurve,
    GateGeometry,
    LevelPoolError,
    LevelPoolErrorCode,
    LevelPoolMode,
    LevelPoolResult,
    StorageCurve,
    gate_discharge,
    gate_opening,
    reverse_inflow,
    route_constant_target,
    route_gate_fault,
    route_level_pool,
    route_pre_release,
    route_schedule,
    route_to_q_safe,
    route_to_z_max,
    route_ungated,
    route_ungated_spill,
    stage_for_storage,
    storage_volume,
    ungated_discharge,
)

_AS_OF = "2020-06-01T00:00:00Z"
_DAY = datetime(2020, 6, 1, tzinfo=UTC)
_SILL = 82.7
_FLOOD_LIMIT = 84.2
_Q_SAFE = 80.0
_G = 9.81


def _stamp(instant: datetime) -> str:
    return instant.astimezone(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def _hourly(values: tuple[float, ...] | list[float], *, start: datetime) -> list[tuple[str, float]]:
    return [
        (_stamp(start + timedelta(hours=index)), value) for index, value in enumerate(values)
    ]


def _wide_storage() -> StorageCurve:
    stages = tuple(float(stage) for stage in range(70, 121, 2))
    return StorageCurve(
        stages_m=stages,
        volumes_m3=tuple(500_000.0 * (stage - 70.0) for stage in stages),
    )


def _wide_gates() -> GateGeometry:
    return GateGeometry(
        n_gates=2,
        width_m=2000.0,
        sill_elevation_m=_SILL,
        max_opening_m=5.0,
        weir_coefficient=0.5,
        orifice_coefficient=0.65,
    )


def _flood(start_hour: int) -> list[tuple[str, float]]:
    """Steep wave: zero through ``start_hour``, then 200, 400, 200, and zeros."""

    values = [0.0] * 49
    values[start_hour + 1] = 200.0
    values[start_hour + 2] = 400.0
    values[start_hour + 3] = 200.0
    return _hourly(values, start=_DAY)


def _balance(result: LevelPoolResult) -> tuple[float, float]:
    future = [point for point in result.points if point.hours_from_as_of >= 0.0]
    dt_s = result.timestep_hours * 3600.0
    worst = 0.0
    net = 0.0
    for left, right in zip(future, future[1:], strict=False):
        residual = (
            0.5 * (left.inflow_m3s + right.inflow_m3s) * dt_s
            - 0.5 * (left.outflow_m3s + right.outflow_m3s) * dt_s
            - (right.storage_m3 - left.storage_m3)
        )
        worst = max(worst, abs(residual))
        net += residual
    return worst, net


def test_storage_curve_is_strictly_monotone_and_bidirectional() -> None:
    curve = StorageCurve(stages_m=(80.0, 82.0, 86.0), volumes_m3=(0.0, 1_000_000.0, 3_000_000.0))
    assert storage_volume(curve, 84.0) == pytest.approx(2_000_000.0)
    assert stage_for_storage(curve, 2_000_000.0) == pytest.approx(84.0)
    assert stage_for_storage(curve, storage_volume(curve, 83.2)) == pytest.approx(83.2)
    with pytest.raises(LevelPoolError) as error:
        StorageCurve(stages_m=(80.0, 82.0), volumes_m3=(10.0, 10.0))
    assert error.value.code is LevelPoolErrorCode.NON_MONOTONE_STORAGE


def test_ungated_discharge_is_zero_at_and_below_the_crest() -> None:
    curve = DischargeCurve(stages_m=(9.0, 10.0, 12.0), discharges_m3s=(0.0, 0.0, 40.0))
    assert curve.crest_stage_m == pytest.approx(10.0)
    assert ungated_discharge(curve, 8.0) == 0.0
    assert ungated_discharge(curve, 10.0) == 0.0
    assert ungated_discharge(curve, 11.0) == pytest.approx(20.0)


def test_ungated_modified_puls_conserves_mass_and_peaks_on_the_recession() -> None:
    stages = tuple(float(stage) for stage in range(0, 31))
    storage = StorageCurve(
        stages_m=stages,
        volumes_m3=tuple(1_000_000.0 * stage for stage in stages),
    )
    discharge = DischargeCurve(
        stages_m=(0.0, 10.0, 12.0, 16.0, 20.0, 30.0),
        discharges_m3s=(0.0, 0.0, 40.0, 160.0, 320.0, 800.0),
    )
    inflows = [0.0] * 36
    for hour, value in enumerate((20.0, 50.0, 90.0, 140.0, 180.0, 140.0, 90.0, 40.0, 10.0)):
        inflows[4 + hour] = value
    result = route_ungated(
        storage,
        discharge,
        _hourly(inflows, start=_DAY),
        z0_m=10.0,
        as_of=_AS_OF,
        timestep_hours=1.0,
    )
    worst, net = _balance(result)
    assert worst < 1e-6
    assert abs(net) < 1e-4
    assert result.points[0].hours_from_as_of == pytest.approx(-12.0)
    assert result.points[-1].hours_from_as_of >= 12.0
    future = [point for point in result.points if point.hours_from_as_of >= 0.0]
    peak_stage = max(future, key=lambda point: point.stage_m)
    peak_out = max(future, key=lambda point: point.outflow_m3s)
    assert peak_stage.time == peak_out.time
    assert result.t_z_max == result.t_q_out_max
    peak_index = future.index(peak_out)
    assert peak_index > 0
    before = future[peak_index - 1]
    after = future[peak_index]
    assert before.inflow_m3s >= before.outflow_m3s
    assert after.outflow_m3s >= after.inflow_m3s
    assert any(point.inflow_m3s > point.outflow_m3s for point in future[:peak_index])


def test_weir_orifice_scaling_and_opening_inverse() -> None:
    gates = GateGeometry(
        n_gates=2,
        width_m=5.0,
        sill_elevation_m=100.0,
        max_opening_m=4.0,
        weir_coefficient=0.48,
        orifice_coefficient=0.62,
    )
    orifice = gate_discharge(gates, 108.0, 0.4)
    assert orifice == pytest.approx(31.0703521705, abs=1e-6)
    assert gate_discharge(gates, 108.0, 0.8) == pytest.approx(2.0 * orifice, rel=1e-9)
    assert gate_discharge(gates, 102.0, 0.4) == pytest.approx(orifice / math.sqrt(4.0), rel=1e-9)

    weir_head_2 = gate_discharge(gates, 102.0, gates.max_opening_m)
    weir_head_4 = gate_discharge(gates, 104.0, gates.max_opening_m)
    assert weir_head_2 == pytest.approx(60.1361654913, abs=1e-6)
    assert weir_head_4 / weir_head_2 == pytest.approx(2.0**1.5, rel=1e-9)

    # Fully open uses the weir even when e/H0 is still in the orifice range.
    low_gate = GateGeometry(
        n_gates=2,
        width_m=5.0,
        sill_elevation_m=100.0,
        max_opening_m=0.3,
        weir_coefficient=0.48,
        orifice_coefficient=0.62,
    )
    forced = gate_discharge(low_gate, 110.0, 0.3)
    orifice_if_not_weir = 0.62 * 2 * 5.0 * 0.3 * math.sqrt(2.0 * _G * 10.0)
    assert forced == pytest.approx(weir_head_2 * (5.0**1.5), rel=1e-9)
    assert forced > orifice_if_not_weir * 10.0

    command = gate_opening(gates, 108.0, orifice)
    assert command.limited_by_capacity is False
    assert command.opening_m == pytest.approx(0.4, abs=1e-6)
    assert gate_discharge(gates, 108.0, command.opening_m) == pytest.approx(orifice, abs=1e-6)
    capacity = gate_discharge(gates, 102.0, gates.max_opening_m)
    limited = gate_opening(gates, 102.0, capacity + 25.0)
    assert limited.opening_m == pytest.approx(gates.max_opening_m)
    assert limited.discharge_m3s == pytest.approx(capacity, abs=1e-6)
    assert limited.limited_by_capacity is True


def test_pre_release_holds_the_sill_then_the_flood_limit() -> None:
    storage = _wide_storage()
    gates = _wide_gates()
    inflow = _flood(8)
    pre = route_pre_release(
        storage,
        gates,
        inflow,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
    )
    closed = route_to_q_safe(
        storage,
        gates,
        inflow,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
        wait_closed_until_z=_FLOOD_LIMIT,
    )
    for result in (pre, closed):
        worst, net = _balance(result)
        assert worst < 1e-6
        assert abs(net) < 1e-4
        assert result.points[0].hours_from_as_of == pytest.approx(-12.0)
        assert max(point.hours_from_as_of for point in result.points) >= 12.0

    history = [point for point in pre.points if point.hours_from_as_of <= 0.0]
    assert all(point.stage_m == pytest.approx(_FLOOD_LIMIT) for point in history)
    assert all(point.outflow_m3s == pytest.approx(0.0) for point in history)

    draw = [
        point
        for point in pre.points
        if 0.0 < point.hours_from_as_of <= 2.0
    ]
    assert draw
    assert all(point.outflow_m3s == pytest.approx(_Q_SAFE) for point in draw)
    draw_stages = [point.stage_m for point in draw]
    assert draw_stages == sorted(draw_stages, reverse=True)

    assert pre.time_reached_target is not None
    held = [
        point
        for point in pre.points
        if pre.time_reached_target is not None
        and point.time >= pre.time_reached_target
        and point.hours_from_as_of <= 8.0
    ]
    assert held
    assert all(point.stage_m == pytest.approx(_SILL, abs=1e-4) for point in held)
    assert all(point.outflow_m3s == pytest.approx(0.0, abs=1e-4) for point in held)
    assert pre.z_end_m == pytest.approx(_FLOOD_LIMIT, abs=1e-4)
    assert pre.z_max_m < closed.z_max_m - 1.0

    before_flood = [point for point in closed.points if 0.0 <= point.hours_from_as_of <= 8.0]
    assert all(point.stage_m == pytest.approx(_FLOOD_LIMIT, abs=1e-6) for point in before_flood)
    assert all(point.outflow_m3s == pytest.approx(0.0, abs=1e-6) for point in before_flood)
    future_closed = [point for point in closed.points if point.hours_from_as_of >= 0.0]
    assert all(point.stage_m >= _FLOOD_LIMIT - 1e-6 for point in future_closed)
    release = [
        point
        for point in closed.points
        if point.stage_m > _FLOOD_LIMIT + 1e-3
    ]
    assert release
    assert all(point.outflow_m3s == pytest.approx(_Q_SAFE, abs=1e-4) for point in release)
    assert closed.z_end_m == pytest.approx(_FLOOD_LIMIT, abs=1e-4)
    _assert_gate_does_not_reopen(closed)


def test_pre_release_holds_flood_limit_when_the_sill_was_never_reached() -> None:
    result = route_pre_release(
        _wide_storage(),
        _wide_gates(),
        _flood(1),
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
    )
    assert result.time_reached_target is None
    assert result.z_end_m == pytest.approx(_FLOOD_LIMIT, abs=1e-4)
    future = [point for point in result.points if point.hours_from_as_of >= 0.0]
    peak_index = max(range(len(future)), key=lambda index: future[index].stage_m)
    assert all(point.stage_m >= _FLOOD_LIMIT - 1e-4 for point in future[peak_index:])
    worst, _net = _balance(result)
    assert worst < 1e-6


def test_pre_release_target_below_the_sill_is_clamped() -> None:
    storage = _wide_storage()
    gates = _wide_gates()
    inflow = _flood(8)
    lowered = route_pre_release(
        storage,
        gates,
        inflow,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
        z_target_m=80.0,
    )
    default = route_pre_release(
        storage,
        gates,
        inflow,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
    )
    assert [point.stage_m for point in lowered.points] == pytest.approx(
        [point.stage_m for point in default.points]
    )


def test_reverse_inflow_recovers_a_hidden_segment() -> None:
    stages = tuple(float(stage) for stage in range(0, 31))
    storage = StorageCurve(
        stages_m=stages,
        volumes_m3=tuple(1_000_000.0 * stage for stage in stages),
    )
    discharge = DischargeCurve(
        stages_m=(0.0, 10.0, 14.0, 18.0, 24.0, 30.0),
        discharges_m3s=(0.0, 0.0, 80.0, 200.0, 420.0, 700.0),
    )
    values = [0.0] * 30
    for hour, value in enumerate((30.0, 70.0, 120.0, 160.0, 130.0, 80.0, 40.0, 15.0)):
        values[2 + hour] = value
    inflow = _hourly(values, start=_DAY)
    routed = route_ungated(
        storage,
        discharge,
        inflow,
        z0_m=10.0,
        as_of=_AS_OF,
        timestep_hours=1.0,
    )
    t_fail = _stamp(_DAY + timedelta(hours=4))
    hidden_as_of = _stamp(_DAY + timedelta(hours=12))
    observed = [
        (point.time, point.stage_m, point.outflow_m3s)
        for point in routed.points
        if t_fail <= point.time <= hidden_as_of
    ]
    recovered = reverse_inflow(
        storage,
        [(time, stage) for time, stage, _outflow in observed],
        [(time, outflow) for time, _stage, outflow in observed],
        t_fail=t_fail,
        as_of=hidden_as_of,
        timestep_hours=1.0,
        observed_inflow=[(time, discharge_m3s) for time, discharge_m3s in inflow if time <= t_fail],
        forecast_inflow=[
            (time, discharge_m3s) for time, discharge_m3s in inflow if time > hidden_as_of
        ],
    )
    truth = {time: discharge_m3s for time, discharge_m3s in inflow}
    gap = [point for point in recovered.points if t_fail < point.time <= hidden_as_of]
    assert gap
    for point in gap:
        assert point.discharge_m3s == pytest.approx(truth[point.time], abs=1e-6)
    assert recovered.points[0].hours_from_as_of <= -12.0
    assert recovered.points[-1].hours_from_as_of >= 12.0


def test_cross_midnight_observations_stay_on_the_trace() -> None:
    as_of = datetime(2020, 1, 1, 6, tzinfo=UTC)
    start = datetime(2019, 12, 31, 18, tzinfo=UTC)
    observed_stage: list[tuple[str, float]] = []
    observed_outflow: list[tuple[str, float]] = []
    observed_opening: list[tuple[str, float]] = []
    instant = start
    index = 0
    while instant <= as_of:
        observed_stage.append((_stamp(instant), 83.0 + 0.05 * index))
        observed_outflow.append((_stamp(instant), 4.0 + index))
        observed_opening.append((_stamp(instant), 0.15 + 0.02 * index))
        instant += timedelta(hours=1)
        index += 1
    z0 = observed_stage[-1][1]
    inflow = _hourly((0.0, 0.0, 0.0), start=as_of)
    result = route_to_q_safe(
        _wide_storage(),
        _wide_gates(),
        inflow,
        z0_m=z0,
        as_of=_stamp(as_of),
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
        wait_closed_until_z=z0,
        observed_stage=observed_stage,
        observed_outflow=observed_outflow,
        observed_gate_opening=observed_opening,
    )
    by_time = {point.time: point for point in result.points}
    assert _stamp(start) in by_time
    assert by_time[_stamp(start)].time.startswith("2019-12-31T")
    for (stage_time, stage), (_out_time, outflow), (_open_time, opening) in zip(
        observed_stage,
        observed_outflow,
        observed_opening,
        strict=True,
    ):
        if stage_time >= _stamp(as_of):
            continue
        point = by_time[stage_time]
        assert point.stage_m == pytest.approx(stage)
        assert point.outflow_m3s == pytest.approx(outflow)
        assert point.gate_opening_m == pytest.approx(opening)
    assert by_time[_stamp(as_of)].stage_m == pytest.approx(z0)
    assert min(point.hours_from_as_of for point in result.points) == pytest.approx(-12.0)


def test_stuck_gates_cut_the_peak_discharge_and_raise_the_stage() -> None:
    stages = tuple(float(stage) for stage in range(70, 121, 2))
    storage = StorageCurve(
        stages_m=stages,
        volumes_m3=tuple(400_000.0 * (stage - 70.0) for stage in stages),
    )
    gates = GateGeometry(
        n_gates=4,
        width_m=8.0,
        sill_elevation_m=80.0,
        max_opening_m=4.0,
        weir_coefficient=0.45,
        orifice_coefficient=0.62,
    )
    values = [0.0] * 36
    for hour, value in enumerate((40.0, 120.0, 220.0, 300.0, 220.0, 120.0, 40.0)):
        values[2 + hour] = value
    inflow = _hourly(values, start=_DAY)
    open_spill = route_ungated_spill(
        storage,
        gates,
        inflow,
        z0_m=82.0,
        as_of=_AS_OF,
        timestep_hours=1.0,
    )
    faulted = route_gate_fault(
        storage,
        gates,
        inflow,
        z0_m=82.0,
        as_of=_AS_OF,
        timestep_hours=1.0,
        n_fault=1,
    )
    assert faulted.q_out_max_m3s < open_spill.q_out_max_m3s
    assert faulted.z_max_m > open_spill.z_max_m
    ratio = gate_discharge(gates, 86.0, 4.0, n_fault=1) / gate_discharge(gates, 86.0, 4.0)
    assert ratio == pytest.approx(3.0 / 4.0)


def test_outflow_controls_are_mutually_exclusive() -> None:
    storage = _wide_storage()
    gates = _wide_gates()
    discharge = DischargeCurve(stages_m=(80.0, 90.0), discharges_m3s=(0.0, 100.0))
    inflow = _flood(8)
    with pytest.raises(LevelPoolError) as both:
        route_level_pool(
            storage,
            inflow,
            mode=LevelPoolMode.UNGATED,
            z0_m=_FLOOD_LIMIT,
            as_of=_AS_OF,
            timestep_hours=1.0,
            discharge=discharge,
            gates=gates,
        )
    assert both.value.code is LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT
    with pytest.raises(LevelPoolError) as neither:
        route_level_pool(
            storage,
            inflow,
            mode=LevelPoolMode.PRE_RELEASE,
            z0_m=_FLOOD_LIMIT,
            as_of=_AS_OF,
            timestep_hours=1.0,
            q_safe_m3s=_Q_SAFE,
        )
    assert neither.value.code is LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT
    with pytest.raises(LevelPoolError) as extra:
        route_pre_release(
            storage,
            gates,
            inflow,
            z0_m=_FLOOD_LIMIT,
            as_of=_AS_OF,
            timestep_hours=1.0,
            q_safe_m3s=_Q_SAFE,
            discharge=discharge,
        )
    assert extra.value.code is LevelPoolErrorCode.OUTFLOW_CONTROL_CONFLICT


def test_schedule_constant_target_and_stage_cap() -> None:
    storage = _wide_storage()
    gates = _wide_gates()
    quiet = _hourly((0.0,) * 8, start=_DAY)
    scheduled = route_schedule(
        storage,
        gates,
        quiet,
        [(_AS_OF, 0.2)],
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
    )
    assert any(point.gate_opening_m == pytest.approx(0.2) for point in scheduled.points)
    targeted = route_constant_target(
        storage,
        gates,
        quiet,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_target_m3s=_Q_SAFE,
    )
    held = route_to_q_safe(
        storage,
        gates,
        quiet,
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        q_safe_m3s=_Q_SAFE,
    )
    assert targeted.z_end_m < _FLOOD_LIMIT - 0.5
    assert held.z_end_m == pytest.approx(_FLOOD_LIMIT, abs=1e-4)
    capped = route_to_z_max(
        storage,
        gates,
        _flood(2),
        z0_m=_FLOOD_LIMIT,
        as_of=_AS_OF,
        timestep_hours=1.0,
        z_lim_m=85.5,
        q_safe_m3s=_Q_SAFE,
        z_hold_m=_FLOOD_LIMIT,
    )
    assert capped.z_max_m == pytest.approx(85.5, abs=0.05)
    assert capped.z_end_m == pytest.approx(_FLOOD_LIMIT, abs=1e-3)
    worst, _net = _balance(capped)
    assert worst < 1e-6


def _assert_gate_does_not_reopen(result: LevelPoolResult) -> None:
    openings = [
        0.0 if point.gate_opening_m is None else point.gate_opening_m
        for point in result.points
        if point.hours_from_as_of >= 0.0
    ]
    seen_open = False
    seen_shut = False
    for opening in openings:
        if opening > 1e-4:
            if seen_shut:
                raise AssertionError("gate closed and opened again")
            seen_open = True
        elif seen_open:
            seen_shut = True

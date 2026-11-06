"""Public behavior of the SCS-CN excess rainfall module."""

from __future__ import annotations

import math
import random

import pytest

from hydro_agent.core import (
    AmcCondition,
    ImperviousTreatment,
    ScsCnCase,
    ScsCnError,
    ScsCnErrorCode,
    ScsCnParameters,
    ScsCnResult,
    TimedRainfall,
    compute_excess,
    compute_excess_batch,
    excess_after,
)

# Spec acceptance for the HEC-HMS hand calculation.
_HAND_CALC_TOLERANCE_MM = 1e-6
_BALANCE_TOLERANCE_MM = 1e-9

_SIX_HOUR_MM = (5.0, 10.0, 20.0, 25.0, 15.0, 5.0)
# Cumulative P, cumulative R, incremental excess, cumulative infiltration.
# Independent evaluation of the cumulative SCS equation at CN = 75, λ = 0.2
# (S = 84.666… mm, Ia = 16.933… mm).
_SIX_HOUR_EXPECTED = (
    (5.0, 0.0, 0.0, 0.0),
    (15.0, 0.0, 0.0, 0.0),
    (35.0, 3.177200951763, 3.177200951763, 14.889465714904),
    (60.0, 14.520389700765, 11.343188749003, 28.546276965901),
    (75.0, 23.622637396855, 9.102247696090, 34.444029269812),
    (80.0, 26.922864019254, 3.300226622399, 36.143802647413),
)

_CN_75 = ScsCnParameters(cn_ii=75)


def test_six_hour_storm_matches_hec_hms_cumulative_hand_calc() -> None:
    result = compute_excess(_SIX_HOUR_MM, _CN_75)

    assert round(result.potential_retention_mm, 2) == 84.67
    assert round(result.initial_abstraction_mm, 2) == 16.93
    assert result.curve_number == 75.0
    assert len(result.steps) == len(_SIX_HOUR_EXPECTED)

    paired = zip(_SIX_HOUR_MM, result.steps, _SIX_HOUR_EXPECTED, strict=True)
    for rainfall, step, expected in paired:
        cumulative_p, cumulative_r, incremental_r, infiltration = expected
        assert step.rainfall_mm == rainfall
        assert step.time is None
        assert step.cumulative_rainfall_mm == pytest.approx(
            cumulative_p, abs=_HAND_CALC_TOLERANCE_MM
        )
        assert step.cumulative_excess_mm == pytest.approx(
            cumulative_r, abs=_HAND_CALC_TOLERANCE_MM
        )
        assert step.excess_mm == pytest.approx(incremental_r, abs=_HAND_CALC_TOLERANCE_MM)
        assert step.cumulative_infiltration_mm == pytest.approx(
            infiltration, abs=_HAND_CALC_TOLERANCE_MM
        )
        assert step.loss_mm == pytest.approx(
            rainfall - incremental_r, abs=_HAND_CALC_TOLERANCE_MM
        )

    assert result.total_excess_mm == pytest.approx(26.922864019254, abs=_HAND_CALC_TOLERANCE_MM)
    assert result.total_rainfall_mm == pytest.approx(80.0)
    assert result.runoff_coefficient == pytest.approx(result.total_excess_mm / 80.0)


def test_no_excess_while_cumulative_rain_is_within_initial_abstraction() -> None:
    result = compute_excess((4.0, 4.0, 4.0), _CN_75)

    assert result.total_rainfall_mm == 12.0
    assert result.total_rainfall_mm < result.initial_abstraction_mm
    assert result.total_excess_mm == 0.0
    assert result.initial_abstraction_consumed_mm == 12.0
    assert result.cumulative_infiltration_mm == 0.0
    assert result.runoff_coefficient == 0.0
    assert all(step.excess_mm == 0.0 for step in result.steps)


def test_excess_is_zero_when_cumulative_rain_equals_initial_abstraction() -> None:
    initial_abstraction = 0.2 * (25400.0 / 75.0 - 254.0)
    result = compute_excess((initial_abstraction,), _CN_75)

    assert result.total_excess_mm == 0.0
    assert result.initial_abstraction_consumed_mm == pytest.approx(initial_abstraction)
    assert result.cumulative_infiltration_mm == 0.0


def test_crossing_initial_abstraction_puts_only_the_surplus_into_the_scs_equation() -> None:
    initial_abstraction = 0.2 * (25400.0 / 75.0 - 254.0)
    result = compute_excess((initial_abstraction + 1.0,), _CN_75)

    assert result.steps[0].excess_mm > 0.0
    assert result.steps[0].excess_mm < 1.0
    assert result.total_excess_mm == pytest.approx(0.0116731517509727, abs=_HAND_CALC_TOLERANCE_MM)


def test_cn_100_converts_every_millimetre_to_excess() -> None:
    result = compute_excess((10.0, 0.0, 3.0), ScsCnParameters(cn_ii=100, impervious_fraction=0.25))

    assert result.curve_number == 100.0
    assert result.potential_retention_mm == 0.0
    assert result.initial_abstraction_mm == 0.0
    assert [step.excess_mm for step in result.steps] == [10.0, 0.0, 3.0]
    assert result.initial_abstraction_consumed_mm == 0.0
    assert result.cumulative_infiltration_mm == 0.0
    assert result.runoff_coefficient == 1.0


def test_cn_100_remains_fully_excess_for_every_antecedent_moisture() -> None:
    rainfall = (4.0, 6.0)
    for amc in AmcCondition:
        result = compute_excess(rainfall, ScsCnParameters(cn_ii=100, amc_condition=amc))
        assert result.curve_number == 100.0
        assert result.potential_retention_mm == 0.0
        assert [step.excess_mm for step in result.steps] == [4.0, 6.0]


def test_small_curve_number_produces_no_excess_for_an_ordinary_storm() -> None:
    result = compute_excess((10.0, 20.0, 30.0, 40.0), ScsCnParameters(cn_ii=1))

    assert result.total_rainfall_mm == 100.0
    assert result.total_rainfall_mm < result.initial_abstraction_mm
    assert result.total_excess_mm == 0.0
    assert result.initial_abstraction_consumed_mm == 100.0
    assert result.cumulative_infiltration_mm == 0.0


def test_wetter_antecedent_moisture_produces_more_excess() -> None:
    totals = [
        compute_excess(_SIX_HOUR_MM, ScsCnParameters(cn_ii=75, amc_condition=amc)).total_excess_mm
        for amc in (AmcCondition.AMC_I, AmcCondition.AMC_II, AmcCondition.AMC_III)
    ]

    assert totals[0] < totals[1] < totals[2]


def test_antecedent_moisture_curve_numbers_follow_hawkins() -> None:
    dry = compute_excess(
        _SIX_HOUR_MM,
        ScsCnParameters(cn_ii=75, amc_condition=AmcCondition.AMC_I),
    )
    normal = compute_excess(_SIX_HOUR_MM, ScsCnParameters(cn_ii=75, amc_condition="AMC_II"))
    wet = compute_excess(
        _SIX_HOUR_MM,
        ScsCnParameters(cn_ii=75, amc_condition=AmcCondition.AMC_III),
    )

    assert dry.curve_number == pytest.approx(55.752212389381, abs=1e-9)
    assert normal.curve_number == 75.0
    assert wet.curve_number == pytest.approx(87.341772151899, abs=1e-9)
    assert dry.curve_number < normal.curve_number < wet.curve_number


def test_split_impervious_area_yields_excess_before_initial_abstraction() -> None:
    result = compute_excess(
        (5.0,),
        ScsCnParameters(cn_ii=75, impervious_fraction=0.5),
    )

    assert result.curve_number == 75.0
    assert result.initial_abstraction_mm > 5.0
    assert result.steps[0].excess_mm == 2.5
    assert result.initial_abstraction_consumed_mm == 2.5
    assert result.cumulative_infiltration_mm == 0.0


def test_split_impervious_area_weights_pervious_scs_excess() -> None:
    result = compute_excess(
        (10.0, 30.0),
        ScsCnParameters(cn_ii=75, impervious_fraction=0.5),
    )

    assert result.steps[0].excess_mm == 5.0
    assert result.steps[1].cumulative_rainfall_mm == pytest.approx(40.0)
    assert result.steps[1].cumulative_excess_mm == pytest.approx(
        22.469389438944, abs=_HAND_CALC_TOLERANCE_MM
    )
    assert result.steps[1].excess_mm == pytest.approx(17.469389438944, abs=_HAND_CALC_TOLERANCE_MM)
    assert result.cumulative_infiltration_mm == pytest.approx(
        9.063943894389, abs=_HAND_CALC_TOLERANCE_MM
    )
    assert result.initial_abstraction_consumed_mm == pytest.approx(8.466666666667, abs=1e-9)


def test_composite_cn_100_with_half_impervious_becomes_99() -> None:
    result = compute_excess(
        (10.0,),
        ScsCnParameters(
            cn_ii=100,
            impervious_fraction=0.5,
            impervious_treatment=ImperviousTreatment.COMPOSITE_CN,
        ),
    )

    assert result.curve_number == 99.0
    assert result.total_excess_mm < 10.0


def test_composite_curve_number_uses_the_tr55_adjustment_after_amc() -> None:
    adjusted = compute_excess(
        (30.0,),
        ScsCnParameters(
            cn_ii=75,
            impervious_fraction=0.2,
            impervious_treatment=ImperviousTreatment.COMPOSITE_CN,
        ),
    )
    wet_adjusted = compute_excess(
        (30.0,),
        ScsCnParameters(
            cn_ii=75,
            amc_condition=AmcCondition.AMC_III,
            impervious_fraction=0.2,
            impervious_treatment=ImperviousTreatment.COMPOSITE_CN,
        ),
    )

    assert adjusted.curve_number == 79.6
    assert adjusted.total_excess_mm == pytest.approx(3.513204550414, abs=_HAND_CALC_TOLERANCE_MM)
    assert wet_adjusted.curve_number == pytest.approx(89.473417721519, abs=1e-9)
    assert adjusted.total_excess_mm < wet_adjusted.total_excess_mm


def test_zero_lambda_starts_excess_immediately() -> None:
    result = compute_excess((10.0,), ScsCnParameters(cn_ii=75, lambda_coeff=0))

    assert result.initial_abstraction_mm == 0.0
    assert result.total_excess_mm == pytest.approx(1.056338028169, abs=_HAND_CALC_TOLERANCE_MM)


def test_revised_lambda_changes_the_initial_abstraction_threshold() -> None:
    result = compute_excess((10.0,), ScsCnParameters(cn_ii=75, lambda_coeff=0.05))

    assert result.initial_abstraction_mm == pytest.approx(4.233333333333, abs=1e-9)
    assert result.total_excess_mm == pytest.approx(0.367723307532, abs=_HAND_CALC_TOLERANCE_MM)


def test_lambda_one_holds_rainfall_until_retention_is_filled() -> None:
    result = compute_excess((10.0, 20.0, 30.0), ScsCnParameters(cn_ii=75, lambda_coeff=1))

    assert result.initial_abstraction_mm == pytest.approx(result.potential_retention_mm)
    assert result.total_rainfall_mm < result.initial_abstraction_mm
    assert result.total_excess_mm == 0.0


def test_excess_after_t_cut_does_not_restart_initial_abstraction() -> None:
    storm = (
        TimedRainfall("2026-07-01T00:00:00Z", 20.0),
        TimedRainfall("2026-07-01T01:00:00Z", 10.0),
        TimedRainfall("2026-07-01T02:00:00Z", 15.0),
        TimedRainfall("2026-07-01T03:00:00Z", 5.0),
    )
    t_cut = "2026-07-01T01:00:00Z"

    remainder = excess_after(storm, _CN_75, t_cut)
    restarted = compute_excess((15.0, 5.0), _CN_75)

    assert [step.time for step in remainder.steps] == [
        "2026-07-01T02:00:00Z",
        "2026-07-01T03:00:00Z",
    ]
    assert remainder.steps[0].cumulative_rainfall_mm == pytest.approx(45.0)
    assert remainder.remainder_excess_mm == pytest.approx(
        7.540151319683, abs=_HAND_CALC_TOLERANCE_MM
    )
    assert restarted.total_excess_mm == pytest.approx(0.107193515704, abs=_HAND_CALC_TOLERANCE_MM)
    assert remainder.remainder_excess_mm > restarted.total_excess_mm
    assert remainder.storm.total_excess_mm > remainder.remainder_excess_mm


def test_excess_after_compares_instants_and_keeps_labels() -> None:
    storm = (
        TimedRainfall("2026-07-01T00:00:00Z", 20.0),
        TimedRainfall("2026-07-01T01:00:00Z", 10.0),
        TimedRainfall("2026-07-01T02:00:00Z", 15.0),
        TimedRainfall("2026-07-01T03:00:00Z", 5.0),
    )

    by_z = excess_after(storm, _CN_75, "2026-07-01T01:00:00Z")
    by_offset = excess_after(storm, _CN_75, "2026-07-01T09:00:00+08:00")

    assert [step.time for step in by_offset.steps] == [step.time for step in by_z.steps]
    assert by_offset.remainder_excess_mm == pytest.approx(by_z.remainder_excess_mm)
    assert by_offset.t_cut == "2026-07-01T09:00:00+08:00"

    labelled = compute_excess((TimedRainfall("2026-07-01T08:00:00+08:00", 12.0),), _CN_75)
    assert labelled.steps[0].time == "2026-07-01T08:00:00+08:00"


def test_timestamps_are_ordered_as_instants() -> None:
    compute_excess(
        (
            TimedRainfall("2026-07-01T10:00:00+08:00", 5.0),
            TimedRainfall("2026-07-01T03:00:00Z", 5.0),
        ),
        _CN_75,
    )

    with pytest.raises(ScsCnError) as caught:
        compute_excess(
            (
                TimedRainfall("2026-07-01T03:00:00Z", 5.0),
                TimedRainfall("2026-07-01T10:00:00+08:00", 5.0),
            ),
            _CN_75,
        )

    assert caught.value.code is ScsCnErrorCode.UNORDERED_TIMESTAMPS


def test_water_balance_holds_across_varied_storms() -> None:
    rng = random.Random(20260601)
    curve_numbers = (1.0, 30.0, 55.0, 75.0, 90.0, 100.0)
    lambdas = (0.0, 0.05, 0.2, 1.0)
    fractions = (0.0, 0.25, 0.5, 1.0)

    for _ in range(40):
        length = rng.randint(1, 80)
        rainfall = [
            0.0 if rng.random() < 0.1 else rng.random() * 40.0 for _ in range(length)
        ]
        parameters = ScsCnParameters(
            cn_ii=rng.choice(curve_numbers),
            amc_condition=rng.choice(tuple(AmcCondition)),
            lambda_coeff=rng.choice(lambdas),
            impervious_fraction=rng.choice(fractions),
            impervious_treatment=rng.choice(tuple(ImperviousTreatment)),
        )
        result = compute_excess(rainfall, parameters)
        _assert_water_balance(result)


def test_zero_rainfall_has_no_runoff_coefficient() -> None:
    result = compute_excess((0.0, 0.0, 0.0), _CN_75)

    assert result.total_rainfall_mm == 0.0
    assert result.total_excess_mm == 0.0
    assert result.runoff_coefficient is None
    assert all(step.excess_mm == 0.0 and step.loss_mm == 0.0 for step in result.steps)


def test_batch_preserves_case_order_and_higher_cn_yields_more_excess() -> None:
    storm = (10.0, 25.0, 15.0)
    cases = tuple(
        ScsCnCase(storm, ScsCnParameters(cn_ii=curve_number)) for curve_number in (55, 75, 90)
    )

    results = compute_excess_batch(cases)

    assert len(results) == 3
    assert [result.total_rainfall_mm for result in results] == [50.0, 50.0, 50.0]
    assert results[0].total_excess_mm < results[1].total_excess_mm < results[2].total_excess_mm
    assert [result.curve_number for result in results] == [55.0, 75.0, 90.0]


@pytest.mark.parametrize("curve_number", [0, -5, 100.01, 101, float("nan"), float("inf")])
def test_rejects_illegal_curve_number(curve_number: float) -> None:
    with pytest.raises(ScsCnError) as caught:
        ScsCnParameters(cn_ii=curve_number)

    assert caught.value.code is ScsCnErrorCode.INVALID_CURVE_NUMBER
    assert "cn_ii" in caught.value.message


def test_rejects_boolean_curve_number() -> None:
    with pytest.raises(ScsCnError) as caught:
        ScsCnParameters(cn_ii=True)

    assert caught.value.code is ScsCnErrorCode.INVALID_CURVE_NUMBER


@pytest.mark.parametrize("ratio", [-0.01, 1.01, float("nan")])
def test_rejects_illegal_initial_abstraction_ratio(ratio: float) -> None:
    with pytest.raises(ScsCnError) as caught:
        ScsCnParameters(cn_ii=75, lambda_coeff=ratio)

    assert caught.value.code is ScsCnErrorCode.INVALID_LAMBDA


@pytest.mark.parametrize("fraction", [-0.01, 1.01])
def test_rejects_illegal_impervious_fraction(fraction: float) -> None:
    with pytest.raises(ScsCnError) as caught:
        ScsCnParameters(cn_ii=75, impervious_fraction=fraction)

    assert caught.value.code is ScsCnErrorCode.INVALID_IMPERVIOUS_FRACTION


def test_rejects_illegal_antecedent_moisture_name() -> None:
    with pytest.raises(ScsCnError) as caught:
        ScsCnParameters(cn_ii=75, amc_condition="wet")

    assert caught.value.code is ScsCnErrorCode.INVALID_AMC_CONDITION


@pytest.mark.parametrize(
    "rainfall",
    [(-1.0,), (1.0, -0.1, 2.0), (float("nan"),), (float("-inf"),)],
)
def test_rejects_negative_or_non_finite_rainfall(rainfall: tuple[float, ...]) -> None:
    with pytest.raises(ScsCnError) as caught:
        compute_excess(rainfall, _CN_75)

    assert caught.value.code is ScsCnErrorCode.INVALID_RAINFALL


def test_rejects_empty_rainfall() -> None:
    with pytest.raises(ScsCnError) as caught:
        compute_excess((), _CN_75)

    assert caught.value.code is ScsCnErrorCode.EMPTY_RAINFALL


@pytest.mark.parametrize(
    "stamp",
    [
        "2026-07-01T00:00:00",
        "2026-07-01 00:00:00Z",
        "not-a-time",
        "2026-13-01T00:00:00Z",
    ],
)
def test_rejects_timestamps_that_are_not_rfc3339(stamp: str) -> None:
    with pytest.raises(ScsCnError) as caught:
        compute_excess((TimedRainfall(stamp, 1.0),), _CN_75)

    assert caught.value.code is ScsCnErrorCode.INVALID_TIMESTAMP


def test_rejects_duplicate_timestamps() -> None:
    with pytest.raises(ScsCnError) as caught:
        compute_excess(
            (
                TimedRainfall("2026-07-01T00:00:00Z", 1.0),
                TimedRainfall("2026-07-01T00:00:00+00:00", 1.0),
            ),
            _CN_75,
        )

    assert caught.value.code is ScsCnErrorCode.UNORDERED_TIMESTAMPS


def test_excess_after_rejects_a_cut_that_is_not_rfc3339() -> None:
    storm = (TimedRainfall("2026-07-01T00:00:00Z", 5.0),)
    with pytest.raises(ScsCnError) as caught:
        excess_after(storm, _CN_75, "2026-07-01T00:00:00")

    assert caught.value.code is ScsCnErrorCode.INVALID_TIMESTAMP


def _assert_water_balance(result: ScsCnResult) -> None:
    summed_rain = math.fsum(step.rainfall_mm for step in result.steps)
    summed_excess = math.fsum(step.excess_mm for step in result.steps)
    summed_loss = math.fsum(step.loss_mm for step in result.steps)
    components = (
        result.total_excess_mm
        + result.initial_abstraction_consumed_mm
        + result.cumulative_infiltration_mm
    )

    assert abs(summed_rain - result.total_rainfall_mm) < _BALANCE_TOLERANCE_MM
    assert abs(summed_excess - result.total_excess_mm) < _BALANCE_TOLERANCE_MM
    assert abs(summed_excess + summed_loss - summed_rain) < _BALANCE_TOLERANCE_MM
    assert abs(components - result.total_rainfall_mm) < _BALANCE_TOLERANCE_MM

    previous_rain = 0.0
    previous_excess = 0.0
    for step in result.steps:
        assert step.excess_mm >= 0.0
        assert step.loss_mm >= 0.0
        assert step.excess_mm <= step.rainfall_mm + _BALANCE_TOLERANCE_MM
        assert step.cumulative_rainfall_mm + _BALANCE_TOLERANCE_MM >= previous_rain
        assert step.cumulative_excess_mm + _BALANCE_TOLERANCE_MM >= previous_excess
        previous_rain = step.cumulative_rainfall_mm
        previous_excess = step.cumulative_excess_mm

    if result.total_rainfall_mm == 0.0:
        assert result.runoff_coefficient is None
    else:
        assert result.runoff_coefficient == pytest.approx(
            result.total_excess_mm / result.total_rainfall_mm
        )

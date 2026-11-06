"""Shared hydrology domain logic: SCS-CN excess, curves, and level-pool routing."""

from hydro_agent.core.scs_cn import (
    AmcCondition,
    ExcessStep,
    ImperviousTreatment,
    RainfallSeries,
    RemainderExcess,
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

__all__ = [
    "AmcCondition",
    "ExcessStep",
    "ImperviousTreatment",
    "RainfallSeries",
    "RemainderExcess",
    "ScsCnCase",
    "ScsCnError",
    "ScsCnErrorCode",
    "ScsCnParameters",
    "ScsCnResult",
    "TimedRainfall",
    "compute_excess",
    "compute_excess_batch",
    "excess_after",
]

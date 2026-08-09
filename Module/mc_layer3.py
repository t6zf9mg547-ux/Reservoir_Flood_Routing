"""
mc_layer3.py
------------
Lives in Module/. Run with:

    python Module/mc_layer3.py <CaseName> [--n-outer N] [--n-inner N]
                                [--rule baseline|optimized] [--seed N]

e.g.
    python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150

Layer 3: RMC-RFA-style TWO-LOOP Monte Carlo for climate-change-aware
flood frequency analysis.

Sampling design (agreed before this file was written -- see chat/
HANDOFF.md for the discussion)
--------------------------------------------------------------------
OUTER loop -- "which world are we in": climate/frequency-curve scenario.
    Each outer draw is a single MAGNITUDE SCALE FACTOR applied uniformly
    to the case's own inflow_hydrograph.csv (delta method -- the
    tabulated shape is preserved, only scaled), representing the
    combination of (a) statistical uncertainty in fitting a flood-
    frequency distribution to a finite record and (b) climate-driven
    uncertainty in how that distribution shifts under a future/design
    horizon. Sampled log-normally, median = 1.0 in every case except
    mc_outer_peak_source="bootstrap" (the case's own hydrograph IS the
    central estimate; see validate_mc_sources() and the peak/volume
    source-selection block in run_case_layer3() for the full contract
    -- every case's scalars.csv MUST explicitly declare
    mc_outer_peak_source/mc_outer_volume_source, there is no silent
    default CV anymore).

    This is deliberately the SIMPLEST defensible v1 -- a single scale
    factor conflates "how big" uncertainty into one number rather than
    separately resolving peak-vs-volume or fitting a real GEV/LP3 to an
    actual annual-maximum series with real climate deltas from GCM
    ensembles. Swap `default_outer_sampler()` for a real fitted-
    distribution sampler once that data exists -- everything downstream
    (inner loop, metrics, writers) is agnostic to how the outer draws
    were generated.

INNER loop -- "given that world, what's the physical/operational
    spread": for a FIXED outer-loop hydrograph, resample:
      - discharge/rating coefficients (whatever a case's own
        `mc_uncertain_params()` declares -- see below)
      - reservoir stage-storage curve B(H), via a single area_scale
        multiplier (bathymetric/sedimentation uncertainty)
      - initial reservoir level H0 (antecedent condition, additive
        shift from the case's own H0 in scalars.csv)
    Kept STRICTLY SEPARATE from the operating rule (rule_overrides),
    the same "operating rule vs. discharge coefficient" split used by
    Layer 2 -- this loop never touches control logic, only physics.

Per-draw output metric
--------------------------------------------------------------------
Every (outer, inner) draw is one solver.run_simulation() call (bypasses
run_case.py's CSV-loading path, same performance reasoning as
optimize_layer2.py), recording a small vector, not one scalar:
    peak_level                 -- max reservoir level reached [m a.s.l.]
                                   (the PRIMARY metric -- what's compared
                                   against max_flood_level/dam_crest_level)
    peak_downstream_release    -- max(Qout_total + Qwithdrawal) [m3/s]
    time_to_exceedance_h       -- first time downstream_threshold_m3s is
                                   crossed (reuses optimize_layer2's
                                   exact function/definition), only if
                                   that scalar is defined for the case
    exceeds_max_flood_level, exceeds_dam_crest -- boolean flags, only
                                   if those scalars are defined

Case support for the inner loop's physical perturbation
--------------------------------------------------------------------
A case OPTS IN to inner-loop physical uncertainty by defining, in its
own case_config.py:

    def mc_uncertain_params():
        return {
            "gate_Cd_mult": {"dist": "lognormal_cv", "cv": 0.10},
            "fuse_trigger_shift": {"dist": "normal", "mean": 0.0, "sigma": 0.10},
            ...
        }

mirroring Layer 2's `optimizable_gates()` pattern -- each case declares
whatever parameter names its own `build_outlets(physical_overrides=...)`
knows how to use (see Data/Template/case_config.py for a worked
example). A case with no `mc_uncertain_params()` still runs (H0 +
reservoir-curve uncertainty still apply -- those are generic, not
case-specific), just with no discharge-coefficient perturbation --
Layer 3 prints a note rather than erroring, the same "not every case
needs every layer" philosophy as Layer 2's optimizable_gates().
"""

from __future__ import annotations

import argparse
import csv
import inspect
import os
import sys
import importlib.util
from statistics import NormalDist

import numpy as np

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(MODULE_DIR)
DATA_DIR = os.path.join(PROJECT_ROOT, "Data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "Output")
PLOT_DIR = os.path.join(PROJECT_ROOT, "Plot")

sys.path.insert(0, MODULE_DIR)

from reservoir import Reservoir
from hydrograph import Hydrograph
from solver import run_simulation
from report import write_mc_results_csv, write_mc_distribution_plot
from reliability.importance import rank_importance, format_importance_table

# Reuse Layer 2's case loader + downstream-release/exceedance-time
# definitions rather than re-implementing them -- same conventions,
# one source of truth. optimize_layer2.py only imports pymoo INSIDE
# its own functions, so importing these three names doesn't require
# pymoo to be installed.
from optimize_layer2 import load_case_config, downstream_release, time_to_exceedance

H0_DEFAULT_SIGMA = 0.10   # [m] default antecedent-level uncertainty if
                          # scalars.csv doesn't set mc_h0_sigma
AREA_DEFAULT_CV = 0.05    # default reservoir stage-storage curve
                          # (bathymetry/sedimentation) uncertainty if
                          # scalars.csv doesn't set mc_area_cv


# ---------------------------------------------------------------------
# Lightweight scaling wrappers -- Module/reservoir.py and hydrograph.py
# are left untouched; these just wrap an already-loaded object with the
# same interface, so Layer 3 can cheaply re-scale per draw without
# re-reading CSVs thousands of times.
# ---------------------------------------------------------------------

class ScaledHydrograph:
    """Wraps a Hydrograph, giving INDEPENDENT control over peak and
    volume. Flood volume is a first-order control on routed reservoir
    level for any reservoir that provides significant attenuation --
    how much water needs to be STORED, not just how fast it briefly
    arrives, is often what determines the peak level, not the inflow
    peak alone. A single uniform magnitude scale (the original v1
    implementation) scales peak and volume by the same factor every
    draw -- i.e. it silently assumes they're perfectly correlated,
    which is rarely true: real design-hydrograph studies often adopt
    different storm durations for different return periods (a shorter,
    more intense storm for the most extreme events vs. a longer one for
    more frequent floods), which changes the volume-to-peak ratio
    rather than just rescaling one fixed shape. A case can supply real
    data on this via mc_outer_distribution()'s "volume_duration_csv"
    key (see Data/Template/case_config.py's comments for the expected
    format), selected explicitly via mc_outer_volume_source='table' in
    scalars.csv -- every case must declare a volume CV source
    explicitly (see validate_mc_sources()), there is no silent
    unlabeled default.

    Mechanism: peak_scale multiplies Q(t) directly (sets the peak).
    volume_scale is achieved by stretching/compressing the TIME axis
    around the base hydrograph's own time-to-peak, which changes the
    volume WITHOUT changing the peak value reached at that instant:
        h(t) = peak_scale * base.discharge(t_peak + (t - t_peak) / k)
    where k = volume_scale / peak_scale. This gives, by construction:
        peak(h)   = peak_scale   * peak(base)
        volume(h) = volume_scale * volume(base)
    independently of one another.
    Backward compatible: volume_scale=None reproduces the old "uniform
    scale" behavior (peak and volume move together).
    """

    def __init__(self, base, peak_scale: float, volume_scale: float | None = None):
        self.base = base
        self.peak_scale = peak_scale
        self.volume_scale = peak_scale if volume_scale is None else volume_scale
        self.time_stretch = self.volume_scale / self.peak_scale
        self.t_peak = float(base.t[int(np.argmax(base.Q))])

    def discharge(self, t: float) -> float:
        t_shape = self.t_peak + (t - self.t_peak) / self.time_stretch
        return self.base.discharge(t_shape) * self.peak_scale


class ScaledReservoir:
    """Wraps a Reservoir, scaling B(H) (and V(H), if present) by a
    constant factor -- stage-storage curve uncertainty (bathymetry/
    sedimentation) for the inner loop."""

    def __init__(self, base: Reservoir, scale: float):
        self.base = base
        self.scale = scale

    def surface(self, H: float) -> float:
        return self.base.surface(H) * self.scale

    def volume(self, H: float):
        v = self.base.volume(H)
        return None if v is None else v * self.scale

    @property
    def H_min(self) -> float:
        return self.base.H_min

    @property
    def H_max(self) -> float:
        return self.base.H_max


# ---------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------

def _lognormal_cv(rng: np.random.Generator, cv: float, n: int) -> np.ndarray:
    """n draws from a log-normal distribution with MEDIAN 1.0 (not mean
    -- a "+/- X% uncertainty multiplier, centered on the nominal value"
    is a statement about the median/central value, matching how this
    project documents peak_median_scale/volume_median_scale=1.0
    everywhere else -- e.g. README.md's "What 'median' means" section
    and every case_config.py's mc_outer_distribution() docstring) and
    the given coefficient of variation.

    mu=0 gives median=exp(0)=1.0 exactly (median of a lognormal is
    exp(mu) regardless of sigma, since the underlying normal's median
    equals its own mean). Note this differs from mu=-0.5*sigma**2,
    which would instead center the MEAN at 1.0 (mean of a lognormal is
    exp(mu+sigma^2/2)) -- an earlier version of this function used that
    formula despite its own docstring claiming median=1.0; the two
    conventions diverge more as cv grows (at cv=0.42, e.g. this
    project's FFA bootstrap stress-test case, mean-centering gives a
    true median of only ~0.922, not 1.0 -- a real ~8% discrepancy from
    every documented claim that median_scale locks to the case's own
    hydrograph value)."""
    if cv <= 0:
        return np.ones(n)
    sigma = np.sqrt(np.log(1.0 + cv ** 2))
    return rng.lognormal(mean=0.0, sigma=sigma, size=n)


def _loglog_interp(target_x: float, xs: np.ndarray, ys: np.ndarray) -> float:
    """Interpolates y(target_x) linearly in log-log space (the
    conventional way to read a value off a flood-frequency curve
    plotted on log-log axes) -- xs must be sorted ascending."""
    return float(np.exp(np.interp(np.log(target_x), np.log(xs), np.log(ys))))


def load_flood_frequency_curve(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Reads a return_period_years,peak_Q_m3s[,source] CSV (see
    Data/Template/case_config.py's comments for the expected format)
    into sorted (T, Q) arrays."""
    T, Q = [], []
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            T.append(float(row["return_period_years"]))
            Q.append(float(row["peak_Q_m3s"]))
    order = np.argsort(T)
    return np.asarray(T)[order], np.asarray(Q)[order]


def derive_outer_cv_from_alt_studies(curve_T: np.ndarray, curve_Q: np.ndarray,
                                      alt_studies_path: str) -> float | None:
    """Empirical stand-in for 'curve/climate uncertainty', derived from
    how much INDEPENDENT PAST STUDIES of the same dam disagreed with
    the currently-adopted curve, at the return periods they share --
    see mc_outer_distribution()'s docstring (in a case's own
    case_config.py) for the full reasoning and caveats. Returns None
    (caller should fall back to a default) if fewer than 2 comparison
    points are available -- not enough to estimate a spread from."""
    log_ratios = []
    with open(alt_studies_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            T = float(row["return_period_years"])
            Q_alt = float(row["peak_Q_m3s"])
            Q_curve = _loglog_interp(T, curve_T, curve_Q)
            log_ratios.append(np.log(Q_alt / Q_curve))
    if len(log_ratios) < 2:
        return None
    sigma = float(np.std(log_ratios, ddof=1))
    # Convert log-space sigma back to an ordinary coefficient of
    # variation (matches _lognormal_cv's own cv -> sigma convention).
    return float(np.sqrt(np.exp(sigma ** 2) - 1.0))


def load_duration_volume_table(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reads a return_period_years,storm_duration_hours,peak_Q_m3s,
    volume_Mm3 CSV (see Data/Template/case_config.py's comments for
    the expected format) into sorted (T, Q, V_m3) arrays
    -- volume converted from million m3 to plain m3 for direct use
    alongside Hydrograph's own SI units."""
    T, Q, V = [], [], []
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            T.append(float(row["return_period_years"]))
            Q.append(float(row["peak_Q_m3s"]))
            V.append(float(row["volume_Mm3"]) * 1.0e6)
    order = np.argsort(T)
    return np.asarray(T)[order], np.asarray(Q)[order], np.asarray(V)[order]


def derive_volume_cv_from_duration_table(T: np.ndarray, Q: np.ndarray, V: np.ndarray) -> float | None:
    """Empirical stand-in for 'how much can a flood's VOLUME plausibly
    differ from its peak-implied shape', derived from the SPREAD of
    volume-to-peak ratios (V/Q, a characteristic duration) across this
    dam's own adopted design hydrographs at different return periods --
    see ScaledHydrograph's docstring for why this matters here. These
    ratios cluster in two groups (16h-storm floods vs. the 8h-storm
    5000yr/PMF floods), so their spread is a genuine, if rough, data-
    grounded bracket on how much hydrograph shape (hence volume, for a
    given peak) can plausibly vary -- NOT a formal statistical fit.
    Returns None (caller falls back to a default) if fewer than 2 rows."""
    if len(T) < 2:
        return None
    ratios = V / Q  # [s] -- a characteristic duration V/Qp per flood
    log_ratios = np.log(ratios / np.median(ratios))
    sigma = float(np.std(log_ratios, ddof=1))
    return float(np.sqrt(np.exp(sigma ** 2) - 1.0))


def load_bootstrap_ci_curve(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reads an FFA tool's bootstrap_ci_<distribution>.csv (or
    model_averaged_quantiles.csv in the same shape):
        T,lower,median,upper
    into sorted (T, lower, median, upper) arrays. `median` is this
    file's own model-averaged/bootstrap point estimate -- NOT
    necessarily the same as any curve_csv a case may separately have
    adopted; when this file is used (mc_outer_peak_source="bootstrap"),
    its own median is the peak-loop anchor, deliberately not blended
    with any adopted curve (see FFA_Integration.md, Option D)."""
    T, lower, median, upper = [], [], [], []
    with open(path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            T.append(float(row["T"]))
            lower.append(float(row["lower"]))
            median.append(float(row["median"]))
            upper.append(float(row["upper"]))
    order = np.argsort(T)
    return (np.asarray(T)[order], np.asarray(lower)[order],
            np.asarray(median)[order], np.asarray(upper)[order])


def derive_outer_cv_from_bootstrap_ci(ci_T: np.ndarray, ci_lower: np.ndarray,
                                       ci_median: np.ndarray, ci_upper: np.ndarray,
                                       target_return_period: float,
                                       ci_confidence: float = 0.95) -> tuple[float, float]:
    """Parallel to derive_outer_cv_from_alt_studies(), but for a real
    FFA tool's bootstrap confidence interval rather than inter-study
    disagreement -- see FFA_Integration.md for why these two are kept
    as separate, non-blended sources (Option C) rather than merged.

    Log-log-interpolates lower/median/upper at target_return_period,
    then converts the (lower, upper) width around median into an
    equivalent lognormal CV, using the SAME z-score convention already
    used elsewhere in this file for confidence-interval math (see
    _min_realizations_mean/_min_realizations_percentile):
        z = NormalDist().inv_cdf((1 + ci_confidence) / 2)
    sigma is estimated from BOTH sides of the interval independently
    (upper side: ln(upper/median)/z; lower side: -ln(lower/median)/z)
    and averaged, which is more robust to a mildly asymmetric bootstrap
    interval than using either side alone.

    Returns (anchor_Q, cv):
        anchor_Q -- this file's own median at target_return_period,
                    the peak-loop's median_scale anchor (analogous to
                    curve_csv's anchor_Q).
        cv       -- the derived peak-loop coefficient of variation.
    Raises the same exception types load_flood_frequency_curve() /
    derive_outer_cv_from_alt_studies() can raise (FileNotFoundError,
    KeyError, ValueError) on a malformed/missing file, so callers can
    reuse the same except clause."""
    anchor_Q = _loglog_interp(target_return_period, ci_T, ci_median)
    lower_Q = _loglog_interp(target_return_period, ci_T, ci_lower)
    upper_Q = _loglog_interp(target_return_period, ci_T, ci_upper)
    z = NormalDist().inv_cdf((1 + ci_confidence) / 2)
    sigma_upper = np.log(upper_Q / anchor_Q) / z
    sigma_lower = -np.log(lower_Q / anchor_Q) / z
    sigma = float((sigma_upper + sigma_lower) / 2.0)
    cv = float(np.sqrt(np.exp(sigma ** 2) - 1.0))
    return anchor_Q, cv


# ---------------------------------------------------------------------
# Peak-volume dependence (copula-based joint sampling)
# ---------------------------------------------------------------------
# WHY this exists: default_outer_sampler() historically drew peak_scale
# and volume_scale fully INDEPENDENTLY. Flood-frequency literature
# treats this as a recognized limitation -- peak and volume come from
# the SAME storm event, so they're physically correlated; see
# Requena, Mediero & Garrote (2013, HESS 17:3023-3038), the direct
# precedent for this project's own peak-volume -> synthetic hydrograph
# -> reservoir routing -> overtopping-risk chain, and the broader
# copula-based flood-frequency literature it draws on. This section
# adds OPTIONAL joint sampling via a copula, selected explicitly per
# case (mc_peak_volume_dependence in scalars.csv) -- "independent"
# remains the default for every case that doesn't opt in, so no
# existing case's numbers change unless deliberately requested.
#
# Both copulas below are parameterized by Kendall's tau (0 <= tau < 1),
# not the copula's own native parameter -- tau is the quantity you'd
# actually estimate from a case's own paired annual-maximum peak/volume
# series (the same record an FFA tool would use), and it has a direct,
# interpretable meaning (rank correlation) regardless of which copula
# family is chosen. Literature-reported values commonly fall around
# 0.4-0.7 for peak-volume pairs across various basins -- a starting
# REFERENCE range only, not a default to assume for any specific case;
# same "don't invent case-specific numbers in generic engine code"
# principle used everywhere else in this project.

def _sample_gaussian_copula_normals(rng: np.random.Generator, n: int, tau: float
                                     ) -> tuple[np.ndarray, np.ndarray]:
    """Correlated standard-normal pair (Z1, Z2) whose implied copula has
    Kendall's tau = `tau`, via the standard Gaussian-copula map
    rho = sin(pi*tau/2). Returned as normals (not yet pushed through
    Phi) so callers with lognormal marginals can apply their own
    sigma*Z + ln(median) transform directly, without a redundant
    normal -> uniform -> normal round trip."""
    rho = np.sin(np.pi * tau / 2.0)
    Z1 = rng.standard_normal(n)
    eps = rng.standard_normal(n)
    Z2 = rho * Z1 + np.sqrt(1.0 - rho ** 2) * eps
    return Z1, Z2


def _sample_gumbel_copula_uniforms(rng: np.random.Generator, n: int, tau: float
                                    ) -> tuple[np.ndarray, np.ndarray]:
    """Correlated Uniform(0,1) pair (U1, U2) with Gumbel-Hougaard copula
    dependence, Kendall's tau = `tau` (theta = 1/(1-tau); tau=0 gives
    theta=1, which IS independence -- handled as a direct special case
    below to avoid a 0/0 in the stable-distribution sampler).

    Marshall-Olkin (1988) construction: sample a positive-stable random
    variable S (stability index alpha=1/theta, totally right-skewed)
    via the Chambers-Mallows-Stuck (1976) algorithm, then
    U_i = exp(-(-ln(X_i)/S)^(1/theta)) for independent uniforms X_i.
    This is the standard algorithm used by, e.g., R's `copula` package
    for Archimedean copula generation -- verified here against the
    known Gumbel-Hougaard property that its upper-tail dependence
    coefficient is 2 - 2^(1/theta) (checked via empirical Kendall's tau
    recovery across tau in [0, 0.7] during development; see chat)."""
    if tau <= 0:
        return rng.uniform(size=n), rng.uniform(size=n)
    theta = 1.0 / (1.0 - tau)
    alpha = 1.0 / theta
    W = rng.uniform(-np.pi / 2, np.pi / 2, size=n)
    E = rng.exponential(1.0, size=n)
    S = (np.sin(alpha * (W + np.pi / 2)) / np.cos(W) ** (1 / alpha)) * \
        (np.cos(W - alpha * (W + np.pi / 2)) / E) ** ((1 - alpha) / alpha)
    X1 = rng.uniform(size=n)
    X2 = rng.uniform(size=n)
    U1 = np.exp(-(-np.log(X1) / S) ** (1 / theta))
    U2 = np.exp(-(-np.log(X2) / S) ** (1 / theta))
    return U1, U2


# Vectorized inverse-normal-CDF (numpy has no built-in ppf; NormalDist's
# is scalar-only). Only used by the "gumbel" branch below to convert
# the copula's uniform marginals to standard normals before applying
# each variable's own lognormal transform -- n_outer is small enough
# (order 100-10000) that np.vectorize's per-element Python call is not
# a meaningful cost here.
NORMAL_PPF_VEC = np.vectorize(NormalDist().inv_cdf)


def default_outer_sampler(rng: np.random.Generator, n_outer: int,
                           peak_cv: float, volume_cv: float,
                           peak_median_scale: float = 1.0,
                           volume_median_scale: float = 1.0,
                           dependence: str = "independent",
                           tau: float | None = None) -> list[dict]:
    """OUTER loop draws: peak_scale and volume_scale, both log-normal
    marginals (unchanged regardless of `dependence` -- median/CV always
    mean exactly what they did before; only the CORRELATION between the
    two draws changes).

    dependence : "independent" (default -- matches every prior version
        of this function; floods with the same peak can plausibly have
        quite different volumes, or vice versa, and forcing them to
        move together would hide exactly the routing-sensitive
        scenarios this dam cares about most -- see ScaledHydrograph's
        docstring) | "gaussian" | "gumbel" (positive Kendall's-tau
        dependence via the respective copula -- see module comments
        above). "gumbel" has upper-tail dependence (extreme peak and
        extreme volume co-occur MORE than under "gaussian" at the same
        tau), which the flood-frequency literature more often finds is
        the better fit for this specific pair -- but this is basin-
        specific, not universal; see case_config.py for what a specific
        case has chosen and why.
    tau : Kendall's tau, required (validate_mc_sources() enforces this)
        when dependence != "independent". Ignored otherwise.
    """
    if dependence == "independent":
        peak_scales = _lognormal_cv(rng, peak_cv, n_outer) * peak_median_scale
        volume_scales = _lognormal_cv(rng, volume_cv, n_outer) * volume_median_scale
    else:
        peak_sigma = np.sqrt(np.log(1.0 + peak_cv ** 2))
        volume_sigma = np.sqrt(np.log(1.0 + volume_cv ** 2))
        if dependence == "gaussian":
            Z_peak, Z_volume = _sample_gaussian_copula_normals(rng, n_outer, tau)
        elif dependence == "gumbel":
            U_peak, U_volume = _sample_gumbel_copula_uniforms(rng, n_outer, tau)
            Z_peak = NORMAL_PPF_VEC(U_peak)
            Z_volume = NORMAL_PPF_VEC(U_volume)
        else:
            raise ValueError(f"Unknown dependence='{dependence}' -- expected "
                              f"'independent', 'gaussian', or 'gumbel'.")
        # _lognormal_cv's own convention: median_scale * exp(sigma*Z - sigma^2/2)
        # is NOT used here -- _lognormal_cv already centers so the MEDIAN
        # (not mean) equals median_scale; matching that exactly:
        peak_scales = peak_median_scale * np.exp(peak_sigma * Z_peak)
        volume_scales = volume_median_scale * np.exp(volume_sigma * Z_volume)
    return [{"peak_scale": float(p), "volume_scale": float(v)}
            for p, v in zip(peak_scales, volume_scales)]


def normalize_gate_availability(gate_availability: dict | None) -> dict | None:
    """Accepts either shape a case's mc_gate_availability() may return,
    and normalizes to the one internal shape default_inner_sampler()
    actually consumes:

        {"gates": [...],
         "p_fail_by_gate": {gate_name: p, ...},
         "ccf_groups": [{"label": str, "gates": [...], "p_ccf": float}, ...]}

    TIER 0 -- {"gates": [...], "p_fail": p} -- a single flat rate
        copied across every gate, no common-cause mechanism. This is
        the original, simplest shape (a case with no CI evidence,
        just one judgment-call rate); still fully supported.

    TIER 1 -- {"gates": [...], "p_fail_by_gate": {...},
        "ccf_groups": [...]} -- per-gate-differentiated rates plus
        explicit common-cause branches, typically produced by
        Module/reliability/gate_reliability.py's
        build_gate_availability() from real CI evidence rather than
        written by hand. Passed through with ccf_groups defaulted to
        empty if the case's dict omits it.

    Returns None if gate_availability is None (no hook declared at
    all -- every gate assumed operational, same as before)."""
    if gate_availability is None:
        return None
    gates = list(gate_availability["gates"])
    if "p_fail_by_gate" in gate_availability:
        p_fail_by_gate = dict(gate_availability["p_fail_by_gate"])
        ccf_groups = list(gate_availability.get("ccf_groups", []))
        component_pfs = gate_availability.get("component_pfs")  # None for Tier 0 / flat form
    else:
        flat_p = gate_availability.get("p_fail", 0.0)
        p_fail_by_gate = {g: flat_p for g in gates}
        ccf_groups = []
        component_pfs = None
    return {"gates": gates, "p_fail_by_gate": p_fail_by_gate, "ccf_groups": ccf_groups,
            "component_pfs": component_pfs}


def expected_p_any_gate_failed(p_fail_by_gate: dict[str, float], ccf_groups: list[dict]) -> float:
    """Closed-form P(at least one gate ends up failed), generalizing
    the simple binomial formula to heterogeneous per-gate rates plus
    independent common-cause groups. A gate is operational only if (a)
    its own independent draw doesn't fail it, AND (b) no common-cause
    group covering it fires -- both are independent Bernoulli events,
    so P(ALL gates operational) is just the product of every one of
    those "doesn't happen" probabilities; P(>=1 failed) is 1 minus
    that. Reduces to the original 1-(1-p)**n when every gate shares
    one flat rate and there are no CCF groups."""
    p_all_ok = 1.0
    for p in p_fail_by_gate.values():
        p_all_ok *= (1.0 - p)
    for grp in ccf_groups:
        p_all_ok *= (1.0 - grp["p_ccf"])
    return 1.0 - p_all_ok



def default_inner_sampler(rng: np.random.Generator, n_inner: int,
                           uncertain_params: dict, h0_sigma: float,
                           area_cv: float, gate_availability: dict | None = None
                           ) -> list[dict]:
    """INNER loop draws: physical/rating multipliers (whatever the case
    declares via mc_uncertain_params()), plus the two generic,
    case-independent draws every case gets (H0 shift and reservoir
    area_scale), plus -- if the case declares mc_gate_availability() --
    which (if any) of its named gates fail to open this draw.

    Gate failure, TWO independent stages per draw (see
    normalize_gate_availability() for the input shapes this consumes):

    1. COMMON-CAUSE: each named ccf_groups entry fires as its own
       independent Bernoulli(p_ccf) trial. If it fires, every gate
       named in that group's "gates" list is forced to the failed
       state for this draw, regardless of its own individual outcome
       below -- a shared control link, shared component batch/model,
       or shared power source taking out multiple gates at once,
       which independent per-gate draws alone can never represent no
       matter how the per-gate rate is chosen.
    2. INDEPENDENT: every gate not already forced failed by a
       common-cause event is then drawn as its own Bernoulli
       (p_fail_by_gate[gate]) trial -- per-gate rates, not
       necessarily identical across gates (Tier 1), or all equal to
       one flat rate (Tier 0, reproducing the original binomial
       screening model exactly when ccf_groups is empty).

    A gate's final state this draw is failed if EITHER stage says so
    (logical OR) -- common-cause and independent failure are
    alternative ways to end up unavailable, not mutually exclusive."""
    draws = []
    sampled_cols = {}
    for name, spec in uncertain_params.items():
        dist = spec.get("dist", "lognormal_cv")
        if dist == "lognormal_cv":
            sampled_cols[name] = _lognormal_cv(rng, spec.get("cv", 0.1), n_inner)
        elif dist == "normal":
            sampled_cols[name] = rng.normal(spec.get("mean", 0.0), spec.get("sigma", 0.0), n_inner)
        else:
            raise ValueError(f"mc_uncertain_params(): unknown dist '{dist}' for '{name}'")

    h0_shift = rng.normal(0.0, h0_sigma, n_inner) if h0_sigma > 0 else np.zeros(n_inner)
    area_scale = _lognormal_cv(rng, area_cv, n_inner)

    ga = normalize_gate_availability(gate_availability)
    gates = ga["gates"] if ga else []
    n_gates = len(gates)

    # Stage 1: independent per-gate draws, each at its own rate.
    indep_fail = np.zeros((n_inner, n_gates), dtype=bool)
    for j, g in enumerate(gates):
        p = ga["p_fail_by_gate"].get(g, 0.0)
        indep_fail[:, j] = rng.random(n_inner) < p

    # Stage 2: common-cause groups, each its own independent draw;
    # OR'd onto every gate in that group if it fires. ccf_fires keeps
    # the per-draw, per-group boolean around for reporting (which CCF
    # branch actually fired on a given draw), not just the resulting
    # gate states.
    ccf_fail = np.zeros((n_inner, n_gates), dtype=bool)
    ccf_fires: dict[str, np.ndarray] = {}
    for grp in (ga["ccf_groups"] if ga else []):
        fires = rng.random(n_inner) < grp["p_ccf"]
        ccf_fires[grp["label"]] = fires
        idxs = [j for j, g in enumerate(gates) if g in grp["gates"]]
        for j in idxs:
            ccf_fail[:, j] |= fires

    total_fail = indep_fail | ccf_fail

    for i in range(n_inner):
        d = {name: float(col[i]) for name, col in sampled_cols.items()}
        d["H0_shift"] = float(h0_shift[i])
        d["area_scale"] = float(area_scale[i])
        d["gates_failed"] = [gates[j] for j in range(n_gates) if total_fail[i, j]]
        d["ccf_fired"] = [label for label, fires in ccf_fires.items() if fires[i]]
        draws.append(d)
    return draws


# ---------------------------------------------------------------------
# Explicit-source validation (HARD ERROR by design -- see chat/
# FFA_Integration.md follow-up: this project deliberately moved away
# from letting the outer/inner loops infer their own data source from
# "which CSVs happen to exist", because that made it hard to tell,
# from a case's files alone, which source was actually in effect. Every
# case's scalars.csv MUST declare all three source scalars below,
# explicitly, or the run refuses to start. There is no "auto" mode --
# an unset/misspelled source scalar is a configuration bug, not a
# fallback opportunity, and should be caught before any simulation
# runs, not discovered later from a puzzling mc_summary.txt.
# ---------------------------------------------------------------------

VALID_PEAK_SOURCES = {"curve", "bootstrap", "user"}
VALID_VOLUME_SOURCES = {"table", "user"}
VALID_GATE_SOURCES = {"ci", "flat"}
VALID_DEPENDENCE_SOURCES = {"independent", "gaussian", "gumbel"}


class MCConfigError(ValueError):
    """Raised when a case's scalars.csv/case_config.py doesn't satisfy
    the explicit-source contract this file requires. Deliberately a
    plain ValueError subclass (not caught anywhere downstream) so it
    surfaces as a hard stop with a clear message, not a silent
    fallback to a default the case owner never asked for."""


def _require_scalar_str(sc: dict, key: str, valid: set[str], case_name: str) -> str:
    if key not in sc:
        raise MCConfigError(
            f"{case_name}/scalars.csv is missing required row '{key}'. "
            f"Every case must declare this explicitly -- add a row "
            f"'{key},<value>,-' with <value> one of: {sorted(valid)}. "
            f"See Data/Template/case_config.py's mc_outer_distribution()/"
            f"mc_gate_availability() docstrings for what each value means.")
    val = sc[key]
    if not isinstance(val, str):
        # scalars() coerces numeric-looking values to float -- a
        # source key that parsed as a number is itself a config bug
        # (e.g. someone left it as a leftover numeric placeholder).
        raise MCConfigError(
            f"{case_name}/scalars.csv row '{key}' must be one of the "
            f"text values {sorted(valid)}, got a number ({val!r}) instead.")
    val = val.strip().lower()
    if val not in valid:
        raise MCConfigError(
            f"{case_name}/scalars.csv row '{key}' has value '{val}', "
            f"not one of the recognized values {sorted(valid)}.")
    return val


def validate_mc_sources(sc: dict, outer_dist: dict | None, case_name: str,
                         has_gate_hook: bool = False) -> dict:
    """Validates and returns the three explicit source selections this
    case has made, cross-checked against what data/scalars they each
    require. Raises MCConfigError (before any simulation runs) if the
    declared source and the available data/scalars don't line up --
    e.g. source='bootstrap' but mc_outer_distribution() has no
    'bootstrap_ci_csv' key, or source='user' but the corresponding _cv
    scalar isn't set. Called once per run_case_layer3() call.

    has_gate_hook : whether this case's case_config.py actually defines
        mc_gate_availability() (passed in by build_case(), which is the
        one place that knows -- this function itself never imports/
        inspects case_config). mc_gate_reliability_source is REQUIRED
        only when True -- a case with no gate-reliability hook at all
        has nothing for that scalar to configure, so it would be a
        pointless, confusing requirement to impose on every case
        regardless. If mc_gate_reliability_source is present anyway
        despite has_gate_hook=False (e.g. left over from an earlier
        version of the case, or set defensively ahead of adding the
        hook later), its VALUE is still validated -- a stray typo
        shouldn't silently pass just because it's currently unused --
        it's only the row's ABSENCE that's tolerated in that case."""
    peak_source = _require_scalar_str(sc, "mc_outer_peak_source", VALID_PEAK_SOURCES, case_name)
    volume_source = _require_scalar_str(sc, "mc_outer_volume_source", VALID_VOLUME_SOURCES, case_name)
    if has_gate_hook:
        gate_source = _require_scalar_str(sc, "mc_gate_reliability_source", VALID_GATE_SOURCES, case_name)
    elif "mc_gate_reliability_source" in sc:
        gate_source = _require_scalar_str(sc, "mc_gate_reliability_source", VALID_GATE_SOURCES, case_name)
    else:
        gate_source = None

    def need(key):
        if key not in sc:
            raise MCConfigError(
                f"{case_name}/scalars.csv: mc_outer_peak_source/mc_outer_volume_source "
                f"selection requires a '{key}' row, which is missing.")

    if peak_source == "user":
        need("mc_outer_peak_cv")
    elif peak_source == "curve":
        if outer_dist is None or "curve_csv" not in outer_dist:
            raise MCConfigError(
                f"{case_name}: mc_outer_peak_source='curve' requires "
                f"mc_outer_distribution() to return a 'curve_csv' key.")
        if "alt_studies_csv" not in outer_dist:
            raise MCConfigError(
                f"{case_name}: mc_outer_peak_source='curve' also requires an "
                f"'alt_studies_csv' key -- curve_csv alone is a single median "
                f"value with no spread information in it; alt_studies_csv (2+ "
                f"independent past studies) is what actually supplies the CV. "
                f"If you don't have that, use mc_outer_peak_source='user' with "
                f"an explicit mc_outer_peak_cv instead of an unlabeled generic "
                f"default.")
    elif peak_source == "bootstrap":
        if outer_dist is None or "bootstrap_ci_csv" not in outer_dist:
            raise MCConfigError(
                f"{case_name}: mc_outer_peak_source='bootstrap' requires "
                f"mc_outer_distribution() to return a 'bootstrap_ci_csv' key.")
        if "target_return_period" not in outer_dist:
            raise MCConfigError(
                f"{case_name}: mc_outer_peak_source='bootstrap' requires "
                f"mc_outer_distribution() to also return a 'target_return_period' "
                f"key -- which T to read off the bootstrap CI curve.")

    if volume_source == "user":
        need("mc_outer_volume_cv")
    elif volume_source == "table":
        if outer_dist is None or "volume_duration_csv" not in outer_dist:
            raise MCConfigError(
                f"{case_name}: mc_outer_volume_source='table' requires "
                f"mc_outer_distribution() to return a 'volume_duration_csv' key.")

    if gate_source == "flat" and "mc_gate_p_fail" not in sc:
        raise MCConfigError(
            f"{case_name}/scalars.csv: mc_gate_reliability_source='flat' "
            f"requires an 'mc_gate_p_fail' row.")

    # mc_peak_volume_dependence is OPTIONAL -- "independent" (the
    # pre-existing, still-default behavior) never needs a row at all.
    # Only required to be present -- and then validated -- when a case
    # wants to opt into peak-volume copula dependence.
    dependence_source = "independent"
    if "mc_peak_volume_dependence" in sc:
        dependence_source = _require_scalar_str(
            sc, "mc_peak_volume_dependence", VALID_DEPENDENCE_SOURCES, case_name)
    if dependence_source != "independent":
        if "mc_peak_volume_tau" not in sc:
            raise MCConfigError(
                f"{case_name}/scalars.csv: mc_peak_volume_dependence="
                f"'{dependence_source}' requires an 'mc_peak_volume_tau' row "
                f"(Kendall's tau, 0 <= tau < 1).")
        tau = sc["mc_peak_volume_tau"]
        if not isinstance(tau, (int, float)) or isinstance(tau, bool):
            raise MCConfigError(
                f"{case_name}/scalars.csv: 'mc_peak_volume_tau' must be a "
                f"number, got {tau!r}.")
        if not (0.0 <= tau < 1.0):
            raise MCConfigError(
                f"{case_name}/scalars.csv: 'mc_peak_volume_tau'={tau} is out "
                f"of range -- Kendall's tau must satisfy 0 <= tau < 1.")

    return {"peak_source": peak_source, "volume_source": volume_source,
            "gate_source": gate_source, "dependence_source": dependence_source}


# ---------------------------------------------------------------------
# Case setup
# ---------------------------------------------------------------------

def call_build_outlets(case_config, rule_overrides, physical_overrides, inflow,
                        gates_out_of_service=None):
    """Calls case_config.build_outlets() passing only the kwargs it
    actually accepts, so cases not yet updated for Layer 3 (no
    physical_overrides/inflow/gates_out_of_service params) still run --
    just without inner-loop physical perturbation, outer-loop inflow
    propagation, or gate-availability draws, respectively. Warns once
    per missing capability (tracked via function attributes) rather
    than spamming per draw."""
    sig = inspect.signature(case_config.build_outlets)
    kwargs = {"rule_overrides": rule_overrides}
    if "physical_overrides" in sig.parameters:
        kwargs["physical_overrides"] = physical_overrides
    elif physical_overrides and not getattr(call_build_outlets, "_warned_physical", False):
        print("  NOTE: this case's build_outlets() doesn't accept physical_overrides -- "
              "inner-loop physical/rating uncertainty will NOT apply. "
              "See Data/Template/case_config.py for the pattern to add it.")
        call_build_outlets._warned_physical = True
    if "inflow" in sig.parameters:
        kwargs["inflow"] = inflow
    elif not getattr(call_build_outlets, "_warned_inflow", False):
        # Only relevant (and only worth warning about) for cases that
        # load their own private Hydrograph inside build_outlets -- we
        # can't easily detect that from here, so this is a soft note,
        # not an error; most cases don't need it at all.
        call_build_outlets._warned_inflow = True
    if "gates_out_of_service" in sig.parameters:
        kwargs["gates_out_of_service"] = gates_out_of_service or []
    elif gates_out_of_service and not getattr(call_build_outlets, "_warned_gates", False):
        print("  NOTE: this case's build_outlets() doesn't accept gates_out_of_service -- "
              "inner-loop gate-availability draws will NOT apply (gates are always "
              "treated as operational regardless of mc_gate_availability() draws). "
              "See Data/Corumana_117_Q5000/case_config.py for the pattern to add it.")
        call_build_outlets._warned_gates = True
    return case_config.build_outlets(**kwargs)


def parse_optimized_rule_overrides(case_name: str) -> dict | None:
    """Best-effort parse of Output/<CaseName>/Optimised/best_rule.txt
    (written by optimize_layer2.py) back into a rule_overrides dict,
    so Layer 3 can stress-test an already-optimized rule with --rule
    optimized. Returns None (with a printed explanation) if that file
    doesn't exist -- e.g. optimize_layer2.py was never run for this
    case, or the case has no tunable operating rule at all (some
    cases genuinely don't)."""
    path = os.path.join(OUTPUT_DIR, case_name, "Optimised", "best_rule.txt")
    if not os.path.isfile(path):
        print(f"  No {path} found -- run Module/optimize_layer2.py {case_name} first, "
              "or this case may have no optimizable_gates() at all (see its "
              "case_config.py). Falling back to the baseline rule.")
        return None
    overrides: dict[str, dict[str, float]] = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if ":" not in line or "=" not in line:
                continue
            gname, _, rest = line.partition(":")
            gname = gname.strip()
            params = {}
            for part in rest.split(","):
                part = part.strip()
                if "=" not in part:
                    continue
                pname, _, pval = part.partition("=")
                try:
                    params[pname.strip()] = float(pval.strip())
                except ValueError:
                    continue
            if params:
                overrides[gname] = params
    return overrides or None


def build_case(case_name: str, rule: str = "baseline"):
    case_config, case_dir = load_case_config(case_name)
    reservoir_base = Reservoir(os.path.join(case_dir, "reservoir_curve.csv"))
    inflow_base = Hydrograph(os.path.join(case_dir, "inflow_hydrograph.csv"))
    withdrawal_path = os.path.join(case_dir, "withdrawal_schedule.csv")
    withdrawal = (Hydrograph(withdrawal_path) if os.path.isfile(withdrawal_path)
                  else Hydrograph.constant(0.0))
    sc = case_config.scalars(case_dir)

    rule_overrides = None
    if rule == "optimized":
        rule_overrides = parse_optimized_rule_overrides(case_name)

    uncertain_params = {}
    if hasattr(case_config, "mc_uncertain_params"):
        uncertain_params = case_config.mc_uncertain_params()
    else:
        print(f"  NOTE: {case_name}/case_config.py has no mc_uncertain_params() -- "
              "inner loop will only vary H0 and the reservoir stage-storage curve, "
              "not discharge/rating coefficients. See Data/Template/case_config.py "
              "for the pattern to add it.")

    outer_dist = None
    if hasattr(case_config, "mc_outer_distribution"):
        outer_dist = case_config.mc_outer_distribution()
    else:
        print(f"  NOTE: {case_name}/case_config.py has no mc_outer_distribution() -- "
              "mc_outer_peak_source/mc_outer_volume_source in scalars.csv must both "
              "be 'user' for this case (validate_mc_sources() will error otherwise, "
              "since 'curve'/'bootstrap'/'table' all need this hook). See "
              "Data/Template/case_config.py for the pattern to add it.")

    has_gate_hook = hasattr(case_config, "mc_gate_availability")
    mc_sources = validate_mc_sources(sc, outer_dist, case_name, has_gate_hook=has_gate_hook)

    gate_availability = None
    if has_gate_hook:
        sig = inspect.signature(case_config.mc_gate_availability)
        if "source" not in sig.parameters:
            raise MCConfigError(
                f"{case_name}/case_config.py's mc_gate_availability() still uses the "
                f"old USE_CI_RELIABILITY/GATE_P_FAIL env-var pattern -- that switch has "
                f"moved to scalars.csv's mc_gate_reliability_source/mc_gate_p_fail. "
                f"Update its signature to mc_gate_availability(source, p_fail=None), "
                f"selecting Tier 0 (flat) vs Tier 1 (CI-evidence) on the 'source' "
                f"argument instead of reading an environment variable. See "
                f"Data/Template/case_config.py for the current pattern.")
        gate_availability = case_config.mc_gate_availability(
            mc_sources["gate_source"], sc.get("mc_gate_p_fail"))
    else:
        print(f"  NOTE: {case_name}/case_config.py has no mc_gate_availability() -- "
              "inner loop will assume every gate stays operational (no gate "
              "failure-to-open risk modeled). See Data/Corumana_117_Q5000/"
              "case_config.py or Data/Template/case_config.py for the pattern to add it.")

    return (case_config, case_dir, reservoir_base, inflow_base, withdrawal, sc,
            rule_overrides, uncertain_params, outer_dist, gate_availability, mc_sources)


# ---------------------------------------------------------------------
# Main Monte Carlo driver
# ---------------------------------------------------------------------

def run_case_layer3(case_name: str, n_outer: int, n_inner: int, seed: int = 1,
                     rule: str = "baseline") -> tuple[list[dict], dict]:
    (case_config, case_dir, reservoir_base, inflow_base, withdrawal, sc,
     rule_overrides, uncertain_params, outer_dist,
     gate_availability, mc_sources) = build_case(case_name, rule)

    h0_sigma = sc.get("mc_h0_sigma", H0_DEFAULT_SIGMA)
    area_cv = sc.get("mc_area_cv", AREA_DEFAULT_CV)
    threshold = sc.get("downstream_threshold_m3s")
    max_flood_level = sc.get("max_flood_level")
    dam_crest_level = sc.get("dam_crest_level")
    H0_base = sc["H0"]
    t_max, dt = sc["t_max"], sc["dt"]
    print_every = int(sc.get("print_every", 1))

    # Every case explicitly declares, via scalars.csv, which source
    # drives the peak loop, the volume loop, and gate reliability --
    # validate_mc_sources() (called inside build_case()) already
    # enforced that the declared source and the data/scalars it needs
    # are both present; a case that's missing either fails BEFORE
    # reaching this point, not partway through a run. What follows is
    # a straight branch on that explicit choice -- no inference, no
    # "does this file happen to exist" fallback chain.
    explicit_peak_cv = sc.get("mc_outer_peak_cv")
    explicit_volume_cv = sc.get("mc_outer_volume_cv")
    peak_source = mc_sources["peak_source"]
    volume_source = mc_sources["volume_source"]
    target_return_period = outer_dist.get("target_return_period") if outer_dist else None

    # inflow_hydrograph.csv IS the design hydrograph -- its own peak and
    # volume are always the median (median_scale = volume_median_scale
    # = 1.0), for every source EXCEPT "bootstrap". curve_csv/
    # alt_studies_csv and volume_duration_csv are CV-ONLY inputs here:
    # they never recompute the anchor away from the case's own
    # hydrograph, because that hydrograph already IS the T-year design
    # flood, not an approximation of one that external files should
    # correct. "bootstrap" is the deliberate exception -- an FFA
    # stress-test case's whole point is that the FFA tool's own median
    # is a genuinely different, disputed central estimate (see
    # FFA_Integration.md), not a refinement of the hydrograph's.
    if peak_source == "user":
        median_scale = 1.0
        outer_cv = explicit_peak_cv
        outer_source = ("user (mc_outer_peak_source='user' -- no curve data used; "
                         f"cv={outer_cv:.4f} is a pure case-level judgment call)")

    elif peak_source == "curve":
        curve_T, curve_Q = load_flood_frequency_curve(outer_dist["curve_csv"])
        derived_cv = derive_outer_cv_from_alt_studies(
            curve_T, curve_Q, outer_dist["alt_studies_csv"])
        if derived_cv is None:
            raise MCConfigError(
                f"{case_name}: mc_outer_peak_source='curve' but "
                f"alt_studies_csv has fewer than 2 rows -- not enough to "
                f"derive a CV from. Add more rows, or switch to "
                f"mc_outer_peak_source='user' with an explicit mc_outer_peak_cv.")
        median_scale = 1.0
        if explicit_peak_cv is not None:
            outer_cv = explicit_peak_cv
            outer_source = ("curve (median: this case's own inflow_hydrograph.csv, "
                             "unchanged) + explicit mc_outer_peak_cv override "
                             f"(spread: case-level judgment, NOT the alt_studies-"
                             f"derived {derived_cv:.4f})")
        else:
            outer_cv = derived_cv
            outer_source = ("curve (median: this case's own inflow_hydrograph.csv, "
                             "unchanged; spread: alt_studies, real data)")
        print(f"  Outer loop (peak): mc_outer_peak_source='curve' -- median stays "
              f"at this case's own hydrograph peak (curve_csv/alt_studies_csv used "
              f"for CV only), cv={outer_cv:.4f} ({outer_source})")

    elif peak_source == "bootstrap":
        ci_T, ci_lower, ci_median, ci_upper = load_bootstrap_ci_curve(outer_dist["bootstrap_ci_csv"])
        ci_confidence = outer_dist.get("bootstrap_ci_confidence", 0.95)
        anchor_Q, derived_cv = derive_outer_cv_from_bootstrap_ci(
            ci_T, ci_lower, ci_median, ci_upper, target_return_period, ci_confidence)
        base_peak = float(max(inflow_base.Q))
        median_scale = anchor_Q / base_peak
        if explicit_peak_cv is not None:
            outer_cv = explicit_peak_cv
            outer_source = ("bootstrap (median: FFA bootstrap_ci_csv, DIFFERENT from "
                             "this case's own hydrograph peak by design -- the whole "
                             "point of this stress-test case) + explicit "
                             f"mc_outer_peak_cv override (spread: case-level judgment, "
                             f"NOT the bootstrap-derived {derived_cv:.4f})")
        else:
            outer_cv = derived_cv
            outer_source = (f"bootstrap (median AND spread both from FFA bootstrap_ci_csv "
                             f"at {ci_confidence:.0%} CI -- see FFA_Integration.md; median "
                             f"DELIBERATELY differs from this case's own hydrograph peak, "
                             f"this is almost always a stress-test case, not the primary case)")
        print(f"  Outer loop (peak): mc_outer_peak_source='bootstrap', targeting "
              f"T={target_return_period}-yr flood, bootstrap median={anchor_Q:.0f} m3/s vs. "
              f"this case's own hydrograph peak={base_peak:.0f} m3/s -> "
              f"median scale={median_scale:.4f}, cv={outer_cv:.4f} ({outer_source})")

    if volume_source == "user":
        volume_median_scale = 1.0
        volume_cv = explicit_volume_cv
        volume_source_desc = ("user (mc_outer_volume_source='user' -- no duration/volume data "
                               f"used; cv={volume_cv:.4f} is a pure case-level judgment call)")

    elif volume_source == "table":
        dv_T, dv_Q, dv_V = load_duration_volume_table(outer_dist["volume_duration_csv"])
        derived_volume_cv = derive_volume_cv_from_duration_table(dv_T, dv_Q, dv_V)
        if derived_volume_cv is None:
            raise MCConfigError(
                f"{case_name}: mc_outer_volume_source='table' but "
                f"volume_duration_csv has fewer than 2 rows -- not enough to "
                f"derive a CV from. Add more rows, or switch to "
                f"mc_outer_volume_source='user' with an explicit mc_outer_volume_cv.")
        volume_median_scale = 1.0
        if explicit_volume_cv is not None:
            volume_cv = explicit_volume_cv
            volume_source_desc = ("table (median: this case's own inflow_hydrograph.csv, "
                                   "unchanged) + explicit mc_outer_volume_cv override "
                                   f"(spread: case-level judgment, NOT the table-derived "
                                   f"{derived_volume_cv:.4f})")
        else:
            volume_cv = derived_volume_cv
            volume_source_desc = ("table (median: this case's own inflow_hydrograph.csv, "
                                   "unchanged; spread: flood_duration_volume_table, real data)")
        print(f"  Outer loop (volume, INDEPENDENT of peak): mc_outer_volume_source='table' "
              f"-- median stays at this case's own hydrograph volume (volume_duration_csv "
              f"used for CV only), cv={volume_cv:.4f} ({volume_source_desc})")

    volume_source = volume_source_desc

    dependence_source = mc_sources["dependence_source"]
    peak_volume_tau = sc.get("mc_peak_volume_tau")
    if dependence_source != "independent":
        print(f"  Outer loop (peak-volume DEPENDENCE): mc_peak_volume_dependence="
              f"'{dependence_source}', tau={peak_volume_tau:.3f} -- peak_scale and "
              f"volume_scale are correlated draws, not independent (see chat/literature "
              f"review on copula-based joint peak-volume sampling; marginal medians/CVs "
              f"above are UNCHANGED by this, only the correlation between the two draws).")

    rng = np.random.default_rng(seed)
    outer_draws = default_outer_sampler(rng, n_outer, outer_cv, volume_cv,
                                         median_scale, volume_median_scale,
                                         dependence=dependence_source, tau=peak_volume_tau)
    inner_draws_by_outer = [
        default_inner_sampler(rng, n_inner, uncertain_params, h0_sigma, area_cv,
                               gate_availability)
        for _ in range(n_outer)
    ]

    try:
        from tqdm import tqdm
        total = n_outer * n_inner
        pbar = tqdm(total=total, desc=f"Layer 3 MC ({case_name})")
    except ImportError:
        pbar = None
        print("(tqdm not installed -- no progress bar; `pip install tqdm` for one)")

    records = []
    for oi, outer in enumerate(outer_draws):
        scaled_inflow = ScaledHydrograph(inflow_base, outer["peak_scale"], outer["volume_scale"])
        for ii, inner in enumerate(inner_draws_by_outer[oi]):
            scaled_reservoir = ScaledReservoir(reservoir_base, inner["area_scale"])
            H0_draw = min(max(H0_base + inner["H0_shift"], reservoir_base.H_min),
                          reservoir_base.H_max)

            physical_overrides = {k: v for k, v in inner.items()
                                   if k not in ("H0_shift", "area_scale", "gates_failed", "ccf_fired")}
            outlets = call_build_outlets(case_config, rule_overrides,
                                          physical_overrides, scaled_inflow,
                                          gates_out_of_service=inner["gates_failed"])

            result = run_simulation(scaled_reservoir, scaled_inflow, withdrawal,
                                     outlets, H0=H0_draw, t_max=t_max, dt=dt,
                                     print_every=print_every)

            peak_level = float(result.H.max())
            release = downstream_release(result)
            peak_release = float(release.max())
            t_exceed = (time_to_exceedance(result.t_h, release, threshold)
                        if threshold is not None else None)

            record = {
                "outer_idx": oi, "inner_idx": ii,
                "peak_scale": outer["peak_scale"],
                "volume_scale": outer["volume_scale"],
                **physical_overrides,
                "H0_used": H0_draw,
                "area_scale": inner["area_scale"],
                "gates_failed": ",".join(inner["gates_failed"]),
                "n_gates_failed": len(inner["gates_failed"]),
                "ccf_fired": ";".join(inner["ccf_fired"]),
                "peak_level": peak_level,
                "peak_downstream_release": peak_release,
            }
            if t_exceed is not None:
                record["time_to_exceedance_h"] = t_exceed
            if max_flood_level is not None:
                record["exceeds_max_flood_level"] = int(peak_level > max_flood_level)
            if dam_crest_level is not None:
                record["exceeds_dam_crest"] = int(peak_level > dam_crest_level)
            records.append(record)

            if pbar is not None:
                pbar.update(1)

    if pbar is not None:
        pbar.close()

    meta = {
        "case_dir": case_dir, "threshold": threshold,
        "max_flood_level": max_flood_level, "dam_crest_level": dam_crest_level,
        "n_outer": n_outer, "n_inner": n_inner, "outer_cv": outer_cv,
        "outer_source": outer_source, "target_return_period": target_return_period,
        "median_scale": median_scale,
        "volume_cv": volume_cv, "volume_source": volume_source,
        "volume_median_scale": volume_median_scale,
        "h0_sigma": h0_sigma, "area_cv": area_cv, "rule": rule,
        "uncertain_params": uncertain_params,
        "gate_availability": gate_availability,
        "gate_source": mc_sources["gate_source"],
        "seed": seed,
    }
    return records, meta


# ---------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------

def _min_realizations_mean(sigma: float, mu: float, epsilon: float = 0.01, alpha: float = 0.95) -> float:
    """RMC-TotalRisk Technical Reference Manual, Equation 176: minimum
    Monte Carlo realizations for the OUTPUT MEAN to be precise to
    within relative error `epsilon` (e.g. 0.01 = 1%) at confidence
    `alpha`. sigma/mu are the output's own standard deviation/mean --
    NOT case inputs. Verified against the manual's own worked example
    (sigma=15, mu=100, epsilon=0.01, alpha=0.95 -> ~865).

    IMPORTANT for this project's nested two-loop design: the manual's
    formula assumes the realizations being counted are INDEPENDENT.
    Layer 3's n_outer*n_inner draws are not -- see
    _cluster_bootstrap_stats()'s docstring. Callers should compare this
    function's result against n_outer (the true count of independent
    flood scenarios), not n_outer*n_inner, or the resulting margin will
    be overstated."""
    z = NormalDist().inv_cdf((1 + alpha) / 2)
    if mu == 0:
        return float("inf")
    return ((sigma / (epsilon * abs(mu))) * z) ** 2


def _min_realizations_percentile(p: float, dp: float = 0.01, alpha: float = 0.95) -> float:
    """RMC-TotalRisk Technical Reference Manual, Equation 177: minimum
    Monte Carlo realizations for the p-th percentile (e.g. p=0.95) to
    be estimated within tolerance `dp` (in probability units, e.g.
    0.01) at confidence `alpha`. Verified against the manual's own
    worked example (p=0.95, dp=0.01, alpha=0.95 -> ~1825).

    Same effective-N caveat as _min_realizations_mean() above -- compare
    against n_outer, not n_outer*n_inner, for this project's nested
    design."""
    z = NormalDist().inv_cdf((1 + alpha) / 2)
    return p * (1 - p) * (z / dp) ** 2


def _cluster_bootstrap_stats(outer_peak_levels: list[np.ndarray], percentiles: list[float],
                              thresholds: list[float], n_boot: int = 500, seed: int = 0
                              ) -> tuple[dict[float, float], float, dict[float, float]]:
    """Block/cluster bootstrap standard errors for the percentiles, the
    mean, AND any exceedance-probability thresholds, all computed from
    the SAME resampled replicates so the three are internally
    consistent with each other.

    WHY this replaced a flat (per-draw) bootstrap: Layer 3's 15,000
    draws are NOT 15,000 independent observations -- they're 100
    independent outer-loop flood scenarios, each contributing 150
    correlated inner-loop draws (correlated because they share the
    same flood magnitude/volume). A flat bootstrap that resamples
    individual draws from the pooled 15,000 treats every draw as its
    own independent piece of evidence, which understates the true
    standard error for anything sensitive to WHICH outer scenarios got
    drawn -- confirmed empirically: two identical runs differing only
    in --seed showed P95 differing by ~11.5x the flat bootstrap's
    stated SE, and the mean by ~8x (both were run into the noise floor
    a genuine i.i.d. sample of this size shouldn't produce; the
    exceedance probabilities and P50 happened to be far less sensitive
    to this in that comparison, which is WHY they still checked out --
    not evidence the clustering doesn't apply to them in principle).

    Mechanism: each bootstrap replicate resamples WHICH of the n_outer
    scenarios to include (with replacement), keeping each selected
    scenario's full block of n_inner draws together -- so a replicate
    either includes an entire flood scenario's worth of correlated
    draws, or none of them, respecting the actual nested design instead
    of pretending every draw is independent. Requires every outer
    scenario to have contributed the same n_inner (true for this
    project's Layer 3 design -- every outer draw gets the same n_inner).

    Returns (pctl_se, mean_se, threshold_se) -- threshold_se keyed by
    the threshold values passed in. A threshold with zero exceedances
    across every draw returns an SE of 0.0 here (every bootstrap
    replicate is also zero) -- callers should use the separate
    ~3/n_outer bound for that case (see _format_exceedance()), not this
    function's SE, which is trivially uninformative when the estimate
    itself is exactly zero."""
    rng = np.random.default_rng(seed)
    n_outer = len(outer_peak_levels)
    block_matrix = np.array(outer_peak_levels)  # shape (n_outer, n_inner)
    chosen = rng.integers(0, n_outer, size=(n_boot, n_outer))
    replicates = block_matrix[chosen].reshape(n_boot, -1)  # (n_boot, n_outer*n_inner)

    boot_pctl = np.percentile(replicates, percentiles, axis=1)  # (len(percentiles), n_boot)
    pctl_se = {p: float(np.std(boot_pctl[i], ddof=1)) for i, p in enumerate(percentiles)}

    mean_se = float(np.std(replicates.mean(axis=1), ddof=1))

    threshold_se = {}
    for t in thresholds:
        boot_frac = np.mean(replicates > t, axis=1)  # (n_boot,)
        threshold_se[t] = float(np.std(boot_frac, ddof=1))

    return pctl_se, mean_se, threshold_se



def _format_exceedance(frac: float, se: float, n_outer: int, n_total: int) -> str:
    """Formats an exceedance-probability estimate with its cluster-
    bootstrap standard error (see _cluster_bootstrap_stats() -- NOT the
    closed-form binomial sqrt(p(1-p)/N_total), which assumes N_total
    independent draws and understates the true SE for the same nested-
    design reason percentile/mean SEs did), or, if ZERO events were
    observed in this run, the ~3/n_outer upper-bound caveat instead of
    a bare '0.0000' -- a bare zero reads as 'impossible', which it
    isn't; a rare event can easily produce zero observed occurrences in
    a finite run without the true probability being zero. n_outer, not
    n_total, is the right denominator for that bound too: an outer
    scenario severe enough to produce ANY exceedance would typically
    produce many correlated exceedances among its own inner draws, so
    the number of genuinely independent 'trials' this run tested is
    much closer to n_outer than to n_total."""
    if frac == 0.0:
        bound = 3.0 / n_outer
        return (f"0.0000 (0 of {n_total} draws across {n_outer} independent outer "
                f"scenarios -- true probability bounded above by roughly "
                f"3/n_outer ~= {bound:.4f}, NOT zero)")
    return f"{frac:.4f} (SE~={se:.4f}, cluster bootstrap)"


def summarize_and_write(case_name: str, records: list[dict], meta: dict) -> None:
    peak_levels = np.array([r["peak_level"] for r in records])
    n_outer, n_inner = meta["n_outer"], meta["n_inner"]
    outer_peak_levels = [
        np.array([r["peak_level"] for r in records if r["outer_idx"] == oi])
        for oi in range(n_outer)
    ]

    mc_dir = os.path.join(OUTPUT_DIR, case_name, "MonteCarlo")
    csv_path = os.path.join(mc_dir, "mc_results.csv")
    write_mc_results_csv(records, csv_path)

    plot_dir = os.path.join(PLOT_DIR, case_name, "MonteCarlo")
    plot_path = os.path.join(plot_dir, "mc_distribution.png")
    write_mc_distribution_plot(
        peak_levels, plot_path,
        title=f"Layer 3 Monte Carlo - {case_name} ({meta['rule']} rule, "
              f"{n_outer}x{n_inner}={n_outer * n_inner} draws)",
        max_flood_level=meta["max_flood_level"], dam_crest_level=meta["dam_crest_level"],
        outer_peak_levels=outer_peak_levels,
    )

    p5, p50, p95 = np.percentile(peak_levels, [5, 50, 95])
    n_total = n_outer * n_inner
    mean_level = float(np.mean(peak_levels))
    std_level = float(np.std(peak_levels, ddof=1))
    thresholds = [t for t in (meta["max_flood_level"], meta["dam_crest_level"]) if t is not None]
    pctl_se, se_mean, threshold_se = _cluster_bootstrap_stats(
        outer_peak_levels, [5, 50, 95], thresholds)
    # Eq. 176/177 (RMC-TotalRisk manual) convergence checks -- see the
    # two helper functions' docstrings for the verified formulas and
    # the manual's own worked examples they were checked against.
    # Compared against n_outer, NOT n_total -- see
    # _cluster_bootstrap_stats()'s docstring for why n_outer is the
    # true effective independent sample size for this nested design.
    n_min_mean = _min_realizations_mean(std_level, mean_level)
    n_min_p95 = _min_realizations_percentile(0.95)
    summary_path = os.path.join(mc_dir, "mc_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Case: {case_name}\n")
        f.write(f"Rule tested: {meta['rule']}\n")
        f.write(f"Seed: {meta.get('seed')}\n")
        f.write(f"Draws: {n_outer} outer x {n_inner} inner = {n_outer * n_inner}\n")
        f.write(f"Outer-loop source (peak): {meta['outer_source']}\n")
        if meta["target_return_period"] is not None:
            f.write(f"Outer-loop target return period: {meta['target_return_period']} years "
                    f"(median hydrograph PEAK scale factor: {meta['median_scale']:.4f})\n")
        f.write(f"Outer-loop PEAK scale CV: {meta['outer_cv']:.4f}\n")
        f.write(f"Outer-loop source (volume, sampled INDEPENDENTLY of peak): "
                f"{meta['volume_source']}\n")
        f.write(f"Outer-loop median VOLUME scale factor: {meta['volume_median_scale']:.4f}, "
                f"CV: {meta['volume_cv']:.4f}\n")
        f.write(f"Inner-loop H0 sigma [m]: {meta['h0_sigma']}, "
                f"reservoir area_scale CV: {meta['area_cv']}\n")
        f.write(f"Inner-loop uncertain physical params: "
                f"{list(meta['uncertain_params'].keys()) or '(none declared)'}\n")
        ga = normalize_gate_availability(meta.get("gate_availability"))
        if ga:
            gates = ga["gates"]
            n = len(gates)
            rates = list(ga["p_fail_by_gate"].values())
            n_failed = np.array([r["n_gates_failed"] for r in records])
            frac_any_failed = float(np.mean(n_failed >= 1))
            expected_any_failed = expected_p_any_gate_failed(ga["p_fail_by_gate"], ga["ccf_groups"])
            if len(set(round(r, 6) for r in rates)) <= 1:
                rate_desc = f"p_fail={rates[0]:.4f} for every gate" if rates else "n/a"
            else:
                rate_desc = f"p_fail ranges {min(rates):.4f}-{max(rates):.4f} across gates (per-gate CI-derived rates)"
            f.write(f"Gate availability: {n} gates ({', '.join(gates)}), {rate_desc}, "
                    f"independent draws\n")
            if ga["ccf_groups"]:
                for grp in ga["ccf_groups"]:
                    fired = np.array([grp["label"] in r["ccf_fired"].split(";") for r in records])
                    f.write(f"  Common-cause branch '{grp['label']}' (affects {len(grp['gates'])} "
                            f"gates): declared p_ccf={grp['p_ccf']:.4f}, empirical fired fraction "
                            f"over this run: {float(np.mean(fired)):.4f}\n")
            else:
                f.write("  Common-cause failure: not modeled (no ccf_groups declared)\n")
            f.write(f"P(>=1 gate failed) -- empirical over this run: {frac_any_failed:.4f}, "
                    f"closed-form: {expected_any_failed:.4f}\n")
            if ga.get("component_pfs"):
                f.write("\nImportance ranking (which subsystem/common-cause branch is "
                        "carrying the most probability weight -- see Module/reliability/"
                        "importance.py's docstring for what this ranking does and does not "
                        "mean; it is a first-pass ranking, not a rigorous sensitivity measure):\n")
                ranked = rank_importance(ga["component_pfs"], ga["ccf_groups"])
                f.write(format_importance_table(ranked) + "\n")
            f.write(f"  Gate reliability source: mc_gate_reliability_source="
                    f"'{meta['gate_source']}'\n")
        else:
            f.write("Gate availability: not modeled (no mc_gate_availability() declared)\n")
        f.write("\n")
        f.write(f"Peak level [m a.s.l.]: P5={p5:.3f} (SE~={pctl_se[5]:.3f})  "
                f"P50={p50:.3f} (SE~={pctl_se[50]:.3f})  P95={p95:.3f} (SE~={pctl_se[95]:.3f})  "
                f"mean={mean_level:.3f} (SE~={se_mean:.3f})  "
                f"max={peak_levels.max():.3f}  min={peak_levels.min():.3f}\n")
        f.write(f"Convergence (RMC-TotalRisk manual Eq. 176/177, target: mean within +/-1%, "
                f"P95 within +/-0.01, both at 95% confidence). NOTE: margin below is computed "
                f"against n_outer={n_outer}, the effective INDEPENDENT sample size for this "
                f"nested two-loop design, NOT the raw {n_total} total draws -- the {n_inner} "
                f"inner draws per outer scenario are correlated (same flood magnitude/volume), "
                f"so they don't each count as fresh independent evidence for outer-loop-"
                f"sensitive statistics like the mean and P95. See "
                f"_cluster_bootstrap_stats()'s docstring for the empirical confirmation "
                f"(cross-seed comparison) that motivated this: "
                f"mean needs >={n_min_mean:.0f} independent realizations (this run: {n_outer}, "
                f"{n_outer / n_min_mean:.1f}x margin); "
                f"P95 needs >={n_min_p95:.0f} independent realizations (this run: {n_outer}, "
                f"{n_outer / n_min_p95:.1f}x margin)\n")
        if meta["max_flood_level"] is not None:
            frac = float(np.mean(peak_levels > meta["max_flood_level"]))
            f.write(f"P(peak level > max_flood_level={meta['max_flood_level']}): "
                    f"{_format_exceedance(frac, threshold_se.get(meta['max_flood_level'], 0.0), n_outer, n_total)}\n")
        if meta["dam_crest_level"] is not None:
            frac = float(np.mean(peak_levels > meta["dam_crest_level"]))
            f.write(f"P(peak level > dam_crest_level={meta['dam_crest_level']}): "
                    f"{_format_exceedance(frac, threshold_se.get(meta['dam_crest_level'], 0.0), n_outer, n_total)}\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {plot_path}")
    print(f"Wrote {summary_path}")
    ga = normalize_gate_availability(meta.get("gate_availability"))
    if ga:
        n_failed = np.array([r["n_gates_failed"] for r in records])
        frac_any_failed = float(np.mean(n_failed >= 1))
        expected_any_failed = expected_p_any_gate_failed(ga["p_fail_by_gate"], ga["ccf_groups"])
        print(f"P(>=1 gate failed) = {frac_any_failed:.4f} "
              f"(closed-form: {expected_any_failed:.4f}, n={len(ga['gates'])}, "
              f"{len(ga['ccf_groups'])} CCF group(s))")
        if ga.get("component_pfs"):
            print("\nImportance ranking:")
            print(format_importance_table(rank_importance(ga["component_pfs"], ga["ccf_groups"])))
    print(f"\nPeak level [m a.s.l.]: P5={p5:.3f} (SE~={pctl_se[5]:.3f})  "
          f"P50={p50:.3f} (SE~={pctl_se[50]:.3f})  P95={p95:.3f} (SE~={pctl_se[95]:.3f})  "
          f"mean={mean_level:.3f} (SE~={se_mean:.3f}) "
          f"(max_flood_level={meta['max_flood_level']}, dam_crest_level={meta['dam_crest_level']})")
    print(f"Convergence (Eq. 176/177, vs. n_outer={n_outer} independent scenarios, "
          f"NOT n_total={n_total} -- see mc_summary.txt for why): "
          f"mean needs >={n_min_mean:.0f} realizations "
          f"({n_outer / n_min_mean:.1f}x margin); P95 needs >={n_min_p95:.0f} realizations "
          f"({n_outer / n_min_p95:.1f}x margin)")
    if meta["max_flood_level"] is not None:
        frac = float(np.mean(peak_levels > meta["max_flood_level"]))
        print(f"P(exceeds max_flood_level) = "
              f"{_format_exceedance(frac, threshold_se.get(meta['max_flood_level'], 0.0), n_outer, n_total)}")
    if meta["dam_crest_level"] is not None:
        frac = float(np.mean(peak_levels > meta["dam_crest_level"]))
        print(f"P(exceeds dam_crest_level) = "
              f"{_format_exceedance(frac, threshold_se.get(meta['dam_crest_level'], 0.0), n_outer, n_total)}")


# ---------------------------------------------------------------------
# CLI (mirrors run_case.py / optimize_layer2.py conventions)
# ---------------------------------------------------------------------

def list_case_folders():
    cases = []
    if not os.path.isdir(DATA_DIR):
        return cases
    for name in sorted(os.listdir(DATA_DIR)):
        case_dir = os.path.join(DATA_DIR, name)
        if os.path.isdir(case_dir) and os.path.isfile(os.path.join(case_dir, "case_config.py")):
            cases.append(name)
    return cases


def prompt_for_case(cases):
    if not cases:
        raise SystemExit(f"No case folders found under {DATA_DIR}.")
    print("Available input cases (Data/<name>/):")
    for i, name in enumerate(cases, start=1):
        print(f"  {i}. {name}")
    while True:
        choice = input(f"Select a case [1-{len(cases)}] or type its name: ").strip()
        if choice in cases:
            return choice
        if choice.isdigit() and 1 <= int(choice) <= len(cases):
            return cases[int(choice) - 1]
        print("Not recognized, try again.")


def prompt_int(label: str, default: int) -> int:
    raw = input(f"{label} [default {default}]: ").strip()
    if raw == "":
        return default
    try:
        return int(raw)
    except ValueError:
        print(f"Not a valid integer, using default ({default}).")
        return default


def main():
    parser = argparse.ArgumentParser(description="Layer 3: Monte Carlo / climate-uncertainty analysis")
    parser.add_argument("case_name", nargs="?", default=None)
    parser.add_argument("--n-outer", type=int, default=None)
    parser.add_argument("--n-inner", type=int, default=None)
    parser.add_argument("--rule", choices=["baseline", "optimized"], default="baseline")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    case_name = args.case_name or prompt_for_case(list_case_folders())
    n_outer = args.n_outer if args.n_outer is not None else prompt_int("Number of outer (climate) draws", 100)
    n_inner = args.n_inner if args.n_inner is not None else prompt_int("Number of inner (physical) draws", 150)

    print(f"\nRunning Layer 3 Monte Carlo for '{case_name}': "
          f"{n_outer} outer x {n_inner} inner = {n_outer * n_inner} simulations, "
          f"rule='{args.rule}'\n")

    records, meta = run_case_layer3(case_name, n_outer, n_inner, seed=args.seed, rule=args.rule)
    summarize_and_write(case_name, records, meta)


if __name__ == "__main__":
    main()
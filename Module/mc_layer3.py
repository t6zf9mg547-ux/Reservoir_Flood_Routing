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
    horizon. Sampled log-normally, median = 1.0 (the case's own
    hydrograph IS the central estimate), spread set by
    `mc_outer_peak_cv` in scalars.csv (default OUTER_DEFAULT_CV below).

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

# Reuse Layer 2's case loader + downstream-release/exceedance-time
# definitions rather than re-implementing them -- same conventions,
# one source of truth. optimize_layer2.py only imports pymoo INSIDE
# its own functions, so importing these three names doesn't require
# pymoo to be installed.
from optimize_layer2 import load_case_config, downstream_release, time_to_exceedance

OUTER_DEFAULT_CV = 0.15   # default spread on the outer-loop hydrograph
                          # peak scale factor if scalars.csv doesn't
                          # set mc_outer_peak_cv
VOLUME_DEFAULT_CV = 0.15  # default spread on the outer-loop hydrograph
                          # VOLUME scale factor (independent of peak --
                          # see ScaledHydrograph) if no real duration/
                          # volume data is available for a case
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
    format) -- absent that, a generic placeholder CV is used instead.

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
    """n draws from a log-normal distribution with median 1.0 and the
    given coefficient of variation (matches how a "+/- X% uncertainty
    multiplier, centered on the nominal value" is usually meant)."""
    if cv <= 0:
        return np.ones(n)
    sigma = np.sqrt(np.log(1.0 + cv ** 2))
    mu = -0.5 * sigma ** 2  # so that median = exp(mu) = 1.0
    return rng.lognormal(mean=mu, sigma=sigma, size=n)


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


def default_outer_sampler(rng: np.random.Generator, n_outer: int,
                           peak_cv: float, volume_cv: float,
                           peak_median_scale: float = 1.0,
                           volume_median_scale: float = 1.0) -> list[dict]:
    """OUTER loop draws: peak_scale and volume_scale, sampled fully
    INDEPENDENTLY (per-user decision -- floods with the same peak can
    plausibly have quite different volumes, or vice versa, and forcing
    them to move together would hide exactly the routing-sensitive
    scenarios this dam cares about most; see ScaledHydrograph's
    docstring). Both log-normal. Swap this out for a real fitted-
    distribution sampler with actual parameter-covariance/bootstrap
    uncertainty (and, ideally, a real joint peak-volume distribution)
    once that data exists -- nothing downstream depends on how these
    scales were generated."""
    peak_scales = _lognormal_cv(rng, peak_cv, n_outer) * peak_median_scale
    volume_scales = _lognormal_cv(rng, volume_cv, n_outer) * volume_median_scale
    return [{"peak_scale": float(p), "volume_scale": float(v)}
            for p, v in zip(peak_scales, volume_scales)]


def default_inner_sampler(rng: np.random.Generator, n_inner: int,
                           uncertain_params: dict, h0_sigma: float,
                           area_cv: float, gate_availability: dict | None = None
                           ) -> list[dict]:
    """INNER loop draws: physical/rating multipliers (whatever the case
    declares via mc_uncertain_params()), plus the two generic,
    case-independent draws every case gets (H0 shift and reservoir
    area_scale), plus -- if the case declares mc_gate_availability() --
    which (if any) of its named gates fail to open this draw.

    Gate failure: every gate in gate_availability["gates"] is drawn as
    an INDEPENDENT Bernoulli(p_fail) trial, using the SAME p_fail for
    every gate (per-user decision -- gate unavailability is modeled as
    a shared equipment/maintenance-driven rate across identical gates,
    not per-gate-specific data). This exactly reproduces the standard
    binomial screening model P(K=k) = C(n,k)*p^k*(1-p)^(n-k) for the
    number of gates K that fail to open, with n = len(gates). Common-
    cause (correlated) failure is explicitly NOT modeled here -- these
    draws are independent, matching the binomial model's own stated
    scope (see the reference table this was validated against)."""
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

    gates = list(gate_availability["gates"]) if gate_availability else []
    p_fail = gate_availability.get("p_fail", 0.0) if gate_availability else 0.0
    # n_inner x n_gates independent Bernoulli(p_fail) trials, drawn
    # once here (vectorized, same as every other inner-loop quantity)
    # rather than per-draw in the caller's loop.
    gate_fail_draws = (rng.random((n_inner, len(gates))) < p_fail) if gates else None

    for i in range(n_inner):
        d = {name: float(col[i]) for name, col in sampled_cols.items()}
        d["H0_shift"] = float(h0_shift[i])
        d["area_scale"] = float(area_scale[i])
        d["gates_failed"] = ([g for g, failed in zip(gates, gate_fail_draws[i]) if failed]
                              if gates else [])
        draws.append(d)
    return draws


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
              "outer loop will use the generic placeholder (a flat CV around this "
              "case's own hydrograph peak), not a real flood-frequency curve. See "
              "Data/Template/case_config.py for the pattern to add it.")

    gate_availability = None
    if hasattr(case_config, "mc_gate_availability"):
        gate_availability = case_config.mc_gate_availability()
    else:
        print(f"  NOTE: {case_name}/case_config.py has no mc_gate_availability() -- "
              "inner loop will assume every gate stays operational (no gate "
              "failure-to-open risk modeled). See Data/Corumana_117_Q5000/"
              "case_config.py or Data/Template/case_config.py for the pattern to add it.")

    return (case_config, case_dir, reservoir_base, inflow_base, withdrawal, sc,
            rule_overrides, uncertain_params, outer_dist, gate_availability)


# ---------------------------------------------------------------------
# Main Monte Carlo driver
# ---------------------------------------------------------------------

def run_case_layer3(case_name: str, n_outer: int, n_inner: int, seed: int = 1,
                     rule: str = "baseline") -> tuple[list[dict], dict]:
    (case_config, case_dir, reservoir_base, inflow_base, withdrawal, sc,
     rule_overrides, uncertain_params, outer_dist,
     gate_availability) = build_case(case_name, rule)

    h0_sigma = sc.get("mc_h0_sigma", H0_DEFAULT_SIGMA)
    area_cv = sc.get("mc_area_cv", AREA_DEFAULT_CV)
    threshold = sc.get("downstream_threshold_m3s")
    max_flood_level = sc.get("max_flood_level")
    dam_crest_level = sc.get("dam_crest_level")
    H0_base = sc["H0"]
    t_max, dt = sc["t_max"], sc["dt"]
    print_every = int(sc.get("print_every", 1))

    # A case can set mc_outer_peak_cv/mc_outer_volume_cv in its OWN
    # scalars.csv as an EXPLICIT judgment call on the spread, distinct
    # from just leaving it unset (which falls back to either the
    # generic OUTER_DEFAULT_CV/VOLUME_DEFAULT_CV, or -- if real
    # alt_studies/duration data is available below -- the data-derived
    # spread). Captured as None-if-absent here (not yet defaulted) so
    # the derivation logic below can tell "the case didn't set this"
    # apart from "the case explicitly wants this exact number," and
    # let an explicit override win even when real derivation data is
    # ALSO available -- i.e. "use the real curve/table for the MEDIAN
    # (the anchor value), but I want to set the SPREAD myself" is now
    # a supported combination, not just all-or-nothing.
    explicit_peak_cv = sc.get("mc_outer_peak_cv")
    explicit_volume_cv = sc.get("mc_outer_volume_cv")

    outer_source = "placeholder"
    outer_cv = explicit_peak_cv if explicit_peak_cv is not None else OUTER_DEFAULT_CV
    volume_cv = explicit_volume_cv if explicit_volume_cv is not None else VOLUME_DEFAULT_CV
    median_scale = 1.0
    volume_median_scale = 1.0
    target_return_period = None
    if outer_dist is not None:
        try:
            curve_T, curve_Q = load_flood_frequency_curve(outer_dist["curve_csv"])
            target_return_period = outer_dist["target_return_period"]
            anchor_Q = _loglog_interp(target_return_period, curve_T, curve_Q)
            base_peak = float(max(inflow_base.Q))
            median_scale = anchor_Q / base_peak
            derived_cv = derive_outer_cv_from_alt_studies(
                curve_T, curve_Q, outer_dist["alt_studies_csv"])
            if explicit_peak_cv is not None:
                outer_cv = explicit_peak_cv
                outer_source = ("flood_frequency_curve (median: real data) + explicit "
                                 "mc_outer_peak_cv override (spread: case-level judgment"
                                 + (f", NOT the alt_studies-derived {derived_cv:.4f}"
                                    if derived_cv is not None else "") + ")")
            elif derived_cv is not None:
                outer_cv = derived_cv
                outer_source = "flood_frequency_curve + alt_studies (real data)"
            else:
                outer_source = "flood_frequency_curve (median only; CV still placeholder)"
            print(f"  Outer loop (peak): targeting T={target_return_period}-yr flood, "
                  f"curve peak={anchor_Q:.0f} m3/s vs. this case's own hydrograph peak="
                  f"{base_peak:.0f} m3/s -> median scale={median_scale:.4f}, "
                  f"cv={outer_cv:.4f} ({outer_source})")
        except (FileNotFoundError, KeyError, ValueError, StopIteration) as e:
            target_return_period = None
            median_scale = 1.0
            outer_cv = explicit_peak_cv if explicit_peak_cv is not None else OUTER_DEFAULT_CV
            outer_source = "placeholder (mc_outer_distribution() set, but failed to load -- see warning above)"
            print(f"  WARNING: {case_name}/case_config.py's mc_outer_distribution() points at "
                  f"peak data that couldn't be read ({type(e).__name__}: {e}). Falling back to "
                  f"the generic placeholder for the PEAK loop instead of stopping the run. Check "
                  f"the 'curve_csv'/'alt_studies_csv' paths returned by mc_outer_distribution() "
                  f"actually exist and are formatted as described in Data/Template/case_config.py.")

    volume_source = "placeholder (locked to peak_scale)"
    if outer_dist is not None and "volume_duration_csv" in outer_dist:
        try:
            dv_T, dv_Q, dv_V = load_duration_volume_table(outer_dist["volume_duration_csv"])
            anchor_T = target_return_period if target_return_period is not None \
                else outer_dist["target_return_period"]
            anchor_V = _loglog_interp(anchor_T, dv_T, dv_V)
            base_volume = float(np.sum(np.diff(inflow_base.t) *
                                        (inflow_base.Q[:-1] + inflow_base.Q[1:]) / 2.0))  # m3
            volume_median_scale = anchor_V / base_volume
            derived_volume_cv = derive_volume_cv_from_duration_table(dv_T, dv_Q, dv_V)
            if explicit_volume_cv is not None:
                volume_cv = explicit_volume_cv
                volume_source = ("flood_duration_volume_table (median: real data) + explicit "
                                  "mc_outer_volume_cv override (spread: case-level judgment"
                                  + (f", NOT the table-derived {derived_volume_cv:.4f}"
                                     if derived_volume_cv is not None else "") + ")")
            elif derived_volume_cv is not None:
                volume_cv = derived_volume_cv
                volume_source = "flood_duration_volume_table (real data, independent of peak)"
            print(f"  Outer loop (volume, INDEPENDENT of peak): targeting T="
                  f"{anchor_T}-yr flood, table volume={anchor_V/1e6:.0f} Mm3 vs. "
                  f"this case's own hydrograph volume={base_volume/1e6:.0f} Mm3 -> "
                  f"median scale={volume_median_scale:.4f}, cv={volume_cv:.4f} ({volume_source})")
        except (FileNotFoundError, KeyError, ValueError, StopIteration) as e:
            volume_median_scale = 1.0
            volume_cv = explicit_volume_cv if explicit_volume_cv is not None else VOLUME_DEFAULT_CV
            volume_source = "placeholder (mc_outer_distribution() set, but failed to load -- see warning above)"
            print(f"  WARNING: {case_name}/case_config.py's mc_outer_distribution() points at "
                  f"volume data that couldn't be read ({type(e).__name__}: {e}). Falling back to "
                  f"the generic placeholder for the VOLUME loop instead of stopping the run. Check "
                  f"the 'volume_duration_csv' path returned by mc_outer_distribution() actually "
                  f"exists and is formatted as described in Data/Template/case_config.py.")

    rng = np.random.default_rng(seed)
    outer_draws = default_outer_sampler(rng, n_outer, outer_cv, volume_cv,
                                         median_scale, volume_median_scale)
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
                                   if k not in ("H0_shift", "area_scale", "gates_failed")}
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
    }
    return records, meta


# ---------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------

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
    summary_path = os.path.join(mc_dir, "mc_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Case: {case_name}\n")
        f.write(f"Rule tested: {meta['rule']}\n")
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
        ga = meta.get("gate_availability")
        if ga:
            gates, p_fail = list(ga["gates"]), ga.get("p_fail", 0.0)
            n = len(gates)
            n_failed = np.array([r["n_gates_failed"] for r in records])
            frac_any_failed = float(np.mean(n_failed >= 1))
            # Closed-form P(K>=1) for n independent Bernoulli(p_fail)
            # trials -- the same binomial screening model the gate-
            # availability sampling itself implements -- as a sanity
            # check that the empirical draws actually reproduce it.
            expected_any_failed = 1.0 - (1.0 - p_fail) ** n if n else 0.0
            f.write(f"Gate availability: {n} gates ({', '.join(gates)}), per-gate "
                    f"failure-to-open probability p_fail={p_fail:.4f} (independent "
                    f"draws, SAME rate for every gate; common-cause failure NOT modeled)\n")
            f.write(f"P(>=1 gate failed) -- empirical over this run: {frac_any_failed:.4f}, "
                    f"closed-form binomial: {expected_any_failed:.4f}\n")
        else:
            f.write("Gate availability: not modeled (no mc_gate_availability() declared)\n")
        f.write("\n")
        f.write(f"Peak level [m a.s.l.]: P5={p5:.3f}  P50={p50:.3f}  P95={p95:.3f}  "
                f"max={peak_levels.max():.3f}  min={peak_levels.min():.3f}\n")
        if meta["max_flood_level"] is not None:
            frac = float(np.mean(peak_levels > meta["max_flood_level"]))
            f.write(f"P(peak level > max_flood_level={meta['max_flood_level']}): {frac:.4f}\n")
        if meta["dam_crest_level"] is not None:
            frac = float(np.mean(peak_levels > meta["dam_crest_level"]))
            f.write(f"P(peak level > dam_crest_level={meta['dam_crest_level']}): {frac:.4f}\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {plot_path}")
    print(f"Wrote {summary_path}")
    ga = meta.get("gate_availability")
    if ga:
        n_failed = np.array([r["n_gates_failed"] for r in records])
        n_gates = len(ga["gates"])
        p_fail = ga.get("p_fail", 0.0)
        frac_any_failed = float(np.mean(n_failed >= 1))
        expected_any_failed = 1.0 - (1.0 - p_fail) ** n_gates if n_gates else 0.0
        print(f"P(>=1 gate failed) = {frac_any_failed:.4f} "
              f"(closed-form binomial: {expected_any_failed:.4f}, "
              f"n={n_gates}, p_fail={p_fail:.4f})")
    print(f"\nPeak level [m a.s.l.]: P5={p5:.3f}  P50={p50:.3f}  P95={p95:.3f} "
          f"(max_flood_level={meta['max_flood_level']}, dam_crest_level={meta['dam_crest_level']})")
    if meta["max_flood_level"] is not None:
        frac = float(np.mean(peak_levels > meta["max_flood_level"]))
        print(f"P(exceeds max_flood_level) = {frac:.4f}")
    if meta["dam_crest_level"] is not None:
        frac = float(np.mean(peak_levels > meta["dam_crest_level"]))
        print(f"P(exceeds dam_crest_level) = {frac:.4f}")


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
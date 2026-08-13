"""
mc_layer3_outer_only.py -- outer-loop-ONLY Layer 3 sensitivity run
--------------------------------------------------------------------------
Isolates the OUTER loop's contribution to peak-level spread (climate/
flood-frequency uncertainty -- peak_scale/volume_scale, or a hydrograph
ensemble) from the INNER loop's (physical/rating-coefficient and gate-
reliability uncertainty), by running each outer scenario exactly ONCE,
with every inner-loop quantity held at its deterministic, unperturbed
BASELINE value:

  - H0            : exactly this case's own scalars.csv H0 (no h0_sigma
                    shift)
  - area_scale    : 1.0 (no stage-storage curve perturbation)
  - every physical/rating multiplier declared in mc_uncertain_params():
                    1.0 for "lognormal_cv" params (nominal, no
                    perturbation), or the declared "mean" for "normal"
                    params (nominal, no shift) -- see
                    baseline_inner_draw()'s docstring
  - gates         : ALL operational, no failures at all (deterministic
                    best case) -- this is a deliberate simplification,
                    not a claim that gate failure is impossible; see
                    the module's own chat history for the reasoning

This answers a genuinely different question from the full Layer 3 run:
"how much does peak level vary due to WHICH FLOOD we get, holding every
physical/reliability uncertainty at its nominal value?" -- as opposed
to full Layer 3's "how much does peak level vary from everything at
once?". Comparing the two run's spreads is a direct, empirical way to
see which loop dominates the reported uncertainty for a given case,
rather than inferring it indirectly.

Reuses Module/mc_layer3.py's build_case()/build_outer_draws() (the
SAME outer-loop derivation logic -- curve/bootstrap/user/table sources,
copula dependence, hydrograph ensembles -- as the full Layer 3 run;
this script deliberately does NOT duplicate that logic) and its
_cluster_bootstrap_stats()/_min_realizations_*()/_format_margin()/
_format_exceedance() reporting machinery, which -- with exactly one
deterministic "inner draw" per outer scenario -- naturally degenerates
into a standard, textbook i.i.d. bootstrap/Eq.176/177 check over
n_outer independent scenarios (no nested-loop correction needed, since
there IS no inner-loop nesting here).

Writes to Output/<CaseName>/OuterOnly/ and Plot/<CaseName>/OuterOnly/
-- SEPARATE from full Layer 3's Output/<CaseName>/MonteCarlo/, so
running this never overwrites (or gets confused with) a full run's
results.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(MODULE_DIR)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "Output")
PLOT_DIR = os.path.join(PROJECT_ROOT, "Plot")

sys.path.insert(0, MODULE_DIR)

from solver import run_simulation
from report import write_mc_results_csv, write_mc_distribution_plot
from optimize_layer2 import downstream_release, time_to_exceedance

import mc_layer3
from mc_layer3 import (
    build_case, call_build_outlets, ScaledHydrograph, ScaledReservoir,
    _cluster_bootstrap_stats, _min_realizations_mean,
    _min_realizations_percentile, _format_margin, _format_exceedance,
)


def baseline_inner_draw(uncertain_params: dict) -> dict:
    """The single, deterministic 'inner draw' this script uses for
    EVERY outer scenario -- same dict shape default_inner_sampler()
    produces per-draw (see its own docstring), but with every value
    fixed at its nominal/unperturbed baseline instead of sampled:

    - each declared physical/rating multiplier: 1.0 for
      dist="lognormal_cv" (the median of a CV=X spread around nominal
      IS 1.0 -- no perturbation), or the declared "mean" for
      dist="normal" (a "normal" param is typically a SHIFT, e.g.
      fuse_trigger_shift in meters, not a multiplier -- its nominal/
      unperturbed value is whatever "mean" was declared, often but not
      necessarily 0.0)
    - H0_shift=0.0, area_scale=1.0 (no perturbation)
    - gates_failed=[], ccf_fired=[] (ALL gates operational -- this
      script's deliberate, documented simplification; see module
      docstring)
    """
    d = {}
    for name, spec in uncertain_params.items():
        dist = spec.get("dist", "lognormal_cv")
        if dist == "lognormal_cv":
            d[name] = 1.0
        elif dist == "normal":
            d[name] = float(spec.get("mean", 0.0))
        else:
            raise ValueError(f"mc_uncertain_params(): unknown dist '{dist}' for '{name}'")
    d["H0_shift"] = 0.0
    d["area_scale"] = 1.0
    d["gates_failed"] = []
    d["ccf_fired"] = []
    return d


def run_case_outer_only(case_name: str, n_outer: int, seed: int = 1,
                         rule: str = "baseline") -> tuple[list[dict], dict]:
    """Runs each outer scenario EXACTLY ONCE, inner loop held at its
    deterministic baseline (see baseline_inner_draw()). Reuses
    mc_layer3.build_case()/build_outer_draws() verbatim for the outer-
    loop derivation -- same peak/volume/dependence/ensemble sources, same
    validation, same n_outer-capping behavior in ensemble mode -- so
    this script's outer scenarios are identical to what a full Layer 3
    run with the same scalars.csv would have drawn."""
    (case_config, case_dir, reservoir_base, inflow_base, withdrawal, sc,
     rule_overrides, uncertain_params, outer_dist,
     gate_availability, mc_sources) = build_case(case_name, rule)

    threshold = sc.get("downstream_threshold_m3s")
    max_flood_level = sc.get("max_flood_level")
    dam_crest_level = sc.get("dam_crest_level")
    H0_base = sc["H0"]
    t_max, dt = sc["t_max"], sc["dt"]
    print_every = int(sc.get("print_every", 1))

    rng = np.random.default_rng(seed)

    outer = build_outer_draws_(case_config, case_dir, inflow_base, sc, outer_dist,
                                mc_sources, n_outer, rng, case_name)
    outer_draws = outer["outer_draws"]
    n_outer = outer["n_outer"]

    inner = baseline_inner_draw(uncertain_params)

    print(f"  Inner loop: BASELINE ONLY (this script) -- H0={H0_base}, area_scale=1.0, "
          f"every physical multiplier at its nominal value, ALL gates operational "
          f"(no failures modeled). {len(uncertain_params)} physical params held at "
          f"baseline: {list(uncertain_params.keys())}")

    try:
        from tqdm import tqdm
        pbar = tqdm(total=n_outer, desc=f"Layer 3 outer-only ({case_name})")
    except ImportError:
        pbar = None
        print("(tqdm not installed -- no progress bar; `pip install tqdm` for one)")

    records = []
    for oi, o in enumerate(outer_draws):
        scaled_inflow = o.get("hydrograph") or \
            ScaledHydrograph(inflow_base, o["peak_scale"], o["volume_scale"])

        H0_draw = min(max(H0_base + inner["H0_shift"], reservoir_base.H_min),
                      reservoir_base.H_max)
        physical_overrides = {k: v for k, v in inner.items()
                               if k not in ("H0_shift", "area_scale", "gates_failed", "ccf_fired")}
        outlets = call_build_outlets(case_config, rule_overrides,
                                      physical_overrides, scaled_inflow,
                                      gates_out_of_service=inner["gates_failed"])

        result = run_simulation(ScaledReservoir(reservoir_base, inner["area_scale"]),
                                 scaled_inflow, withdrawal,
                                 outlets, H0=H0_draw, t_max=t_max, dt=dt,
                                 print_every=print_every)

        peak_level = float(result.H.max())
        release = downstream_release(result)
        peak_release = float(release.max())
        t_exceed = (time_to_exceedance(result.t_h, release, threshold)
                    if threshold is not None else None)

        record = {
            "outer_idx": oi,
            "peak_scale": o["peak_scale"],
            "volume_scale": o["volume_scale"],
            **physical_overrides,
            "H0_used": H0_draw,
            "area_scale": inner["area_scale"],
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
        "n_outer": n_outer,
        "outer_source": outer["outer_source"], "target_return_period": outer["target_return_period"],
        "median_scale": outer["median_scale"],
        "volume_cv": outer["volume_cv"], "volume_source": outer["volume_source"],
        "volume_median_scale": outer["volume_median_scale"],
        "outer_cv": outer["outer_cv"],
        "dependence_source": outer["dependence_source"], "peak_volume_tau": outer["peak_volume_tau"],
        "rule": rule, "uncertain_params": uncertain_params, "seed": seed,
    }
    return records, meta


def build_outer_draws_(case_config, case_dir, inflow_base, sc, outer_dist,
                        mc_sources, n_outer, rng, case_name):
    """Thin wrapper -- mc_layer3.build_outer_draws() takes case_name
    first; build_case() here returns case_config/case_dir in a
    different position than that function's own parameter order, so
    this just calls it correctly rather than repeating the whole
    derivation."""
    return mc_layer3.build_outer_draws(case_name, case_dir, inflow_base, sc,
                                        outer_dist, mc_sources, n_outer, rng)


def summarize_and_write(case_name: str, records: list[dict], meta: dict) -> None:
    """Simplified relative to mc_layer3.summarize_and_write() -- no
    gate-reliability empirical-fired-fraction reporting (meaningless
    here: gates are ALWAYS operational by this script's own design, not
    something to report statistics about), no nested-loop SE
    correction needed (there is no inner-loop nesting -- n_outer IS
    n_total here). Percentile/mean SEs and the Eq.176/177 convergence
    check still use mc_layer3's own _cluster_bootstrap_stats()/
    _min_realizations_*() -- with exactly one value per outer scenario,
    that machinery degenerates into a standard i.i.d. bootstrap over
    n_outer scenarios, which is the textbook-correct thing to do here."""
    peak_levels = np.array([r["peak_level"] for r in records])
    n_outer = meta["n_outer"]
    outer_peak_levels = [np.array([lvl]) for lvl in peak_levels]  # one "inner" value each

    p5, p50, p95 = np.percentile(peak_levels, [5, 50, 95])
    mean_level = float(np.mean(peak_levels))
    std_level = float(np.std(peak_levels, ddof=1)) if n_outer > 1 else 0.0
    thresholds = [t for t in (meta["max_flood_level"], meta["dam_crest_level"]) if t is not None]
    pctl_se, se_mean, threshold_se = _cluster_bootstrap_stats(
        outer_peak_levels, [5, 50, 95], thresholds)

    n_min_mean = _min_realizations_mean(std_level, mean_level)
    n_min_p95 = _min_realizations_percentile(0.95)

    mc_dir = os.path.join(OUTPUT_DIR, case_name, "OuterOnly")
    csv_path = os.path.join(mc_dir, "outer_only_results.csv")
    write_mc_results_csv(records, csv_path)

    plot_dir = os.path.join(PLOT_DIR, case_name, "OuterOnly")
    plot_path = os.path.join(plot_dir, "outer_only_distribution.png")
    write_mc_distribution_plot(
        peak_levels, plot_path,
        title=f"Layer 3 OUTER-ONLY - {case_name} ({meta['rule']} rule, "
              f"{n_outer} outer draws, baseline inner loop)",
        max_flood_level=meta["max_flood_level"], dam_crest_level=meta["dam_crest_level"],
        outer_peak_levels=outer_peak_levels,
    )

    summary_path = os.path.join(mc_dir, "outer_only_summary.txt")
    os.makedirs(mc_dir, exist_ok=True)
    with open(summary_path, "w") as f:
        f.write(f"Case: {case_name}\n")
        f.write(f"Rule tested: {meta['rule']}\n")
        f.write(f"Seed: {meta['seed']}\n")
        f.write(f"Draws: {n_outer} outer x 1 inner (DETERMINISTIC BASELINE -- see module "
                f"docstring) = {n_outer}\n")
        f.write(f"Inner loop: BASELINE ONLY -- every physical multiplier at its nominal "
                f"value, H0 unperturbed, area_scale=1.0, ALL gates operational (no failures "
                f"modeled). This isolates the OUTER loop's own contribution to spread -- "
                f"compare against a full Layer 3 run's spread to see how much of the total "
                f"uncertainty comes from climate/flood-frequency vs. physical/reliability "
                f"uncertainty.\n")
        f.write(f"Physical params held at baseline: {list(meta['uncertain_params'].keys())}\n")
        f.write(f"Outer-loop source (peak): {meta['outer_source']}\n")
        if meta["target_return_period"] is not None:
            f.write(f"Outer-loop target return period: {meta['target_return_period']} years "
                    f"(median hydrograph PEAK scale factor: {meta['median_scale']:.4f})\n")
        f.write(f"Outer-loop PEAK scale CV: {meta['outer_cv']:.4f}\n")
        f.write(f"Outer-loop source (volume): {meta['volume_source']}\n")
        f.write(f"Outer-loop median VOLUME scale factor: {meta['volume_median_scale']:.4f}, "
                f"CV: {meta['volume_cv']:.4f}\n")
        if meta["dependence_source"] == "ensemble":
            f.write(f"Outer-loop peak-volume DEPENDENCE: inherent in the ensemble replicates.\n")
        elif meta["dependence_source"] not in ("independent", None):
            f.write(f"Outer-loop peak-volume DEPENDENCE: mc_peak_volume_dependence="
                    f"'{meta['dependence_source']}', tau={meta['peak_volume_tau']:.3f}\n")
        else:
            f.write(f"Outer-loop peak-volume DEPENDENCE: independent (default)\n")

        f.write(f"\nPeak level [m a.s.l.]: P5={p5:.3f} (SE~={pctl_se[5]:.3f})  "
                f"P50={p50:.3f} (SE~={pctl_se[50]:.3f})  P95={p95:.3f} (SE~={pctl_se[95]:.3f})  "
                f"mean={mean_level:.3f} (SE~={se_mean:.3f})  max={peak_levels.max():.3f}  "
                f"min={peak_levels.min():.3f}\n")
        f.write(f"Convergence (RMC-TotalRisk manual Eq. 176/177, target: mean within +/-1%, "
                f"P95 within +/-0.01, both at 95% confidence), vs. n_outer={n_outer} "
                f"(= n_total here, no inner-loop nesting to correct for): "
                f"mean needs >={n_min_mean:.0f} realizations (this run: {n_outer}, "
                f"{_format_margin(n_outer, n_min_mean)} margin); "
                f"P95 needs >={n_min_p95:.0f} realizations (this run: {n_outer}, "
                f"{_format_margin(n_outer, n_min_p95)} margin)\n")
        if meta["max_flood_level"] is not None:
            frac = float(np.mean(peak_levels > meta["max_flood_level"]))
            f.write(f"P(peak level > max_flood_level={meta['max_flood_level']}): "
                    f"{_format_exceedance(frac, threshold_se.get(meta['max_flood_level'], 0.0), n_outer, n_outer)}\n")
        if meta["dam_crest_level"] is not None:
            frac = float(np.mean(peak_levels > meta["dam_crest_level"]))
            f.write(f"P(peak level > dam_crest_level={meta['dam_crest_level']}): "
                    f"{_format_exceedance(frac, threshold_se.get(meta['dam_crest_level'], 0.0), n_outer, n_outer)}\n")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {plot_path}")
    print(f"Wrote {summary_path}")
    print(f"\nPeak level [m a.s.l.]: P5={p5:.3f}  P50={p50:.3f}  P95={p95:.3f}  mean={mean_level:.3f}")


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
    parser = argparse.ArgumentParser(description="Layer 3 outer-loop-only sensitivity run")
    parser.add_argument("case_name", nargs="?", default=None)
    parser.add_argument("--n-outer", type=int, default=None)
    parser.add_argument("--rule", choices=["baseline", "optimized"], default="baseline")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    case_name = args.case_name or mc_layer3.prompt_for_case(mc_layer3.list_case_folders())
    n_outer = args.n_outer if args.n_outer is not None else prompt_int("Number of outer (climate) draws", 100)

    print(f"\nRunning Layer 3 outer-loop-only for '{case_name}': "
          f"{n_outer} outer draws, baseline inner loop, rule='{args.rule}'\n")

    records, meta = run_case_outer_only(case_name, n_outer, seed=args.seed, rule=args.rule)
    summarize_and_write(case_name, records, meta)


if __name__ == "__main__":
    main()
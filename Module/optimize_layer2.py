"""
optimize_layer2.py
-------------------
Lives in Module/. Run with:

    python Module/optimize_layer2.py <CaseName>

e.g.
    python Module/optimize_layer2.py Template

Layer 2: optimizes gate OPERATING RULE parameters -- whatever
parameters the case's own case_config.optimizable_gates() declares
per gate/group -- to protect downstream populations without
jeopardizing dam safety. Physical/rating parameters (crest levels,
widths, discharge coefficients, N-1 counts, ...) are never touched --
only the control logic that decides when/how far a gate opens. See the
"operating rule vs. discharge coefficient" section in outlets.py /
README.md for why that split matters.

The most common case is level_trigger_rule()'s (H_open, H_full, a_min,
a_max) per gate (see Data/Template/case_config.py). But the decision
variables are NOT hard-coded to that shape -- optimizable_gates() can
declare any parameter names a case's own build_outlets(rule_overrides)
knows how to use, e.g. a single "shift" that slides an entire
staged_trigger_rule() step schedule earlier/later in level while
preserving a real staggered multi-gate sequence. See build_problem()'s
docstring below for how this flattening works.

Problem formulation (agreed in HANDOFF_LAYER2.md)
--------------------------------------------------
Decision variables : whatever optimizable_gates() declares per group,
                      flattened into one vector for the optimizer.
Objective           : maximize the time before total downstream release
                      exceeds `downstream_threshold_m3s` (a scalar you
                      set in scalars.csv) -- "delay release as long as
                      possible". Framed for pymoo (a minimizer) as
                      minimizing NEGATIVE delay.
Hard constraint      : reservoir level must never exceed
                      `max_flood_level` from scalars.csv (NOT
                      dam_crest_level itself -- see README.md).
Secondary constraint : for any group that declares BOTH "H_open" and
                      "H_full" (or both "a_min" and "a_max"), enforce
                      the greater-than relationship by at least a small
                      margin as an explicit inequality constraint
                      rather than baking it into the variable bounds,
                      so those bounds can overlap freely. Groups using
                      a different parameterization (e.g. a shared
                      "shift") have nothing checked here -- any value
                      within their own bounds is inherently valid.
Optimizer            : pymoo's NSGA-II. With a single objective this
                      behaves like an elitist GA with constraint
                      handling; the Problem below is written so a
                      second objective (e.g. minimize peak downstream
                      release) is a one-line addition later if you want
                      a genuine Pareto front instead of a single answer.

Downstream release definition (please sanity-check for your dam)
------------------------------------------------------------------
"Downstream release" is taken as Qout_total (sum of all outlet/spillway
discharges) + Qwithdrawal (e.g. turbine offtake, if that water also
continues downstream rather than being diverted elsewhere). If your
withdrawal does NOT rejoin the downstream channel, edit
`downstream_release()` below to drop the Qwithdrawal term.

Performance note
-----------------
This bypasses run_case.py's CSV/case-folder loading (meant for one-off
human-run cases) and instead loads the case ONCE, then calls
solver.run_simulation() directly per optimizer evaluation -- each
evaluation is a fresh, cheap RK4 integration (outlets are rebuilt fresh
each time via case_config.build_outlets(rule_overrides=...), and
run_simulation() already calls outlet.reset() at the start, so no state
leaks between evaluations).
"""

from __future__ import annotations

import argparse
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
from report import write_results_csv, write_plot
PLOT_DIR = os.path.join(PROJECT_ROOT, "Plot")


# ---------------------------------------------------------------------
# Case loading (mirrors run_case.py's loading, done ONCE up front)
# ---------------------------------------------------------------------

def load_case_config(case_name: str):
    case_dir = os.path.join(DATA_DIR, case_name)
    config_path = os.path.join(case_dir, "case_config.py")
    if not os.path.isfile(config_path):
        raise SystemExit(f"'{case_name}' not found under {DATA_DIR}")
    spec = importlib.util.spec_from_file_location(f"case_config_{case_name}", config_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, case_dir


def downstream_release(result) -> np.ndarray:
    """Total flow continuing downstream at each recorded time step.
    See the module docstring's "Downstream release definition" note --
    edit this if withdrawal doesn't rejoin the downstream channel."""
    return result.Qout_total + result.Qwithdrawal


def time_to_exceedance(t_h: np.ndarray, series: np.ndarray, threshold: float) -> float:
    """First time [h] at which `series` exceeds `threshold`, linearly
    interpolated between the bracketing recorded points for a smoother
    (less step-size-dependent) objective. Returns t_h[-1] (i.e. "never
    exceeded within the simulated window") if the threshold is never
    crossed."""
    above = series > threshold
    if not above.any():
        return float(t_h[-1])
    i = int(np.argmax(above))  # first True
    if i == 0:
        return float(t_h[0])
    t0, t1 = t_h[i - 1], t_h[i]
    q0, q1 = series[i - 1], series[i]
    if q1 == q0:
        return float(t1)
    frac = (threshold - q0) / (q1 - q0)
    return float(t0 + frac * (t1 - t0))


# ---------------------------------------------------------------------
# pymoo Problem
# ---------------------------------------------------------------------

def build_problem(case_name: str, min_gate_gap: float = 0.2):
    """Returns (problem, meta) for the named case. min_gate_gap is the
    smallest allowed margin for (H_full - H_open) and (a_max - a_min)
    -- see the gap-constraint note below -- wherever a case's
    optimizable_gates() actually declares those specific parameter
    names; it has no effect on other parameter names (e.g. a case using
    staged_trigger_rule()'s "shift" parameter instead).

    Decision variables are NOT hard-coded to (H_open, H_full, a_min,
    a_max) -- each group in optimizable_gates() can declare WHATEVER
    parameter names its own build_outlets(rule_overrides=...) knows how
    to interpret (e.g. a single "shift" for a staged/stepped schedule).
    This function just flattens whatever's declared into one vector for
    pymoo and back."""
    from pymoo.core.problem import Problem

    case_config, case_dir = load_case_config(case_name)

    if not hasattr(case_config, "optimizable_gates"):
        raise SystemExit(
            f"Data/{case_name}/case_config.py has no optimizable_gates() -- "
            "either add one declaring which gates/groups Layer 2 may tune "
            "and their bounds (see Data/Template/case_config.py for the "
            "level_trigger_rule pattern, or README.md's 'Gated structures' "
            "section for the staged_trigger_rule 'shift' pattern), or this "
            "case's operating rule genuinely has no free parameter left to "
            "tune (e.g. a fully reactive/inflow-tracking rule) -- Layer 2 "
            "doesn't apply to every case."
        )

    gate_bounds = case_config.optimizable_gates()
    if not gate_bounds:
        raise SystemExit(f"optimizable_gates() in case '{case_name}' returned no gates.")

    gate_names = list(gate_bounds.keys())

    # Flatten every (group, param) pair across all groups, IN THE ORDER
    # each group's own dict declares them -- this is what makes the
    # decision vector's shape follow the case's own declaration instead
    # of assuming a fixed 4-parameter-per-gate layout.
    param_index = []  # list of (gname, pname), same order as xl/xu below
    xl, xu = [], []
    for gname in gate_names:
        for pname, (lo, hi) in gate_bounds[gname].items():
            param_index.append((gname, pname))
            xl.append(lo)
            xu.append(hi)
    xl = np.array(xl, dtype=float)
    xu = np.array(xu, dtype=float)
    n_var = len(xl)

    sc = case_config.scalars(case_dir)
    threshold = sc.get("downstream_threshold_m3s")
    if threshold is None:
        raise SystemExit(
            f"scalars.csv for case '{case_name}' has no downstream_threshold_m3s -- "
            "add a row, e.g.: downstream_threshold_m3s,100,m3/s"
        )
    max_flood_level = sc.get("max_flood_level")
    if max_flood_level is None:
        raise SystemExit(
            f"scalars.csv for case '{case_name}' has no max_flood_level -- "
            "Layer 2's hard safety constraint needs it (see README.md)."
        )

    reservoir = Reservoir(os.path.join(case_dir, "reservoir_curve.csv"))
    inflow = Hydrograph(os.path.join(case_dir, "inflow_hydrograph.csv"))
    withdrawal = Hydrograph(os.path.join(case_dir, "withdrawal_schedule.csv"))

    def decode(x: np.ndarray) -> dict:
        """Flat vector -> {group_name: {param_name: value, ...}, ...},
        following whatever parameter names optimizable_gates() declared
        for each group."""
        overrides: dict = {}
        for i, (gname, pname) in enumerate(param_index):
            overrides.setdefault(gname, {})[pname] = float(x[i])
        return overrides

    def simulate(overrides: dict | None):
        """Runs one full RK4 simulation for a given rule_overrides dict
        (or None for this case's own hand-set defaults, i.e. the
        baseline) and returns the full SimulationResult -- used both by
        run_with() (fast path, just needs the summary numbers) and by
        main() afterwards to write CSV/plot for the baseline and the
        optimized rule."""
        outlets = case_config.build_outlets(rule_overrides=overrides)
        return run_simulation(
            reservoir=reservoir, inflow=inflow, withdrawal=withdrawal,
            outlets=outlets, H0=sc["H0"], t_max=sc["t_max"], dt=sc["dt"],
            print_every=int(sc["print_every"]),
        )

    def run_with(x: np.ndarray):
        """One RK4 simulation for one candidate rule-parameter vector.
        Returns (peak_level, delay_h, gap_violation_max)."""
        overrides = decode(x)

        # H_full <= H_open, or a_min >= a_max, are not just
        # "undesirable" -- they're invalid input to level_trigger_rule()
        # (raises ValueError) -- NSGA-II's random sampling/crossover can
        # and will produce such points since the bounds in
        # optimizable_gates() are independent per parameter. Only
        # checked for groups that actually declare BOTH names of a pair
        # -- a group using a different parameterization (e.g. a shared
        # "shift") has nothing to check here, and any value
        # within its own bounds is inherently valid.
        gap_violation = 0.0
        invalid = False
        for gname in gate_names:
            params = overrides[gname]
            if "H_open" in params and "H_full" in params:
                level_gap = params["H_full"] - params["H_open"]
                gap_violation = max(gap_violation, min_gate_gap - level_gap)
                if level_gap <= 0.0:
                    invalid = True
            if "a_min" in params and "a_max" in params:
                opening_gap = params["a_max"] - params["a_min"]
                gap_violation = max(gap_violation, min_gate_gap - opening_gap)
                if opening_gap <= 0.0:
                    invalid = True
        if invalid:
            return max_flood_level + 1e3, 0.0, gap_violation

        result = simulate(overrides)
        peak_level = float(result.H.max())
        delay_h = time_to_exceedance(result.t_h, downstream_release(result), threshold)

        return peak_level, delay_h, gap_violation

    class RuleCurveProblem(Problem):
        """
        Variables : flattened parameters, whatever optimizable_gates()
                    declared per group (see build_problem docstring).
        Objective : minimize -(delay to threshold exceedance) [h]
                    (i.e. maximize delay).
        Constraints (both must be <= 0 to be feasible):
            g1 = peak_level - max_flood_level          (no overtopping
                                                         of the design level)
            g2 = max over all gates of:                 (keeps H_full >
                 min_gate_gap - (H_full - H_open),        H_open and
                 min_gate_gap - (a_max - a_min)            a_max > a_min,
                                                            each with a
                                                            safety margin)
        """

        def __init__(self):
            super().__init__(n_var=n_var, n_obj=1, n_ieq_constr=2,
                              xl=xl, xu=xu)

        def _evaluate(self, X, out, *args, **kwargs):
            n_pop = X.shape[0]
            f = np.zeros((n_pop, 1))
            g = np.zeros((n_pop, 2))
            for i in range(n_pop):
                peak_level, delay_h, gap_violation = run_with(X[i])
                f[i, 0] = -delay_h
                g[i, 0] = peak_level - max_flood_level
                g[i, 1] = gap_violation
            out["F"] = f
            out["G"] = g

    problem = RuleCurveProblem()
    meta = {
        "gate_names": gate_names,
        "param_index": param_index,
        "decode": decode,
        "run_with": run_with,
        "simulate": simulate,
        "threshold": threshold,
        "max_flood_level": max_flood_level,
        "dam_crest_level": sc.get("dam_crest_level"),
        "case_dir": case_dir,
    }
    return problem, meta


def list_case_folders():
    """Case names = subfolders of Data/ that contain a case_config.py.
    Mirrors run_case.py's list_case_folders()."""
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
        raise SystemExit(
            f"No case folders found under {DATA_DIR}. "
            "Each case needs its own subfolder containing case_config.py."
        )
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


# ---------------------------------------------------------------------
# Run NSGA-II and report the result
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Layer 2: gate operating-rule optimization")
    parser.add_argument("case_name", nargs="?", default=None,
                         help="Case folder under Data/, e.g. Template. "
                              "Omit to be prompted interactively.")
    parser.add_argument("--pop-size", type=int, default=None,
                         help="Population size (prompted with a default if omitted)")
    parser.add_argument("--n-gen", type=int, default=None,
                         help="Number of generations (prompted with a default if omitted)")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    case_name = args.case_name or prompt_for_case(list_case_folders())
    pop_size = args.pop_size if args.pop_size is not None else prompt_int("Population size", 40)
    n_gen = args.n_gen if args.n_gen is not None else prompt_int("Number of generations", 30)

    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.optimize import minimize

    problem, meta = build_problem(case_name)

    algorithm = NSGA2(pop_size=pop_size)

    print(f"\nOptimizing gate rule(s) {meta['gate_names']} for case '{case_name}'")
    print(f"  downstream_threshold_m3s = {meta['threshold']}")
    print(f"  max_flood_level          = {meta['max_flood_level']}")
    print(f"  pop_size={pop_size}, n_gen={n_gen} "
          f"({pop_size * n_gen} simulations total)\n")

    res = minimize(problem, algorithm, ("n_gen", n_gen),
                    seed=args.seed, verbose=True)

    if res.X is None:
        print("\nNo feasible solution found (every candidate violated a constraint --"
              " likely max_flood_level, or the search bounds/n_gen are too tight)."
              " Try widening optimizable_gates() bounds or increasing --n-gen.")
        return

    # res.X is 1D if a single best solution was found (n_obj=1 -> pymoo
    # returns the single best feasible individual, not a Pareto front).
    x_best = res.X if res.X.ndim == 1 else res.X[0]
    peak_level, delay_h, gap_violation = meta["run_with"](x_best)
    overrides = meta["decode"](x_best)

    print("\nBest feasible rule found:")
    for gname, params in overrides.items():
        param_str = ", ".join(f"{p}={v:.3f}" for p, v in params.items())
        print(f"  {gname}: {param_str}")
    print(f"  -> peak reservoir level : {peak_level:.3f} m a.s.l. "
          f"(limit {meta['max_flood_level']:.3f} m)")
    print(f"  -> downstream threshold ({meta['threshold']} m3/s) reached at "
          f"t = {delay_h:.2f} h")

    # Baseline numbers (case_config.py's own hand-set rule) are still
    # useful context for the comparison below, but writing the
    # Baseline/ CSV+plot is run_case.py's job, not this script's --
    # run `python Module/run_case.py <CaseName>` to (re)generate it.
    baseline_result = meta["simulate"](None)
    baseline_peak = float(baseline_result.H.max())
    baseline_delay = time_to_exceedance(
        baseline_result.t_h, downstream_release(baseline_result), meta["threshold"])

    # --- Optimised: the best rule NSGA-II found. ---
    optimised_result = meta["simulate"](overrides)
    optimised_out_dir = os.path.join(OUTPUT_DIR, case_name, "Optimised")
    write_results_csv(optimised_result, os.path.join(optimised_out_dir, "results.csv"))
    optimised_plot_dir = os.path.join(PLOT_DIR, case_name, "Optimised")
    write_plot(optimised_result, os.path.join(optimised_plot_dir, "routing_result.png"),
               title=f"Reservoir routing - {case_name} (Optimised)",
               max_flood_level=meta["max_flood_level"],
               dam_crest_level=meta["dam_crest_level"],
               downstream_threshold=meta["threshold"])

    best_rule_path = os.path.join(optimised_out_dir, "best_rule.txt")
    with open(best_rule_path, "w") as fh:
        fh.write(f"Case: {case_name}\n")
        fh.write(f"downstream_threshold_m3s: {meta['threshold']}\n")
        fh.write(f"max_flood_level: {meta['max_flood_level']}\n\n")
        for gname, params in overrides.items():
            param_str = ", ".join(f"{p}={v:.3f}" for p, v in params.items())
            fh.write(f"{gname}: {param_str}\n")
        fh.write(f"\npeak_level: {peak_level:.3f}\n")
        fh.write(f"delay_to_threshold_h: {delay_h:.3f}\n")
        fh.write(f"\n--- Baseline (case_config.py's own hand-set rule) for comparison ---\n")
        fh.write(f"baseline_peak_level: {baseline_peak:.3f}\n")
        fh.write(f"baseline_delay_to_threshold_h: {baseline_delay:.3f}\n")

    print(f"\nBaseline (not re-written here -- run `python run_case.py {case_name}` "
          f"for Output/Plot/{case_name}/Baseline/): "
          f"peak {baseline_peak:.2f} m, threshold reached at t={baseline_delay:.2f} h")
    print(f"Optimised : Output/{case_name}/Optimised/results.csv, "
          f"Plot/{case_name}/Optimised/routing_result.png "
          f"(peak {peak_level:.2f} m, threshold reached at t={delay_h:.2f} h)")
    print(f"Rule summary written to {best_rule_path}")


if __name__ == "__main__":
    main()
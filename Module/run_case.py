"""
run_case.py
-----------
Lives in Module/. Run with:

    python Module/run_case.py

Behavior:
    1. Lists the case folders found under Data/ (e.g. Data/Template/,
       Data/DesignFlood_2050/, ...) and asks you to pick one.
    2. Loads that case's reservoir curve, hydrograph, withdrawal
       schedule, scalars, and outlet definitions (case_config.py) from
       Data/<CaseName>/.
    3. Runs the RK4 routing.
    4. Writes results to Output/<CaseName>/results.csv and a plot to
       Plot/<CaseName>/routing_result.png -- i.e. the SAME case name,
       mirrored under Output/ and Plot/.

You can also skip the interactive prompt by passing the case name
directly, e.g.:
    python Module/run_case.py Template
which is useful later for scripted/batch runs (optimization, Monte
Carlo) without a human in the loop.
"""

import sys
import os
import importlib.util

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


def list_case_folders():
    """Case names = subfolders of Data/ that contain a case_config.py.
    Each case is fully self-contained in its own Data/<name>/ folder
    (CSVs + case_config.py) -- no separate matching-by-name file
    elsewhere, so creating a new case is just: duplicate the folder,
    rename it, edit its contents."""
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


def load_case_config(case_name: str):
    """Dynamically imports Data/<case_name>/case_config.py, whatever the
    case name is -- so cases can be added/renamed without editing any
    import statements anywhere else."""
    case_dir = os.path.join(DATA_DIR, case_name)
    config_path = os.path.join(case_dir, "case_config.py")
    spec = importlib.util.spec_from_file_location(f"case_config_{case_name}", config_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, case_dir


def main():
    if len(sys.argv) > 1:
        case_name = sys.argv[1]
        if case_name not in list_case_folders():
            raise SystemExit(f"'{case_name}' not found under {DATA_DIR}")
    else:
        case_name = prompt_for_case(list_case_folders())

    case_config, case_dir = load_case_config(case_name)

    sc = case_config.scalars(case_dir)
    reservoir = Reservoir(os.path.join(case_dir, "reservoir_curve.csv"))
    inflow = Hydrograph(os.path.join(case_dir, "inflow_hydrograph.csv"))
    withdrawal = Hydrograph(os.path.join(case_dir, "withdrawal_schedule.csv"))
    outlets = case_config.build_outlets()

    result = run_simulation(
        reservoir=reservoir,
        inflow=inflow,
        withdrawal=withdrawal,
        outlets=outlets,
        H0=sc["H0"],
        t_max=sc["t_max"],
        dt=sc["dt"],
        print_every=int(sc["print_every"]),
    )

    # --- Output/<case_name>/Baseline/results.csv ---
    # This runs case_config.py's own hand-set rule (build_outlets() with
    # no rule_overrides) -- the "Baseline" folder name pairs with
    # optimize_layer2.py's "Optimised" folder for direct comparison.
    out_dir = os.path.join(OUTPUT_DIR, case_name, "Baseline")
    out_path = os.path.join(out_dir, "results.csv")
    write_results_csv(result, out_path)

    i_peak = int(result.H.argmax())
    peak_H = result.H[i_peak]
    peak_t = result.t_h[i_peak]
    print(f"\nCase: {case_name}")
    print(f"Peak level: {peak_H:.2f} m a.s.l. at t={peak_t:.2f} h")
    print(f"Results written to {out_path}")

    # --- Safety threshold check ---
    # dam_crest_level: true physical top of the dam (actual overtopping
    #   if reached).
    # max_flood_level: stricter design/regulatory maximum reservoir
    #   level, kept below the crest by a freeboard margin -- this is
    #   the level that should never be reached even under the design
    #   flood, precisely so real overtopping never gets close.
    # Both are optional (sc.get(...) returns None if not defined in
    # scalars.csv), so this doesn't break cases that haven't set them.
    dam_crest_level = sc.get("dam_crest_level")
    max_flood_level = sc.get("max_flood_level")

    if dam_crest_level is not None and max_flood_level is not None:
        if max_flood_level >= dam_crest_level:
            print(f"WARNING: max_flood_level ({max_flood_level:.2f}) is not below "
                  f"dam_crest_level ({dam_crest_level:.2f}) -- check scalars.csv, "
                  f"this leaves no freeboard margin.")

    if dam_crest_level is not None and peak_H > dam_crest_level:
        print(f"*** OVERTOPPING: peak level {peak_H:.2f} m exceeds the dam crest "
              f"({dam_crest_level:.2f} m) by {peak_H - dam_crest_level:.2f} m ***")
    elif max_flood_level is not None and peak_H > max_flood_level:
        print(f"WARNING: peak level {peak_H:.2f} m exceeds the max flood/design level "
              f"({max_flood_level:.2f} m) by {peak_H - max_flood_level:.2f} m "
              f"-- design safety margin violated, even though the dam crest itself "
              f"was not reached.")
    elif max_flood_level is not None:
        print(f"OK: peak level is {max_flood_level - peak_H:.2f} m below the max "
              f"flood/design level ({max_flood_level:.2f} m).")
    elif dam_crest_level is not None:
        print(f"OK: peak level is {dam_crest_level - peak_H:.2f} m below the dam "
              f"crest ({dam_crest_level:.2f} m). Consider also setting "
              f"max_flood_level in scalars.csv for a stricter design margin.")
    else:
        print("Note: dam_crest_level / max_flood_level not set in scalars.csv -- "
              "no safety threshold check performed.")

    # --- Plot/<case_name>/Baseline/routing_result.png ---
    plot_dir = os.path.join(PLOT_DIR, case_name, "Baseline")
    plot_path = os.path.join(plot_dir, "routing_result.png")
    if write_plot(result, plot_path, title=f"Reservoir routing - {case_name} (Baseline)",
                  max_flood_level=max_flood_level, dam_crest_level=dam_crest_level):
        print(f"Plot written to {plot_path}")

    return result


if __name__ == "__main__":
    main()

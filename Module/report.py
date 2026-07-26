"""
report.py
---------
Shared result-writing code (CSV + plot), used by both run_case.py (a
single human-run case) and optimize_layer2.py (baseline vs. optimized
rule comparison), so both write results in the exact same format and
stay comparable rather than drifting apart as two separate
implementations.
"""

from __future__ import annotations
import csv
import os


def write_results_csv(result, out_path: str) -> None:
    """Writes the full time series to out_path, same layout as the
    original run_case.py results.csv: time, level, inflow, withdrawal,
    per-outlet discharge, total discharge, rate of rise."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    outlet_names = list(result.Qout_by_outlet.keys())
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_s", "time_h", "H_masl", "Qin", "Qwithdrawal"]
                    + [f"Qs_{n}" for n in outlet_names] + ["Qs_total", "Vmont_cmh"])
        for i in range(len(result.t_s)):
            w.writerow([
                result.t_s[i], result.t_h[i], result.H[i], result.Qin[i], result.Qwithdrawal[i],
                *[result.Qout_by_outlet[n][i] for n in outlet_names],
                result.Qout_total[i], result.V_rate_cmh[i],
            ])


def write_plot(result, plot_path: str, title: str,
                max_flood_level: float | None = None,
                dam_crest_level: float | None = None,
                downstream_threshold: float | None = None) -> bool:
    """Writes the routing plot to plot_path. Returns False (and prints
    a note) instead of raising if matplotlib isn't installed, matching
    run_case.py's original behavior of not hard-failing on a missing
    optional plotting dependency."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plot")
        return False

    os.makedirs(os.path.dirname(plot_path), exist_ok=True)

    fig, ax1 = plt.subplots(figsize=(9, 6))
    ax1.plot(result.t_h, result.H, color="tab:blue", label="Reservoir level H", linewidth=2)
    if max_flood_level is not None:
        ax1.axhline(max_flood_level, color="darkorange", linestyle="--",
                    linewidth=1.3, label="Max flood/design level")
    if dam_crest_level is not None:
        ax1.axhline(dam_crest_level, color="darkred", linestyle="-",
                    linewidth=1.3, label="Dam crest")
    ax1.set_xlabel("Time [h]")
    ax1.set_ylabel("Level [m a.s.l.]", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    ax2 = ax1.twinx()
    ax2.plot(result.t_h, result.Qin, color="tab:red", label="Qin", linestyle="--")
    ax2.plot(result.t_h, result.Qwithdrawal, color="0.5", label="Qwithdrawal", linestyle=":")

    outlet_names = list(result.Qout_by_outlet.keys())
    cmap = plt.get_cmap("tab10")
    for i, oname in enumerate(outlet_names):
        ax2.plot(result.t_h, result.Qout_by_outlet[oname], label=f"Qs {oname}",
                  color=cmap(i % 10), linewidth=1.3)

    ax2.plot(result.t_h, result.Qout_total, color="black", label="Qs total",
              linestyle="-.", linewidth=1.8)
    if downstream_threshold is not None:
        ax2.axhline(downstream_threshold, color="purple", linestyle="--",
                    linewidth=1.3, label="Downstream threshold")
    ax2.set_ylabel("Discharge [m3/s]", color="tab:red")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    fig.legend(lines1 + lines2, labels1 + labels2, loc="upper center",
               bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False)
    plt.title(title)
    fig.tight_layout()

    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ---------------------------------------------------------------------
# Layer 3 (Monte Carlo / climate-uncertainty) writers
# ---------------------------------------------------------------------
# The writers above assume a SINGLE run (one time series). Layer 3
# produces a DISTRIBUTION of results (thousands of scalar summaries,
# one per Monte Carlo draw), so it needs a different shape of output:
# one row per draw rather than one row per timestep.

def write_mc_results_csv(records: list[dict], out_path: str) -> None:
    """Writes one row per Monte Carlo draw. `records` is a list of flat
    dicts (see Module/mc_layer3.py's DrawResult) -- every dict must
    share the same set of keys (the first record's keys are used as
    the column order)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if not records:
        with open(out_path, "w", newline="") as f:
            f.write("")
        return
    fieldnames = list(records[0].keys())
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(records)


def write_mc_distribution_plot(peak_levels, plot_path: str, title: str,
                                max_flood_level: float | None = None,
                                dam_crest_level: float | None = None,
                                outer_peak_levels: "list | None" = None) -> bool:
    """Two-panel summary of the peak-reservoir-level distribution across
    all Monte Carlo draws:
      left  -- histogram of peak_level across every draw (all outer x
               inner combinations pooled)
      right -- empirical exceedance curve (peak_level vs. probability
               of exceedance), on a log probability axis, the
               conventional way to read off "what level is exceeded
               with probability p" -- with max_flood_level/
               dam_crest_level marked for a direct safety-margin read.
    If `outer_peak_levels` is given (list of arrays, one per outer/
    climate-scenario draw), each outer draw's own exceedance curve is
    drawn as a thin line UNDER the pooled curve, so the spread BETWEEN
    them (climate/frequency-curve uncertainty) is visually distinct
    from the spread WITHIN one (physical/rating uncertainty) -- the
    same separation of uncertainty sources the sampling design keeps
    throughout Layer 3.
    """
    import numpy as np
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plot")
        return False

    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    peak_levels = np.asarray(peak_levels)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    ax1.hist(peak_levels, bins=40, color="tab:blue", alpha=0.8)
    if max_flood_level is not None:
        ax1.axvline(max_flood_level, color="darkorange", linestyle="--",
                    linewidth=1.3, label="Max flood/design level")
    if dam_crest_level is not None:
        ax1.axvline(dam_crest_level, color="darkred", linestyle="-",
                    linewidth=1.3, label="Dam crest")
    ax1.set_xlabel("Peak reservoir level [m a.s.l.]")
    ax1.set_ylabel("Count (all draws pooled)")
    ax1.legend(frameon=False)
    ax1.set_title("Peak level distribution")

    def exceedance_curve(vals):
        s = np.sort(np.asarray(vals))[::-1]
        n = len(s)
        p = (np.arange(1, n + 1)) / (n + 1)  # plotting-position exceedance prob
        return s, p

    if outer_peak_levels:
        for arr in outer_peak_levels:
            if len(arr) == 0:
                continue
            s, p = exceedance_curve(arr)
            ax2.plot(p, s, color="tab:blue", alpha=0.25, linewidth=0.8)

    s, p = exceedance_curve(peak_levels)
    ax2.plot(p, s, color="black", linewidth=2, label="Pooled (all draws)")
    if max_flood_level is not None:
        ax2.axhline(max_flood_level, color="darkorange", linestyle="--",
                    linewidth=1.3, label="Max flood/design level")
    if dam_crest_level is not None:
        ax2.axhline(dam_crest_level, color="darkred", linestyle="-",
                    linewidth=1.3, label="Dam crest")
    ax2.set_xscale("log")
    ax2.set_xlabel("Exceedance probability (within this ensemble)")
    ax2.set_ylabel("Peak reservoir level [m a.s.l.]")
    ax2.legend(frameon=False)
    ax2.set_title("Exceedance curve"
                  + (" -- thin lines = individual climate/outer draws"
                     if outer_peak_levels else ""))

    plt.suptitle(title)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True

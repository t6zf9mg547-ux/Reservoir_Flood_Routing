"""
importance.py
--------------
Lives in Module/reliability/. Answers "what's actually driving this
probability" -- the reason this engine exists at all, per the point
that a screening tool whose only output is a number is far less useful
than one that also says where to focus a limited inspection/repair
budget.

METHOD (deliberately simple, not a full Birnbaum/Fussell-Vesely
importance analysis): for each named subsystem TYPE (e.g. "hydraulic",
"structural", "control_electrical" -- aggregated across every gate
that has one) and each common-cause BRANCH, sum that item's
probability contribution across the fleet, then report each as a
fraction of the total. This is a first-pass ranking, not a rigorous
sensitivity measure -- it answers "which category has the most total
probability mass behind it" rather than "how much would the overall
P_FOD actually drop if this item were fixed" (the latter needs
re-running compute_gate_pfods() with that item's Pf set to zero and
comparing -- a natural follow-on once this ranking identifies which
items are worth that closer look).
"""

from __future__ import annotations
from collections import defaultdict


def rank_importance(component_pfs: dict[str, dict[str, float]],
                     ccf_groups: list[dict]) -> list[dict]:
    """Returns a list of {"label": str, "contribution": float,
    "share": float} dicts, sorted by contribution descending.

    component_pfs : {gate_id: {subsystem: Pf}}, as returned by
        gate_reliability.compute_gate_pfods().
    ccf_groups : list of {"label", "gates", "p_ccf", ...}, as returned
        by common_cause.load_ccf_groups().
    """
    by_subsystem: dict[str, float] = defaultdict(float)
    for gate_id, subsystem_pfs in component_pfs.items():
        for subsystem, pf in subsystem_pfs.items():
            by_subsystem[subsystem] += pf

    items = [{"label": f"Subsystem: {name}", "contribution": total}
              for name, total in by_subsystem.items()]
    items += [{"label": f"Common-cause: {grp['label']}", "contribution": grp["p_ccf"]}
              for grp in ccf_groups]

    total = sum(item["contribution"] for item in items) or 1.0
    for item in items:
        item["share"] = item["contribution"] / total

    return sorted(items, key=lambda x: x["contribution"], reverse=True)


def format_importance_table(ranked: list[dict]) -> str:
    """Plain-text table for mc_summary.txt / console output."""
    lines = [f"{'Item':<45} {'Contribution':>14} {'Share':>8}"]
    for item in ranked:
        lines.append(f"{item['label']:<45} {item['contribution']:>14.4f} "
                      f"{item['share']*100:>7.1f}%")
    return "\n".join(lines)

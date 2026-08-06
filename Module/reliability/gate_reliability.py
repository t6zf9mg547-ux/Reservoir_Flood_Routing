"""
gate_reliability.py
--------------------
Lives in Module/reliability/. Reads a case's components.csv, computes
one failure-on-demand probability (P_FOD) per gate, and assembles the
complete dict that a case's own mc_gate_availability() hook returns to
Module/mc_layer3.py.

This is the ONE function most case_config.py files will actually call:

    from reliability.gate_reliability import build_gate_availability

    def mc_gate_availability():
        case_dir = os.path.dirname(os.path.abspath(__file__))
        return build_gate_availability(
            os.path.join(case_dir, "reliability", "components.csv"),
            os.path.join(case_dir, "reliability", "common_cause_events.csv"),
        )

Everything upstream of this (mc_layer3.py's sampling, mc_results.csv,
mc_summary.txt) is unchanged by whether a case populates this from
real CI evidence for whatever case supplies it, or from a single
flat, illustrative guess (see Data/Template/) -- both are valid inputs
to the SAME mc_gate_availability() contract; see the reliability
README for the two-tier explanation.
"""

from __future__ import annotations
import csv
import os
from collections import defaultdict

from .condition_index import pf_from_ci
from .fault_tree import gate_pfod
from .common_cause import load_ccf_groups


def compute_component_pfs(components_csv: str) -> dict[str, dict[str, float]]:
    """Reads components.csv and returns {gate_id: {subsystem: Pf}},
    i.e. one failure probability per (gate, subsystem) row.

    Expected columns: component_id, gate_id, subsystem, CI_mean,
    CI_sd, CI_failure, notes. One row per subsystem per gate (a
    case that wants finer-grained components -- e.g. splitting
    "hydraulic" into "pump" and "cylinder" separately -- can do so;
    every distinct subsystem value for a gate is combined in series,
    so a finer breakdown is a strictly more detailed input, not a
    different contract).
    """
    by_gate: dict[str, dict[str, float]] = defaultdict(dict)
    with open(components_csv, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gate_id = row["gate_id"].strip()
            subsystem = row["subsystem"].strip()
            pf = pf_from_ci(float(row["CI_mean"]), float(row["CI_sd"]),
                             float(row["CI_failure"]))
            by_gate[gate_id][subsystem] = pf
    return dict(by_gate)


def compute_gate_pfods(components_csv: str) -> tuple[dict[str, float], dict[str, dict[str, float]]]:
    """Returns (p_fail_by_gate, component_pfs):
      p_fail_by_gate -- {gate_id: P_FOD}, one series-combined
        probability per gate, independent of common-cause effects
        (those are layered on separately -- see build_gate_availability).
      component_pfs -- {gate_id: {subsystem: Pf}}, the per-subsystem
        breakdown compute_component_pfs() already produced, passed
        through unchanged so importance.py doesn't have to re-read
        the CSV.
    """
    component_pfs = compute_component_pfs(components_csv)
    p_fail_by_gate = {gate_id: gate_pfod(subsystem_pfs)
                       for gate_id, subsystem_pfs in component_pfs.items()}
    return p_fail_by_gate, component_pfs


def build_gate_availability(components_csv: str, common_cause_csv: str | None = None) -> dict:
    """Assembles the complete dict a case's mc_gate_availability() hook
    should return: {"gates": [...], "p_fail_by_gate": {...},
    "ccf_groups": [...]}.

    components_csv is required. common_cause_csv is optional -- if
    None, or if the file doesn't exist, ccf_groups is simply empty
    (no common-cause mechanism modeled, same as a case that has no
    documented shared vulnerability).
    """
    p_fail_by_gate, _ = compute_gate_pfods(components_csv)
    gates = list(p_fail_by_gate.keys())
    ccf_groups = load_ccf_groups(common_cause_csv, gates) if common_cause_csv else []
    return {
        "gates": gates,
        "p_fail_by_gate": p_fail_by_gate,
        "ccf_groups": ccf_groups,
    }

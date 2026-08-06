"""
common_cause.py
----------------
Lives in Module/reliability/. Turns a case's common_cause_events.csv
into a list of CCF (common-cause failure) groups: named, shared
vulnerabilities that, if triggered, are assumed to affect every gate
in the named group at once -- e.g. a shared control link, a shared
component model/batch across every gate's hoist, a shared power
source.

MODEL: the mixture model from the reference methodology --

    P(all gates in the group fail) = P_CCF + (1 - P_CCF) * P_gate**n

-- but this module only produces the P_CCF term per group (from that
branch's own CI estimate, via condition_index.pf_from_ci). The actual
MIXING -- drawing whether a CCF event fires, and if so forcing every
gate in its group to the failed state regardless of that gate's own
independent draw -- happens once per Monte Carlo inner draw, in
Module/mc_layer3.py's default_inner_sampler(), not here. This module
is pure data transformation: CSV row -> {gates, p_ccf, label}.

Each named branch is independent of the others and of every gate's own
individual failure draw -- if a case has two real, distinct shared
vulnerabilities (e.g. a comms link AND a shared component batch), they
belong as two separate rows/groups, not blended into one number.
"""

from __future__ import annotations
import csv

from .condition_index import pf_from_ci


def load_ccf_groups(common_cause_csv: str, all_gates: list[str]) -> list[dict]:
    """Reads common_cause_events.csv and returns a list of
    {"label": str, "gates": list[str], "p_ccf": float} dicts, one per
    row.

    Expected columns: event_id, label, affected_gates, CI_mean, CI_sd,
    CI_failure, notes. affected_gates is either the literal string
    "ALL" (expands to every gate the case declares) or a semicolon-
    separated list of gate names matching the case's own gate IDs
    (e.g. "gate_1;gate_3;gate_5").

    Returns an empty list if the file doesn't exist or has no rows --
    a case with no documented common-cause vulnerability simply skips
    this, exactly like a case with no mc_uncertain_params() skips
    physical-coefficient perturbation.
    """
    import os
    if not os.path.exists(common_cause_csv):
        return []

    groups = []
    with open(common_cause_csv, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            affected = row["affected_gates"].strip()
            gates = list(all_gates) if affected.upper() == "ALL" else \
                [g.strip() for g in affected.split(";") if g.strip()]
            ci_mean = float(row["CI_mean"])
            ci_sd = float(row["CI_sd"])
            ci_failure = float(row["CI_failure"])
            p_ccf = pf_from_ci(ci_mean, ci_sd, ci_failure)
            groups.append({
                "event_id": row.get("event_id", ""),
                "label": row.get("label", row.get("event_id", "")),
                "gates": gates,
                "p_ccf": p_ccf,
                "ci_mean": ci_mean, "ci_sd": ci_sd, "ci_failure": ci_failure,
                "notes": row.get("notes", ""),
            })
    return groups

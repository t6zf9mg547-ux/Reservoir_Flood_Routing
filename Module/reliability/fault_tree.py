"""
fault_tree.py
-------------
Lives in Module/reliability/. Combines a list of independent subsystem
(or component) failure probabilities into one top-event probability,
using the two combination rules a gate-level fault tree needs.

SERIES logic ("gate fails to open if ANY essential subsystem fails")
is the default and correct choice for most spillway-gate subsystems --
hydraulics, structure, and local control are each individually
necessary for the gate to open, so there's no redundancy to credit
unless a case explicitly models one as genuinely independent and
available (see the reliability README for that caution). PARALLEL
logic is provided for completeness (a case with a real, independent,
tested backup for one specific subsystem), but is not exercised by the
default gate-level combination below.
"""

from __future__ import annotations


def series_pf(pfs: list[float]) -> float:
    """Failure probability of a series (any-one-fails) combination.

    P_fail = 1 - product(1 - Pf_i)

    Use for essential subsystems where each one independently prevents
    the gate from opening if it fails -- the default, conservative
    choice absent evidence of genuine, independent redundancy.
    """
    reliability = 1.0
    for pf in pfs:
        reliability *= (1.0 - pf)
    return 1.0 - reliability


def parallel_pf(pfs: list[float]) -> float:
    """Failure probability of a parallel (all-must-fail) combination.

    P_fail = product(Pf_i)

    Valid ONLY if the elements are genuinely independent -- a shared
    bus, shared cabinet, or shared control path between them means
    this overstates reliability; model that as a common-cause branch
    instead (see common_cause.py).
    """
    pf = 1.0
    for item in pfs:
        pf *= item
    return pf


def gate_pfod(subsystem_pfs: dict[str, float]) -> float:
    """Convenience wrapper: one gate's P_FOD from its named subsystem
    failure probabilities (e.g. {"hydraulic": 0.03, "structural": 0.01,
    "control_electrical": 0.08}), combined in series -- the standard
    "any essential subsystem fails => gate fails to open" logic."""
    return series_pf(list(subsystem_pfs.values()))

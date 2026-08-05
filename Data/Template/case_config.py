"""
case_config.py
---------------
Lives inside Data/<CaseName>/ (e.g. Data/Template/case_config.py).
THIS is the file you edit to describe a specific dam/reservoir case:
which outlets exist, their levels, coefficients, and operating rules.
Everything else (reservoir curve, hydrograph, timestep) is read from
the plain CSV files in this same folder.

To create a new scenario: duplicate the whole Data/Template/ folder,
rename the duplicate (e.g. Data/MyCase/), then edit the CSVs and this
file inside it. The filename never changes -- it's always
"case_config.py" -- run_case.py finds it by folder, not by name
matching, so there's no renaming of this file itself, ever.

WHAT'S BELOW is a worked example of ALL FOUR outlet types available in
outlets.py, with illustrative placeholder numbers -- not calibrated
data for any real dam. Delete/edit/duplicate these to match your own
case; keep whichever combination of outlet types your dam actually has.
"""

import sys
import os

# run_case.py already puts Module/ on sys.path before loading this file,
# so this line is only a fallback for running/importing case_config.py
# standalone. Path depth: Data/<CaseName>/case_config.py -> up twice to
# project root -> Module/.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "Module"))

from outlets import (FreeOverflowSpillway, GatedSpillway, FuseGate,
                      BottomOutlet, level_trigger_rule)


def build_outlets(rule_overrides: dict | None = None,
                   physical_overrides: dict | None = None):
    """Returns the list of Outlet objects for this case. Called fresh
    for each simulation run so outlet internal state (e.g. fuse gates)
    always starts reset.

    rule_overrides : optional dict, {gate_name: {"H_open":.., "H_full":..,
        "a_min":.., "a_max":..}}, used by Module/optimize_layer2.py to
        try different operating-rule parameters WITHOUT touching any
        physical/rating parameter (crest levels, widths, coefficients,
        N-1 counts, etc). Only gates listed in optimizable_gates() below
        are expected to appear here. Leave None (the default) to get
        this case's normal, hand-set operating rule -- nothing changes
        for a plain `python run_case.py` run.

    physical_overrides : optional dict of MULTIPLIERS (1.0 = nominal,
        unchanged) on discharge/rating coefficients -- added for
        Module/mc_layer3.py's inner (physical-uncertainty) loop, kept
        strictly separate from rule_overrides above (same "operating
        rule vs. discharge coefficient" split used throughout this
        tool). Recognized keys, all optional:
            "spillway_1_C_mult", "gate_1_Cd_mult", "gate_1_free_flow_C_mult",
            "fuse_gate_1_pre_C_mult", "fuse_gate_1_post_C_mult",
            "bottom_outlet_1_Cd_mult"
        Leave None for a normal run (no change to any coefficient).
    """
    rule_overrides = rule_overrides or {}
    physical_overrides = physical_overrides or {}
    po = physical_overrides  # short alias, read throughout below

    # -------------------------------------------------------------
    # 1) Free-overflow (ungated) spillway, with pier contraction
    #    Q = C * Leff * (H - crest)^1.5
    #    Leff = L - 2*N_piers*Kp*(H - crest)
    # -------------------------------------------------------------
    spillway_1 = FreeOverflowSpillway(
        name="spillway_1",
        crest_level=100.0,   # [m a.s.l.] -- EDIT to your spillway's crest elevation
        length=20.0,         # [m] GROSS crest length (before pier reduction)
        C=2.0 * po.get("spillway_1_C_mult", 1.0),  # [m^0.5/s] weir
                             # coefficient -- typically ~1.7-2.2 for a
                             # broad-crested weir; refine against your
                             # own hydraulic design/model data
        N_piers=2,           # number of piers spanning the crest (0 = none)
        Kp=0.01,             # pier contraction coefficient -- EDIT to match
                             # your pier nose shape / design data
    )

    # -------------------------------------------------------------
    # 2) Gated spillway (underflow/orifice gate), with a free-flow
    #    transition once the gate is raised clear of the flow path
    #    Q = Cd(H, a) * width * a * sqrt(2*g*H_eff)          [gated]
    #    Q = free_flow_C * Leff * (H - sill)^1.5             [free flow,
    #                                                          once a >= free_flow_ratio*(H-sill)]
    #
    #    Operating rule vs. discharge coefficient -- these are separate:
    #    - operating_rule(H, t) decides the OPENING (control logic)
    #    - discharge_coefficient/free_flow_C decide the FLOW through
    #      that opening in each regime (physical rating data)
    # -------------------------------------------------------------
    gate_1_params = {
        "H_open": 98.0,  # [m a.s.l.] level at which the gate starts to open
        "H_full": 99.5,  # [m a.s.l.] level at which the gate reaches a_max
        "a_min": 0.5,    # [m] SMALLEST opening the gate is allowed once
                         # triggered -- it jumps straight from closed (0)
                         # to a_min rather than creeping through
                         # near-zero openings, which would rest the gate
                         # seal in a partially-open position under head
                         # (accelerated seal wear / cavitation risk).
                         # EDIT to your gate's minimum-safe-opening spec.
        "a_max": 3.0,    # [m] maximum physical gate opening -- make sure
                         # this is large enough relative to your design
                         # flood head if you want the gate to reach and
                         # stay in the free-flow regime (see GatedSpillway
                         # docstring: it can revert to the gated/orifice
                         # formula at very high heads if a_max is too small)
    }
    gate_1_params.update(rule_overrides.get("gate_1", {}))
    gate_1_rule = level_trigger_rule(**gate_1_params)
    gate_1 = GatedSpillway(
        name="gate_1",
        sill_level=95.0,        # [m a.s.l.] elevation of the gate sill
        width=8.0,               # [m] width PER GATE
        operating_rule=gate_1_rule,
        discharge_coefficient=0.6 * po.get("gate_1_Cd_mult", 1.0),
                                      # constant here; replace with a
                                      # Cd(H, a) function if you have a
                                      # rating curve that varies with
                                      # head and/or opening
        free_flow_ratio=0.6,    # transition threshold a/(H-sill) -- 0.6
                                 # is a commonly used rule of thumb; EDIT
                                 # to your own gate's data if you have it
        free_flow_C=2.0 * po.get("gate_1_free_flow_C_mult", 1.0),
                                 # weir coefficient once in free flow
        free_flow_N_piers=0,    # pier contraction for the free-flow
        free_flow_Kp=0.0,       # branch, if this crest has piers
        n_gates_total=4,        # 4 identical gate bays at this spillway,
                                 # assumed OPERATED SYNCHRONOUSLY (same
                                 # opening, same time) -- if your real
                                 # spillway instead opens gates ONE AT A
                                 # TIME in a staggered sequence (common
                                 # in practice, and often given as a
                                 # gate-by-gate operating table), do NOT
                                 # use n_gates_total for that -- it can
                                 # substantially over-predict discharge
                                 # in the partially-open range. Model
                                 # each gate as its own GatedSpillway
                                 # instance instead, each with its own
                                 # staged_trigger_rule() (see outlets.py)
                                 # reading that gate's own trigger
                                 # levels off the table.
        n_gates_operational=3,  # EXAMPLE: the classic "N-1" design rule
                                 # -- assume the single most capacious
                                 # gate is unavailable. Set equal to
                                 # n_gates_total (or omit) to assume all
                                 # gates operational instead.
                                 #
                                 # NOTE: n_gates_operational treats the
                                 # unavailable gate(s) as fully removed
                                 # from the count -- i.e. ZERO discharge
                                 # from them, always. That's right for a
                                 # gate that's been physically removed/
                                 # isolated, but NOT right for a gate
                                 # that's merely STUCK CLOSED (mechanical/
                                 # electrical failure, maintenance) --
                                 # such a gate can usually still be
                                 # passively OVERTOPPED like a fixed weir
                                 # if the reservoir rises high enough over
                                 # the top of the (immobile) leaf. See the
                                 # gate_1_stuck_example below for that
                                 # pattern; Data/Corumana_117_Q5000/
                                 # case_config.py has a fully worked,
                                 # active version of it (GATES_OUT_OF_SERVICE).
    )

    # EXAMPLE (inactive by default -- not appended to the returned list
    # below) showing how to model ONE SPECIFIC gate stuck closed but
    # still overtoppable, rather than simply removed via
    # n_gates_operational. This only makes sense once gate_1 is split
    # into individual per-gate instances (see the n_gates_total comment
    # above on synchronized vs staggered gates) -- you can't single out
    # "one of the 4 bundled bays" from a single n_gates_total=4 object.
    # Uncomment and append to the returned list below (and correspondingly
    # reduce gate_1's own n_gates_total by one) to actually use this.
    #
    # gate_1_stuck_example = FreeOverflowSpillway(
    #     name="gate_1_stuck",
    #     crest_level=95.0 + 4.0,  # top of the closed gate leaf =
    #                              # sill_level + this gate's own a_max
    #     length=8.0,              # same width as one gate_1 bay
    #     C=1.86,                  # discharge coefficient for flow
    #                              # OVER a closed gate leaf -- NOT the
    #                              # same as free_flow_C (2.0 above),
    #                              # which is calibrated for a gate
    #                              # intentionally lifted clear of the
    #                              # flow, a different geometry than
    #                              # water passing over a fixed closed
    #                              # leaf. EDIT to your own gate's data.
    # )

    # -------------------------------------------------------------
    # ALTERNATIVE to the bundled gate_1 above: if you want gate-
    # availability Monte Carlo (mc_gate_availability(), near the
    # bottom of this file) to fail EACH bay independently -- rather
    # than only knowing "how many of the 4 are down," with no notion
    # of WHICH -- split gate_1 into 4 individually-named, independent
    # GatedSpillway instances instead of one bundled n_gates_total=4
    # object. Same physical/rating parameters (they're still identical
    # gates), just not synchronized into one object -- this is exactly
    # the pattern Data/Corumana_117_Q5000/case_config.py uses for its
    # six real gates. Uncomment this block AND remove gate_1 above
    # from the returned list (don't use both -- that would double-
    # count the bays) to switch a case over to this pattern; you'll
    # also need to add a `gates_out_of_service` parameter to this
    # build_outlets() function, exactly as Data/Corumana_117_Q5000/
    # case_config.py does (see its own build_outlets() docstring).
    #
    # GATE_1_BAY_NAMES = ["gate_1a", "gate_1b", "gate_1c", "gate_1d"]
    # gate_1_bays = []
    # for bay_name in GATE_1_BAY_NAMES:
    #     bay_rule = level_trigger_rule(**gate_1_params)  # same rule/
    #                                                       # params for
    #                                                       # every bay
    #                                                       # here --
    #                                                       # edit per-
    #                                                       # bay if
    #                                                       # yours are
    #                                                       # staggered
    #                                                       # (staged_
    #                                                       # trigger_
    #                                                       # rule())
    #     gate_1_bays.append(GatedSpillway(
    #         name=bay_name,
    #         sill_level=95.0, width=8.0, operating_rule=bay_rule,
    #         discharge_coefficient=0.6 * po.get("gate_1_Cd_mult", 1.0),
    #         free_flow_ratio=0.6,
    #         free_flow_C=2.0 * po.get("gate_1_free_flow_C_mult", 1.0),
    #         free_flow_N_piers=0, free_flow_Kp=0.0,
    #         n_gates_total=1,  # ONE physical bay per instance now
    #         on_off=bay_name not in (gates_out_of_service or []),
    #     ))
    # # ...then append gate_1_bays (instead of gate_1) to the returned
    # # list at the bottom of this function.

    # -------------------------------------------------------------
    # 3) Fuse gate (irreversible once triggered), WITH a pre-trip
    #    overflow -- many real fuse gates overtop (often via a
    #    labyrinth/folded crest, to pass more flow at lower head than
    #    a straight crest of the same bay width) well before they
    #    actually tip. Pre-trip and post-trip are separate regimes,
    #    each with its own crest level/length/coefficient.
    # -------------------------------------------------------------
    fuse_gate_1 = FuseGate(
        name="fuse_gate_1",
        trigger_level=101.5,        # [m a.s.l.] level at which the segment
                                    # tips/washes out (irreversibly)
        post_trip_sill=98.0,        # [m a.s.l.] sill exposed once tripped
        post_trip_length=15.0,      # [m] opening length once tripped
        post_trip_C=2.0 * po.get("fuse_gate_1_post_C_mult", 1.0),
        pre_trip_crest_level=100.5, # [m a.s.l.] this fuse gate overtops
                                    # well before it tips -- EDIT to your
                                    # structure's own crest elevation
                                    # (set equal to trigger_level, or
                                    # leave unset, if yours does NOT
                                    # overtop before tipping)
        pre_trip_length=40.0,       # [m] EFFECTIVE crest length -- often
                                    # much longer than the physical bay
                                    # width because of a labyrinth/folded
                                    # crest shape
        pre_trip_C=1.8 * po.get("fuse_gate_1_pre_C_mult", 1.0),
                                    # weir coefficient for the pre-trip
                                    # regime -- EDIT to your design data
    )

    # -------------------------------------------------------------
    # 4) Bottom outlet (submerged low-level orifice)
    #    Q = opening(H,t) * Cd(H, opening) * A * sqrt(2*g*(H - invert))
    #    Same operating-rule-vs-discharge-coefficient split as the
    #    gated spillway above. Given an operating_rule here (rather
    #    than omitting it for a fixed always-open outlet) so it also
    #    respects a minimum-opening constraint, same reasoning as
    #    gate_1 above -- see a_min there.
    #
    #    H_open/H_full are set well BELOW this case's normal operating
    #    range (reservoir_curve.csv covers H=90-105 m, H0=97 m) so this
    #    outlet reaches a_max and stays there for the levels actually
    #    seen in the Template flood -- i.e. behaves like the previous
    #    always-open default in practice, while still closing (and
    #    respecting a_min on the way back up) if the reservoir ever
    #    drops toward its minimum drawdown level. EDIT H_open/H_full to
    #    your own dam's actual minimum operating level if this outlet
    #    should behave differently.
    # -------------------------------------------------------------
    bottom_outlet_1_params = {
        "H_open": 90.0,   # [m a.s.l.]
        "H_full": 95.0,   # [m a.s.l.]
        "a_min": 0.5,     # [-] minimum opening FRACTION once triggered
                          # (same seal-protection reasoning as gate_1)
        "a_max": 1.0,     # [-] fully open
    }
    bottom_outlet_1_params.update(rule_overrides.get("bottom_outlet_1", {}))
    bottom_outlet_1_rule = level_trigger_rule(**bottom_outlet_1_params)
    bottom_outlet_1 = BottomOutlet(
        name="bottom_outlet_1",
        invert_level=80.0,   # [m a.s.l.] outlet invert elevation
        area=4.0,            # [m^2] flow area PER OUTLET
        discharge_coefficient=0.6 * po.get("bottom_outlet_1_Cd_mult", 1.0),
        n_outlets_total=2,           # 2 identical bottom outlets
        n_outlets_operational=1,    # EXAMPLE: N-1 rule -- 1 assumed
                                     # unavailable. Set equal to
                                     # n_outlets_total (or omit) to
                                     # assume all operational instead.
        operating_rule=bottom_outlet_1_rule,
    )

    return [spillway_1, gate_1, fuse_gate_1, bottom_outlet_1]


def optimizable_gates():
    """Declares which gates Module/optimize_layer2.py is allowed to tune,
    and the search bounds for each rule parameter [m a.s.l.] / [m] (or
    [-] for BottomOutlet's fractional a_min/a_max).

    Only list gates whose operating_rule is built from rule_overrides in
    build_outlets() above -- the optimizer edits nothing else (physical
    rating parameters -- crest, width, coefficients, N-1 counts -- are
    never touched).

    Every gate listed needs bounds for ALL FOUR parameters
    (H_open, H_full, a_min, a_max) -- optimize_layer2.py's decision
    vector is the same shape for every gate in one optimization run. If
    a gate's a_min is fixed rather than tunable, give it a degenerate
    bound, e.g. "a_min": (0.5, 0.5).

    Bounds should bracket a physically sensible range for your dam --
    e.g. H_open should stay within the reservoir_curve.csv range and
    below max_flood_level; a_max should not exceed the gate's real
    physical travel; a_min should not exceed the gate seal's actual
    minimum-safe-opening spec.
    """
    return {
        "gate_1": {
            "H_open": (95.5, 99.0),   # [m a.s.l.]
            "H_full": (96.5, 101.5),  # [m a.s.l.] -- must exceed H_open;
                                       # the optimizer enforces this as a
                                       # constraint, not via these bounds
            "a_min": (0.2, 1.0),       # [m] -- must be < a_max, also
                                       # enforced as a constraint
            "a_max": (0.5, 4.0),       # [m]
        },
        "bottom_outlet_1": {
            "H_open": (88.0, 93.0),   # [m a.s.l.]
            "H_full": (93.0, 96.0),   # [m a.s.l.]
            "a_min": (0.2, 0.6),       # [-]
            "a_max": (1.0, 1.0),       # [-] fixed fully-open (not tuned)
        },
    }


def mc_uncertain_params():
    """Declares which physical_overrides keys Module/mc_layer3.py's
    INNER loop is allowed to sample, and their distributions --
    mirrors optimizable_gates()'s role for Layer 2, but for Layer 3's
    physical/rating uncertainty instead of the operating rule (see
    build_outlets()'s own physical_overrides docstring above for the
    full list of keys this case's build_outlets() understands).

    PLACEHOLDER distributions -- illustrative starting points only
    (typical orders of magnitude for weir/gate discharge-coefficient
    calibration uncertainty), NOT calibrated to any real structure.
    EDIT with real values (rating-curve calibration data, if you have
    it) before trusting Layer 3 results for any real decision.

    dist="lognormal_cv": multiplicative, median 1.0, given coefficient
        of variation (e.g. cv=0.05 means the multiplier typically falls
        within roughly +/-5% of 1.0).
    dist="normal": additive, given (mean, sigma) in the target
        parameter's own units.

    A case with NO mc_uncertain_params() at all still runs under Layer
    3 -- H0 and the reservoir stage-storage curve are varied regardless
    (see mc_layer3.py's default_inner_sampler) -- it just skips
    discharge-coefficient perturbation, which Layer 3 will note rather
    than error on.
    """
    return {
        "spillway_1_C_mult": {"dist": "lognormal_cv", "cv": 0.05},
        "gate_1_Cd_mult": {"dist": "lognormal_cv", "cv": 0.05},
        "gate_1_free_flow_C_mult": {"dist": "lognormal_cv", "cv": 0.05},
        "fuse_gate_1_pre_C_mult": {"dist": "lognormal_cv", "cv": 0.05},
        "fuse_gate_1_post_C_mult": {"dist": "lognormal_cv", "cv": 0.05},
        "bottom_outlet_1_Cd_mult": {"dist": "lognormal_cv", "cv": 0.05},
    }


# mc_outer_distribution() -- NOT defined for this template, deliberately.
#
# Module/mc_layer3.py's OUTER loop (climate/flood-frequency uncertainty)
# can use a case's own real flood-frequency and design-hydrograph data
# instead of a generic placeholder spread, if the case provides it via
# an mc_outer_distribution() function returning a dict like:
#
#   def mc_outer_distribution():
#       case_dir = os.path.dirname(os.path.abspath(__file__))
#       return {
#           "target_return_period": 5000,  # which flood this case represents
#           "curve_csv": os.path.join(case_dir, "flood_frequency_curve.csv"),
#           "alt_studies_csv": os.path.join(case_dir, "flood_frequency_alt_studies.csv"),
#           "volume_duration_csv": os.path.join(case_dir, "flood_duration_volume_table.csv"),
#       }
#
# where:
#   flood_frequency_curve.csv       : return_period_years,peak_Q_m3s[,source]
#       -- your adopted design flood-frequency curve (statistical fit
#       or otherwise). Sets the outer loop's PEAK scale factor's median.
#   flood_frequency_alt_studies.csv : return_period_years,study,peak_Q_m3s
#       -- independent past estimates of the same flood(s), if you have
#       more than one hydrology study for this structure. The spread
#       between them and your adopted curve sets the PEAK scale
#       factor's uncertainty (CV) -- a real, if rough, data-grounded
#       number instead of an assumed one. Needs 2+ rows to be used;
#       otherwise Layer 3 falls back to a generic placeholder CV.
#   volume_duration_csv : return_period_years,storm_duration_hours,peak_Q_m3s,volume_Mm3
#       -- your adopted design hydrographs' PAIRED peak and volume by
#       return period. Sets the outer loop's VOLUME scale factor
#       (median AND uncertainty), sampled INDEPENDENTLY of the peak
#       scale factor -- flood volume, not just peak, is often a
#       first-order control on the routed reservoir level for any
#       reservoir that provides meaningful attenuation, and assuming
#       peak and volume always move together (the alternative, if this
#       hook is absent) can hide the routing-sensitive scenarios that
#       matter most.
#
# Without this hook, Layer 3 still runs correctly -- it just samples
# peak and volume scale factors from a generic placeholder log-normal
# spread (median 1.0, a flat assumed CV) around this case's own
# inflow_hydrograph.csv, and prints a note saying so. Add the three
# CSVs above (real data for your own case) and this function once you
# have them.


# mc_gate_availability() -- NOT defined for this template, deliberately
# (same opt-in reasoning as mc_outer_distribution() above).
#
# Module/mc_layer3.py's inner loop can also model gate-failure-to-open
# risk, IF a case declares which gates are eligible to fail and a
# per-gate failure probability:
#
#   def mc_gate_availability():
#       return {
#           "gates": ["gate_1a", "gate_1b", "gate_1c", "gate_1d"],
#           "p_fail": 0.05,   # SAME probability applied to EVERY gate
#       }
#
# Each gate's outcome is drawn INDEPENDENTLY every inner draw (the
# standard binomial screening model P(K=k) = C(n,k)*p^k*(1-p)^(n-k)
# for the number of gates K that fail to open) -- NOT correlated/
# common-cause failure, which this simple model does not attempt to
# capture (a single maintenance lapse taking out several gates AT ONCE
# would need a different, correlated sampling mechanism, not this one).
#
# This only makes sense for a case using INDIVIDUALLY NAMED gates (see
# the "ALTERNATIVE to the bundled gate_1" block above) -- a case using
# the bundled n_gates_total form has no individual gate identities for
# this hook to refer to, and build_outlets() also needs its own
# `gates_out_of_service` parameter added (see Data/Corumana_117_Q5000/
# case_config.py for the full worked pattern, including how a gate
# that's merely stuck closed -- as opposed to fully removed -- can
# still be modeled as passively overtoppable).
#
# Without this hook, Layer 3 assumes every gate stays operational for
# every draw (no failure-to-open risk modeled), and prints a note
# saying so.


def scalars(case_dir: str):
    """Reads scalars.csv into a dict of {name: value}.
    case_dir is this case's own Data/<CaseName>/ folder, passed in by
    run_case.py."""
    import csv
    vals = {}
    scalars_csv = os.path.join(case_dir, "scalars.csv")
    with open(scalars_csv, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            vals[row["name"]] = float(row["value"])
    return vals
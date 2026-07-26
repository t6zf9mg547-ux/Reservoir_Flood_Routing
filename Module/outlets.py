"""
outlets.py
----------
This is the ONE file in the whole tool where you define the physical
discharge law of each outlet works at your dam, and (later, for PLC-type
rules) the gate operating logic. Everything else (solver, reservoir
curve, hydrograph reading) is generic and you should not need to touch it.

Each outlet is a small Python object with a `discharge(H, t) -> Q [m3/s]`
method: the outflow depends only on reservoir level H and time t. You
add, remove, or edit outlets by editing Data/<CaseName>/case_config.py,
without touching this file, UNLESS you need a genuinely new discharge
law (e.g. a different weir shape), in which case you add a new class
here following the same pattern.

IMPORTANT modeling note on statefulness (fuse gates):
    RK4 evaluates the derivative dH/dt at several trial (H, t) points
    within a single timestep (k1..k4), not just at the accepted start
    and end of the step. A fuse gate's "has it tipped yet?" state must
    NOT be updated from these trial evaluations, or the trigger could
    flip spuriously depending on which RK4 stage happens to overshoot
    the trigger level. Instead, the solver calls `outlet.commit(H)`
    exactly ONCE per accepted timestep (after the RK4 step is finalized),
    and outlets that carry state (like FuseGate) update themselves only
    inside `commit`. Stateless outlets simply ignore `commit`.
"""

from __future__ import annotations
from typing import Callable, Optional


class Outlet:
    """Base class. Subclass and implement discharge(H, t)."""

    name: str = "outlet"

    def discharge(self, H: float, t: float) -> float:
        raise NotImplementedError

    def commit(self, H: float) -> None:
        """Called once per accepted RK4 step with the accepted level H.
        Stateless outlets can ignore this (default no-op)."""
        pass

    def reset(self) -> None:
        """Called before a simulation run to reset any internal state
        (so the same outlet object can be reused across runs, e.g. in
        a Monte Carlo loop later)."""
        pass


class FreeOverflowSpillway(Outlet):
    """Ungated broad-crested weir, with optional pier contraction effect:

        Q = C * Leff * (H - crest)^1.5      for H > crest, else 0
        Leff = L - 2 * N * Kp * (H - crest)  (effective length, reduced
                                               by flow contraction at piers)

    Parameters
    ----------
    crest_level : sill/crest elevation [m a.s.l.]
    length      : GROSS crest length L [m] (i.e. the physical total
                  width of the spillway opening, before any reduction
                  for piers)
    C           : weir discharge coefficient [m^0.5/s], typically ~1.7-2.2
                  for a broad-crested weir (adjust to match model calibration)
    N_piers     : number of piers spanning the spillway crest (0 = none,
                  the default -- Leff then just equals L)
    Kp          : pier contraction coefficient [-], user-defined. Reduces
                  the effective crest length as head rises, since flow
                  contracts (and separates) around each pier nose --
                  typical values run roughly 0.01 (well-rounded pier
                  nose) to ~0.02-0.03 (square-nosed pier), but this
                  should be set from your own hydraulic design/model
                  data or reference tables (e.g. USBR), not assumed.
    on_off      : set False to take this outlet out of service without
                  removing it from the case (e.g. to model an outlet
                  closed for maintenance, or to test a "what if this
                  spillway is unavailable" scenario)

    Leff is clipped to >= 0: if (H - crest) grows large enough relative
    to pier spacing that this formula would predict a negative or
    vanishing effective length, that's a sign the head has left the
    range this simple linear contraction model is valid for (piers
    fully drowned/very high relative head) -- treat results with
    caution near that limit and consult a proper rating curve if it
    matters for your case.
    """

    def __init__(self, name: str, crest_level: float, length: float,
                 C: float = 2.0, N_piers: int = 0, Kp: float = 0.0,
                 on_off: bool = True):
        self.name = name
        self.crest_level = crest_level
        self.length = length
        self.C = C
        self.N_piers = N_piers
        self.Kp = Kp
        self.on_off = on_off

    def discharge(self, H: float, t: float) -> float:
        if not self.on_off or H <= self.crest_level:
            return 0.0
        head = H - self.crest_level
        L_eff = self.length - 2 * self.N_piers * self.Kp * head
        L_eff = max(L_eff, 0.0)
        return self.C * L_eff * head ** 1.5


def gated_spillway_unit_discharge(H: float, a: float, sill_level: float, width: float,
                                    discharge_coefficient=0.6, free_flow_ratio: float = 0.6,
                                    free_flow_C: float = 2.0, free_flow_N_piers: int = 0,
                                    free_flow_Kp: float = 0.0) -> float:
    """Discharge through ONE gate bay given an EXPLICIT opening `a`
    (not via an operating_rule) -- the same two-regime law
    GatedSpillway.discharge() uses, factored out into a standalone
    function so it has a single source of truth AND can be called
    directly wherever an opening needs to be SOLVED FOR rather than
    looked up (see sequential_fill_rule() below, which inverts this to
    find the `a` that produces a target discharge). See GatedSpillway's
    own docstring for the physics behind the two regimes."""
    G = 9.81
    if H <= sill_level or a <= 0.0:
        return 0.0
    head = H - sill_level
    if a >= free_flow_ratio * head:
        L_eff = max(width - 2 * free_flow_N_piers * free_flow_Kp * head, 0.0)
        return free_flow_C * L_eff * head ** 1.5
    Cd = discharge_coefficient(H, a) if callable(discharge_coefficient) else discharge_coefficient
    return Cd * width * a * (2 * G * head) ** 0.5


class GatedSpillway(Outlet):
    """Gate-controlled outlet using the standard underflow/orifice gate
    equation WHILE THE GATE CONSTRAINS THE FLOW:

        Q = Cd(H, a) * width * a * sqrt(2 * g * H_eff)

    ONCE THE GATE OPENING BECOMES LARGE RELATIVE TO THE HEAD, the gate
    lifts clear of the flow path and no longer constrains it -- the
    spillway crest itself becomes the control, and discharge instead
    follows the free-overflow weir equation (see FreeOverflowSpillway,
    including its optional pier contraction correction):

        Q = free_flow_C * Leff * (H - sill)^1.5
        Leff = width - 2*free_flow_N_piers*free_flow_Kp*(H - sill)

    The transition is triggered once the gate opening a exceeds
    `free_flow_ratio * (H - sill)`. 0.6 is a commonly used rule of
    thumb for this transition -- adjust `free_flow_ratio` to match your
    own gate/spillway's actual transition point if you have hydraulic
    model or manufacturer data. Note this creates a step change in the
    discharge law at the transition (the orifice and free-flow formulas
    won't generally agree exactly at the boundary) -- that's inherent
    to modeling this as a hard regime switch, matching how it's usually
    presented in hydraulic references; if a smoother blend matters for
    your case, wrap two GatedSpillway/FreeOverflowSpillway instances in
    a SwitchedOutlet with a narrow ramp instead of a hard cutover.

    THIS IS TWO SEPARATE THINGS, NOT ONE -- keep them that way:
        1. `operating_rule(H, t)` decides WHAT THE GATE OPENING SHOULD
           BE right now, given the reservoir level (and later, a
           forecast). This is control/PLC logic -- see
           `level_trigger_rule()` below for the standard level-triggered
           version, and it's what an optimizer would tune in a later
           phase.
        2. `discharge_coefficient` (constant, or Cd(H, a)) and
           `free_flow_C` decide HOW MUCH WATER ACTUALLY FLOWS in each
           regime -- this is physical calibration/rating data for the
           gate and crest themselves.
    Keeping these separate means you can change the control logic
    without touching either regime's physical rating curve, and vice
    versa.

    Use `level_trigger_rule(H_open, H_full, a_max=<max gate opening>,
    a_min=<min allowed opening once triggered>)` as the operating_rule
    to get a gate that JUMPS from fully closed to a_min the instant the
    reservoir rises past H_open, then ramps from a_min to a_max meters
    open as it rises further to H_full -- rather than creeping through
    a continuous range of near-zero openings, which many real gate
    seals shouldn't be left resting in under head. Omit a_min (defaults
    to 0.0) if your gate has no such constraint.

    Modeling several IDENTICAL gates together: if your spillway has
    several identical, interchangeable gate bays operated THE SAME WAY
    AT THE SAME TIME, set `n_gates_total` (how many exist) and
    `n_gates_operational` (how many are assumed available right now --
    defaults to n_gates_total, i.e. all available). Total discharge is
    simply the per-gate discharge times `n_gates_operational`.

    CAUTION -- "operated the same way at the same time" is a real
    assumption, not just a modeling convenience: many real multi-gate
    spillways instead open gates ONE AT A TIME in a staggered sequence
    (each gate starting a bit later than the last, in level), not
    synchronously. Using n_gates_total/n_gates_operational for a
    spillway that's actually staggered can substantially OVER-predict
    discharge in the partially-open range -- non-conservative for a
    flood study, since it understates the peak reservoir level a real,
    slower-opening spillway would allow. If you have an actual
    gate-by-gate operating table (a common way these are specified --
    "at this level, open gate N to this much"), model it as ONE
    GatedSpillway instance PER GATE (same physical/rating parameters,
    since the gates are still identical -- just not synchronized),
    each with its own `staged_trigger_rule()` (below) reading that
    gate's own trigger levels off the table, rather than one shared
    instance with n_gates_total. Cross-check the result's total
    discharge against the source table across the full head range
    before trusting it.

    This is exactly what you need for the common "N-1" design/safety
    rule used in many countries -- assume the single most capacious
    gate is out of service, and verify the rest still pass the design
    flood -- by setting n_gates_operational = n_gates_total - 1.

    If your gates are NOT identical (e.g. different widths/capacities),
    don't use n_gates_total/n_gates_operational -- instead, define one
    GatedSpillway instance per physical gate, and set on_off=False on
    whichever specific gate you want to take out of service for the
    N-1 check.
    """

    G = 9.81

    def __init__(self, name: str, sill_level: float, width: float,
                 operating_rule: Callable[[float, float], float],
                 discharge_coefficient=0.6,
                 free_flow_ratio: float = 0.6, free_flow_C: float = 2.0,
                 free_flow_N_piers: int = 0, free_flow_Kp: float = 0.0,
                 n_gates_total: int = 1, n_gates_operational: Optional[int] = None,
                 on_off: bool = True):
        self.name = name
        self.sill_level = sill_level
        self.width = width
        self.operating_rule = operating_rule
        self.discharge_coefficient = discharge_coefficient
        self.free_flow_ratio = free_flow_ratio
        self.on_off = on_off
        self.n_gates_total = n_gates_total
        self.n_gates_operational = (
            n_gates_operational if n_gates_operational is not None else n_gates_total
        )
        if not 0 <= self.n_gates_operational <= self.n_gates_total:
            raise ValueError("n_gates_operational must be between 0 and n_gates_total")
        # Free-flow regime reuses FreeOverflowSpillway (same sill level
        # and width-as-crest-length, PER GATE), so pier contraction is
        # available here too if the crest has piers.
        self._free_flow_outlet = FreeOverflowSpillway(
            name=f"{name}_free_flow_branch",
            crest_level=sill_level,
            length=width,
            C=free_flow_C,
            N_piers=free_flow_N_piers,
            Kp=free_flow_Kp,
        )

    def _Cd(self, H: float, a: float) -> float:
        if callable(self.discharge_coefficient):
            return self.discharge_coefficient(H, a)
        return self.discharge_coefficient

    def discharge(self, H: float, t: float) -> float:
        if not self.on_off or self.n_gates_operational <= 0:
            return 0.0
        a = max(self.operating_rule(H, t), 0.0)
        q_per_gate = gated_spillway_unit_discharge(
            H, a, self.sill_level, self.width, self.discharge_coefficient,
            self.free_flow_ratio, self._free_flow_outlet.C,
            self._free_flow_outlet.N_piers, self._free_flow_outlet.Kp,
        )
        return self.n_gates_operational * q_per_gate


class BottomOutlet(Outlet):
    """Bottom outlet / low-level orifice: submerged orifice equation
        Q = opening(H,t) * Cd(H, opening) * A * sqrt(2 * g * (H - invert_level))
    valid while the outlet is submerged (H > invert_level); returns 0
    otherwise. If invert_level is below the outlet's own centerline in
    reality, treat invert_level as the effective driving-head reference
    elevation for your rating curve.

    `discharge_coefficient` follows the same pattern as GatedSpillway:
    pass a constant float, or a callable Cd(H, opening) -> float if your
    rating data shows the coefficient varying with head and/or opening.

    operating_rule(H, t) -> opening fraction in [0, 1] (1.0 = fully
    open). Defaults to always fully open if none is given -- pass
    `level_trigger_rule(...)` (with the default a_max=1.0, and a_min
    set if this outlet's gate/valve also shouldn't be left resting at
    a near-zero opening) if this outlet is itself gated rather than a
    fixed always-open orifice.

    Modeling several IDENTICAL bottom outlets together: if your dam has
    several identical, interchangeable bottom outlets operated the same
    way, set `n_outlets_total` (how many exist) and
    `n_outlets_operational` (how many are assumed available right now
    -- defaults to n_outlets_total, i.e. all available). Total
    discharge is simply the per-outlet discharge times
    `n_outlets_operational`.

    This covers the common "N-1" design/safety rule used in many
    countries -- assume the single most capacious outlet is out of
    service, and verify the rest still pass the design flood -- by
    setting n_outlets_operational = n_outlets_total - 1.

    If your outlets are NOT identical, don't use
    n_outlets_total/n_outlets_operational -- instead, define one
    BottomOutlet instance per physical outlet, and set on_off=False on
    whichever specific one you want to take out of service.
    """

    G = 9.81

    def __init__(self, name: str, invert_level: float, area: float,
                 operating_rule: Optional[Callable[[float, float], float]] = None,
                 discharge_coefficient=0.6,
                 n_outlets_total: int = 1, n_outlets_operational: Optional[int] = None,
                 on_off: bool = True):
        self.name = name
        self.invert_level = invert_level
        self.area = area
        self.discharge_coefficient = discharge_coefficient
        self.on_off = on_off
        self.n_outlets_total = n_outlets_total
        self.n_outlets_operational = (
            n_outlets_operational if n_outlets_operational is not None else n_outlets_total
        )
        if not 0 <= self.n_outlets_operational <= self.n_outlets_total:
            raise ValueError("n_outlets_operational must be between 0 and n_outlets_total")
        # default: always fully open if no rule supplied
        self.operating_rule = operating_rule or (lambda H, t: 1.0)

    def _Cd(self, H: float, opening: float) -> float:
        if callable(self.discharge_coefficient):
            return self.discharge_coefficient(H, opening)
        return self.discharge_coefficient

    def discharge(self, H: float, t: float) -> float:
        if not self.on_off or H <= self.invert_level or self.n_outlets_operational <= 0:
            return 0.0
        opening = min(max(self.operating_rule(H, t), 0.0), 1.0)
        if opening <= 0.0:
            return 0.0
        head = H - self.invert_level
        Cd = self._Cd(H, opening)
        q_per_outlet = opening * Cd * self.area * (2 * self.G * head) ** 0.5
        return self.n_outlets_operational * q_per_outlet


class FuseGate(Outlet):
    """Fuse gate (or fusible levee/dike segment): a structure that can
    OVERTOP AND PASS FLOW BEFORE IT TIPS, then tips/washes out
    irreversibly once the reservoir level exceeds trigger_level, after
    which it behaves as a free overflow weir at a lower (post-trip)
    sill for the remainder of the simulation -- even if the level later
    falls back down.

    Two entirely separate discharge regimes, each with its own weir
    geometry, because they're usually physically very different
    structures:

    1. PRE-TRIP (before tipping): the fuse gate's own crest, which many
       real fuse gates overtop well before failing -- often a labyrinth
       (folded/zigzag) shape specifically to pass more flow at a lower
       head than a straight crest of the same physical bay width would.
       This is why `pre_trip_length` is a SEPARATE parameter from
       `post_trip_length`, not the same number -- a labyrinth crest's
       effective length is typically much longer than the straight-line
       bay width, and its coefficient can differ too:
           Q = pre_trip_C * pre_trip_length * (H - pre_trip_crest_level)^1.5
       Set pre_trip_length=0 (the default) if your fuse gate does NOT
       overtop before tipping (a simple non-overtopping plug) -- this
       recovers the original all-or-nothing behavior.

    2. POST-TRIP (after tipping): the structure has washed out/rotated
       clear, exposing a (usually much larger, plain) opening:
           Q = post_trip_C * post_trip_length * (H - post_trip_sill)^1.5

    State (tripped or not) is only updated in `commit()`, once per
    accepted RK4 step -- see the module docstring for why this matters
    with RK4 -- so the trigger can't fire spuriously from an
    intermediate trial sub-step.
    """

    def __init__(self, name: str, trigger_level: float,
                 post_trip_sill: float, post_trip_length: float, post_trip_C: float = 2.0,
                 pre_trip_crest_level: Optional[float] = None,
                 pre_trip_length: float = 0.0, pre_trip_C: float = 2.0):
        self.name = name
        self.trigger_level = trigger_level
        self.post_trip_sill = post_trip_sill
        self.post_trip_length = post_trip_length
        self.post_trip_C = post_trip_C
        # Defaults to trigger_level (with pre_trip_length=0, i.e. no
        # pre-trip flow) if not given -- recovers the original
        # all-or-nothing fuse gate behavior.
        self.pre_trip_crest_level = (
            pre_trip_crest_level if pre_trip_crest_level is not None else trigger_level
        )
        self.pre_trip_length = pre_trip_length
        self.pre_trip_C = pre_trip_C
        self.tripped = False

    def discharge(self, H: float, t: float) -> float:
        if self.tripped:
            if H <= self.post_trip_sill:
                return 0.0
            return self.post_trip_C * self.post_trip_length * (H - self.post_trip_sill) ** 1.5
        if self.pre_trip_length <= 0.0 or H <= self.pre_trip_crest_level:
            return 0.0
        return self.pre_trip_C * self.pre_trip_length * (H - self.pre_trip_crest_level) ** 1.5

    def commit(self, H: float) -> None:
        if not self.tripped and H > self.trigger_level:
            self.tripped = True

    def reset(self) -> None:
        self.tripped = False


class ConstantOutlet(Outlet):
    """Discharges a fixed value regardless of H or t (e.g. Q=0, or a
    capped/saturated orifice flow). Mainly useful as a building block
    inside SwitchedOutlet, to represent a near-step onset that a smooth
    weir/orifice formula can't reproduce on its own."""

    def __init__(self, name: str, value: float = 0.0):
        self.name = name
        self.value = value

    def discharge(self, H: float, t: float) -> float:
        return self.value


class SwitchedOutlet(Outlet):
    """Generic pattern for a gate/outlet whose discharge law itself
    changes character between two states (e.g. a small partial-opening
    weir law before the gate lifts, and a different, larger-crest-length
    weir law once it is fully raised) -- as opposed to GatedSpillway,
    which scales ONE law by an opening fraction.

    discharge(H, t) = (1 - w) * outlet_low.discharge(H, t)
                      +     w * outlet_high.discharge(H, t)
    where w = operating_rule(H, t) in [0, 1] (use level_trigger_rule for
    the standard PLC-style version; a narrow (H_open, H_full) band gives
    a near-discrete switch, a wide one gives a smooth blend).

    commit()/reset() propagate to both child outlets, so a FuseGate (or
    any stateful outlet) can be used as either branch.
    """

    def __init__(self, name: str, outlet_low: Outlet, outlet_high: Outlet,
                 operating_rule: Callable[[float, float], float]):
        self.name = name
        self.outlet_low = outlet_low
        self.outlet_high = outlet_high
        self.operating_rule = operating_rule

    def discharge(self, H: float, t: float) -> float:
        w = min(max(self.operating_rule(H, t), 0.0), 1.0)
        return (1 - w) * self.outlet_low.discharge(H, t) + w * self.outlet_high.discharge(H, t)

    def commit(self, H: float) -> None:
        self.outlet_low.commit(H)
        self.outlet_high.commit(H)

    def reset(self) -> None:
        self.outlet_low.reset()
        self.outlet_high.reset()


# ---------------------------------------------------------------------
# Reusable operating-rule builder for GatedSpillway / BottomOutlet
# ---------------------------------------------------------------------

def level_trigger_rule(H_open: float, H_full: float,
                        a_max: float = 1.0,
                        a_min: float = 0.0) -> Callable[[float, float], float]:
    """Standard PLC-style level-triggered gate rule:
        - returns 0 (fully closed) at/below H_open
        - JUMPS to a_min the instant H exceeds H_open (not a continuous
          ramp from 0) -- see a_min below for why
        - ramps linearly from a_min to a_max between H_open and H_full
        - returns a_max at/above H_full

    a_min : smallest allowed opening once the gate is triggered [same
        units as a_max: a physical opening in meters for GatedSpillway/
        BottomOutlet, or a 0-1 fraction for BottomOutlet's default
        interpretation / SwitchedOutlet's blend weight]. Many real
        gates must never be left barely cracked open -- resting on a
        near-closed seal under high head accelerates seal wear and can
        cause cavitation/vibration -- so the rule is either FULLY
        CLOSED (0) or open AT LEAST a_min, never something in between.
        Default 0.0 recovers the original continuous-ramp-from-zero
        behavior for gates without this constraint.

    Use the default a_max=1.0 when the outlet interprets its rule's
    output as a 0-1 opening FRACTION (BottomOutlet, SwitchedOutlet's
    blend weight). Set a_max=<max physical gate opening in meters> to
    have this directly return a physical opening for GatedSpillway.

    This is the deterministic, forecast-independent rule matching how a
    real PLC executes a fixed operating table. It ignores t entirely,
    but keeps the (H, t) signature so it's interchangeable with a future
    forecast-aware rule without changing the outlet classes.
    """
    if H_full <= H_open:
        raise ValueError("H_full must be greater than H_open")
    if not 0.0 <= a_min < a_max:
        raise ValueError("a_min must be >= 0 and strictly less than a_max")

    def rule(H: float, t: float) -> float:
        if H <= H_open:
            return 0.0
        if H >= H_full:
            return a_max
        frac = (H - H_open) / (H_full - H_open)
        return a_min + frac * (a_max - a_min)

    return rule


def staged_trigger_rule(stages: list[tuple[float, float]]
                         ) -> Callable[[float, float], float]:
    """Discrete-step gate rule for cases where the real operating
    procedure is specified as a fixed TABLE of (level threshold,
    opening) steps -- e.g. "at 117.5 m, set this gate to 0.5 m; at
    118.1 m, to 2 m; ..." -- rather than a continuous ramp. This is how
    many real dam operation manuals / PLC lookup tables are actually
    written (especially for staggered multi-gate sequences where each
    gate follows the same step sequence but offset in level from the
    others -- model each gate as a separate GatedSpillway/BottomOutlet
    instance, each with its own staged_trigger_rule using that gate's
    own trigger levels).

    stages : list of (H_threshold, opening) pairs; needn't be
        pre-sorted. Returns the opening of the last stage whose
        H_threshold <= H (0.0 if H is below every threshold). This is
        an exact reproduction of a step-function schedule -- if your
        real rule is closer to a continuous ramp instead, use
        level_trigger_rule() (with its a_min) rather than approximating
        a ramp with many small steps here.

    Example: reproducing "opens to 0.5 m at 117.5, to 2 m at 118.1, to
    4 m at 118.7, fully open (7 m) at 119.3":
        staged_trigger_rule([(117.5, 0.5), (118.1, 2.0),
                              (118.7, 4.0), (119.3, 7.0)])
    """
    sorted_stages = sorted(stages, key=lambda pair: pair[0])

    def rule(H: float, t: float) -> float:
        opening = 0.0
        for H_threshold, a in sorted_stages:
            if H >= H_threshold:
                opening = a
            else:
                break
        return opening

    return rule


def sequential_fill_rule(gate_index: int, n_gates: int,
                          target_release_fn: Callable[[float, float], float],
                          sill_level: float, width: float,
                          a_min: float, a_max: float,
                          discharge_coefficient=0.6, free_flow_ratio: float = 0.6,
                          free_flow_C: float = 2.0, free_flow_N_piers: int = 0,
                          free_flow_Kp: float = 0.0,
                          bisect_iters: int = 60) -> Callable[[float, float], float]:
    """Builds the operating_rule for gate number `gate_index` (0-based,
    in ACTIVATION ORDER -- gate 0 opens first) out of `n_gates`
    IDENTICAL gate bays, coordinated so that AT MOST ONE gate is ever
    partway open at a time: every gate before the active one is fully
    open (a_max), every gate after it is fully closed (0), and the
    active gate is set to whatever opening makes TOTAL combined
    discharge across all n_gates match `target_release_fn(H, t)` as
    closely as possible WITHOUT EXCEEDING it.

    This is for real operating philosophies specified as "bring gates
    into service one at a time to track a target release" (e.g. "hold
    the reservoir at a target level by releasing no more than what's
    coming in") rather than a fixed level-triggered table -- see
    level_trigger_rule()/staged_trigger_rule() for that simpler case.
    Since it needs a TARGET RELEASE (not just the reservoir level) to
    decide each gate's opening, `target_release_fn(H, t)` is passed in
    as a callable rather than baked into a fixed threshold -- typically
    this reads an inflow Hydrograph's `.discharge(t)`, gated by a level
    check (e.g. `lambda H, t: inflow.discharge(t) if H > FSL else 0.0`
    for a "pass through the flood once above the target level" rule).

    Mechanics per gate (each of the n_gates operating_rule callables
    built from this function repeats this SAME calculation
    independently and just returns its own gate_index's result --
    stateless, so it's safe to call at every RK4 trial sub-step, same
    as any other operating_rule):
        remaining = target_release_fn(H, t)
        for each gate 0..n_gates-1, in order:
            q_min = discharge at a_min ; q_max = discharge at a_max
            if remaining <= 0 or q_min > remaining:
                this gate stays CLOSED (opening it to even a_min would
                either not be needed, or would overshoot the target --
                per the "never release more than the target" rule,
                staying closed is the safe choice here)
            elif remaining >= q_max:
                this gate goes FULLY OPEN; subtract q_max from
                remaining; move on to the next gate with the leftover
            else:
                bisect for the exact `a` in [a_min, a_max] that hits
                `remaining` for THIS gate; every later gate then stays
                closclosed since remaining is now fully consumed

    Physical discharge is assumed MONOTONIC in opening `a` (more open
    -> more or equal flow) for the bisection to be reliable -- true in
    general, but note GatedSpillway's own docstring on the gated/
    free-flow transition creating a STEP CHANGE in the discharge law:
    once a gate clears the flow path, further opening doesn't add any
    more discharge (the free-flow formula depends only on H, not on
    a), so there's typically a JUMP at the transition -- discharge can
    leap from some value Q_transition_gated straight to a HIGHER
    Q_transition_free_flow with NO achievable opening in between. If a
    target release falls in that gap, this correctly settles at the
    top of the gated regime (the highest value that does NOT exceed
    the target) rather than jumping into free flow and overshooting --
    exactly the intended "never exceed the target" behavior, just
    worth knowing why a gate sometimes stops just short of the
    transition rather than continuing further open for no extra flow.
    """
    def solve_all(H: float, t: float) -> list[float]:
        remaining = max(target_release_fn(H, t), 0.0)
        openings = [0.0] * n_gates
        for i in range(n_gates):
            q_min = gated_spillway_unit_discharge(
                H, a_min, sill_level, width, discharge_coefficient,
                free_flow_ratio, free_flow_C, free_flow_N_piers, free_flow_Kp)
            q_max = gated_spillway_unit_discharge(
                H, a_max, sill_level, width, discharge_coefficient,
                free_flow_ratio, free_flow_C, free_flow_N_piers, free_flow_Kp)
            if remaining <= 0.0 or q_min > remaining:
                break  # this and every later gate stay closed (0.0)
            elif remaining >= q_max:
                openings[i] = a_max
                remaining -= q_max
            else:
                lo, hi = a_min, a_max
                for _ in range(bisect_iters):
                    mid = 0.5 * (lo + hi)
                    q_mid = gated_spillway_unit_discharge(
                        H, mid, sill_level, width, discharge_coefficient,
                        free_flow_ratio, free_flow_C, free_flow_N_piers, free_flow_Kp)
                    if q_mid < remaining:
                        lo = mid
                    else:
                        hi = mid
                # Return `lo`, NOT the midpoint (lo+hi)/2 -- the loop
                # invariant guarantees discharge(lo) < remaining held at
                # every step, so `lo` never exceeds the target. `hi` has
                # no such guarantee: at a discontinuity (see docstring),
                # lo and hi converge to opposite SIDES of the jump, and
                # the midpoint can land on either side depending on
                # floating-point rounding -- i.e. it can silently
                # OVERSHOOT the target right when precision matters
                # most. `lo` is always safe, in both the smooth case
                # (converges to the true root, same as the midpoint
                # would, to within float tolerance) and the
                # discontinuous case (stays on the non-exceeding side).
                openings[i] = lo
                break  # exact match found; every later gate stays closed
        return openings

    def rule(H: float, t: float) -> float:
        return solve_all(H, t)[gate_index]

    return rule
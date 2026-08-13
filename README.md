# Reservoir_Flood_Routing

Python reservoir flood-routing kernel: RK4 solution of the reservoir
continuity equation

    dH/dt = (Qin(t) - Qout(H,t)) / B(H)

for a single dam with multiple outlet works (ungated spillways, gated
spillways, fuse gates, bottom outlets).

## How to run

```
python Module/run_case.py
```

This lists the case folders found under `Data/` (e.g. `Data/Template/`)
and asks you to pick one. It then writes `Output/<CaseName>/Baseline/results.csv`
and `Plot/<CaseName>/Baseline/routing_result.png` -- the same case name, mirrored
under both.

To skip the prompt (useful for scripted/batch runs later):
```
python Module/run_case.py Template
```

To compare a run against a reference results file (optional):
```
python Module/validate_against_original.py Template
```
(edit `ORIGINAL_CSV` at the top of that file to point at your reference).

## Folder structure & how to enter data

- **Data/<CaseName>/** -- one subfolder per case, fully self-contained,
  e.g. `Data/Template/` contains:
  - `reservoir_curve.csv`  : H_masl, Surface_km2 (stage-surface curve --
    area in km^2, converted to m^2 internally)
  - `inflow_hydrograph.csv`: time_h, Q_m3s (flood hydrograph)
  - `withdrawal_schedule.csv`: time_h, Q_m3s (withdrawal -- e.g. turbine
    offtake, irrigation, compensation flow, or any other outflow that
    doesn't depend on reservoir level)
  - `scalars.csv`: t_max, dt, print_every, H0, and optionally
    dam_crest_level, max_flood_level (see "Safety thresholds" below).
    Also where Layer 3's mandatory `mc_outer_peak_source`/
    `mc_outer_volume_source`/`mc_gate_reliability_source` selections,
    the optional `mc_peak_volume_dependence`/`mc_peak_volume_tau`
    pair, and the optional `mc_outer_hydrograph_source` switch live
    (see "Layer 3" below) -- required only if you run Layer 3 for this
    case, irrelevant to Layer 1/2
  - `case_config.py` : the ONE file where you define which outlets
    exist and their parameters/operating rules (see below)

  To create a new scenario: duplicate the whole `Data/Template/` folder
  and rename the duplicate (e.g. `Data/MyCase/`), then edit the CSVs
  and `case_config.py` inside it. That file's name never changes --
  it's always `case_config.py` -- `run_case.py` finds it by which
  folder it's in, not by matching filenames, so there's no renaming of
  code involved, ever; only the folder gets renamed.

  Note: a `__pycache__` folder will reappear inside each case folder
  after you run it (Python creates this automatically next to any
  `.py` file it imports). It's harmless and safe to delete anytime --
  this is the trade-off for keeping case setup to a single
  duplicate-and-rename step, instead of splitting data and code across
  two folders.

- **Module/** -- the generic engine (`run_case.py`, `reservoir.py`,
  `hydrograph.py`, `outlets.py`, `solver.py`, `validate_against_original.py`).
  You shouldn't need to edit these unless you need a genuinely new type
  of outlet discharge law (add a new class in `outlets.py` following
  the existing pattern).

- **Output/<CaseName>/Baseline/results.csv** -- results of the last `run_case.py` run of that case (this case's own hand-set rule). `optimize_layer2.py` (Layer 2, below) adds a sibling `Output/<CaseName>/Optimised/` folder for the optimized rule, so the two are directly comparable.
- **Plot/<CaseName>/Baseline/routing_result.png** -- plot of the last `run_case.py` run of that case

## Safety thresholds: dam crest level vs. max flood/design level

Two related but different things, both optional in `scalars.csv`:

- **`dam_crest_level`** [m a.s.l.] -- the actual physical top of the
  dam. The reservoir level reaching this is true overtopping (a
  structural failure risk).
- **`max_flood_level`** [m a.s.l.] -- a stricter design/regulatory
  maximum reservoir level, kept below the crest by a freeboard margin.
  This is the level that should never be reached even under the design
  flood, precisely so real overtopping never gets close. Set this
  *below* `dam_crest_level`, not equal to it.

If either or both are set, `run_case.py` reports after every run
whether the peak level stayed below the design level, exceeded the
design level but not the crest (a design safety margin violation, not
yet true overtopping), or exceeded the crest itself (true overtopping)
-- and draws both as reference lines on the plot. Neither is required;
omit them and the case just runs without a safety check, as before.

These same two levels are exactly what a future rule-curve optimizer
would use as its no-overtopping constraint: `max_flood_level` as the
target the optimized operating rule must keep the reservoir under.

## Where outlet discharge is defined

`Data/<CaseName>/case_config.py`, inside its `build_outlets()` function,
is where every outlet for a case is defined. The actual discharge
formulas are classes in `Module/outlets.py`:

| Class | Discharge law | Typical use |
|---|---|---|
| `FreeOverflowSpillway` | `Q = C * Leff * (H - crest)^1.5`, with `Leff = L - 2*N_piers*Kp*(H-crest)` (pier contraction; set `N_piers=0` for none) | Ungated weir/spillway |
| `GatedSpillway` | `Q = Cd(H,a) * width * a * sqrt(2g*H_eff)` while gated; switches to the free-overflow weir formula once `a >= free_flow_ratio*(H-sill)` (gate lifted clear of the flow, default ratio 0.6). Supports `n_gates_total`/`n_gates_operational` for identical parallel gates and the N-1 design rule. | Underflow/orifice gate |
| `FuseGate` | Pre-trip: `Q = pre_trip_C * pre_trip_length * (H - pre_trip_crest_level)^1.5` (often labyrinth-shaped, overtopping before it tips); once `H > trigger_level`, irreversibly switches to `Q = post_trip_C * post_trip_length * (H - post_trip_sill)^1.5` | Fuse plugs / fusible levee segments |
| `BottomOutlet` | `Q = opening * Cd(H,opening) * A * sqrt(2g*(H-invert))`. Supports `n_outlets_total`/`n_outlets_operational` for identical parallel outlets and the N-1 design rule. | Submerged low-level orifice |

`Data/Template/case_config.py` has one illustrative, worked example of
each type, with placeholder numbers -- not calibrated data for any real
dam. Edit the numbers, delete the types you don't need, or duplicate a
block for a second outlet of the same type.

`GatedSpillway`'s two-regime per-gate formula is also available
standalone as `gated_spillway_unit_discharge(H, a, sill_level, width,
discharge_coefficient, free_flow_ratio, free_flow_C, ...)` -- the same
function `GatedSpillway.discharge()` itself calls internally, exposed
for cases where you need to compute a gate's discharge from an
EXPLICIT opening without an `operating_rule`/outlet object wrapping
it, most notably inverting it to solve "what opening gives this target
discharge" (see `sequential_fill_rule()` below).

### Gated structures: operating rule vs. discharge coefficient

For `GatedSpillway` and `BottomOutlet`, TWO SEPARATE THINGS decide the
discharge, and it's important to keep them separate:

1. **Operating rule** `(H, t) -> gate opening` -- decides WHAT THE
   OPENING SHOULD BE right now, given the reservoir level (and later,
   a forecast). This is control/PLC logic. `outlets.py` provides
   `level_trigger_rule(H_open, H_full, a_max, a_min=0.0)`, the standard
   level-triggered piecewise-linear rule; a future forecast-aware rule
   can be swapped in without touching any outlet class.

   `a_min` (optional, default 0.0): many real gates must never be left
   barely cracked open under head -- resting near-closed accelerates
   seal wear and can cause cavitation/vibration -- so the rule JUMPS
   from fully closed (0) straight to `a_min` the instant `H` exceeds
   `H_open`, rather than creeping through a continuous range of
   near-zero openings. Set it to your gate's minimum-safe-opening spec;
   omit it (or leave 0.0) for a gate without that constraint.

   If your real operating procedure is instead a fixed STEP TABLE
   ("at 117.5 m, open to 0.5 m; at 118.1 m, to 2 m; ...") rather than a
   continuous ramp -- common in real dam operation manuals -- use
   `staged_trigger_rule(stages)` instead, with `stages` a list of
   `(H_threshold, opening)` pairs read straight off that table.

   If the real operating procedure is REACTIVE rather than a fixed
   table at all -- e.g. "hold the reservoir at a target level by
   releasing no more than what's coming in" -- neither of the above
   fits, since they only ever see the reservoir level, not the inflow.
   Use `sequential_fill_rule(gate_index, n_gates, target_release_fn,
   sill_level, width, a_min, a_max, ...)` instead: it inverts the gate
   discharge law (given a target TOTAL release, solve backward for the
   opening) and brings `n_gates` IDENTICAL gates into service ONE AT A
   TIME, in a fixed order, each ramping `a_min` -> `a_max` before the
   next starts -- so at most one gate is ever partway open. Build one
   `operating_rule` per gate (same `gate_index`/`n_gates`/
   `target_release_fn` passed to each, differing only in `gate_index`);
   `target_release_fn(H, t)` is typically something like `lambda H, t:
   inflow.discharge(t) if H > FSL else 0.0` for a "pass through the
   flood once above a target level" rule. To read `Qin(t)` inside
   `case_config.py`, load your own `Hydrograph` from
   `inflow_hydrograph.csv` there -- `operating_rule`'s `(H, t)`
   signature has no other way to see the inflow.

   Note on the underlying discharge law: `gated_spillway_unit_discharge()`
   (the standalone per-gate formula `GatedSpillway` itself is built on)
   has a genuine STEP DISCONTINUITY at the gated/free-flow transition --
   once a gate clears the flow path, opening it further adds no more
   discharge, so achievable discharge can jump from one value straight
   to a higher one with nothing achievable in between.
   `sequential_fill_rule()` handles this conservatively (if a target
   falls in that gap, it settles on the top of the gated regime rather
   than overshooting into free flow), but it's worth knowing about if
   you're writing your own inversion logic against the same discharge
   law.

   **Multiple gates, and whether they're synchronized**: `GatedSpillway`
   has `n_gates_total`/`n_gates_operational` for several IDENTICAL gate
   bays operated THE SAME WAY AT THE SAME TIME (see its docstring for
   the N-1 design-rule use of `n_gates_operational`). But many real
   multi-gate spillways instead open gates ONE AT A TIME in a staggered
   sequence, not synchronously -- using `n_gates_total` for a spillway
   that's actually staggered can substantially OVER-predict discharge
   in the partially-open range (non-conservative for a flood study, since
   it understates the peak level a real, slower-opening spillway would
   allow). If you have a gate-by-gate operating table, model it as ONE
   `GatedSpillway` instance PER GATE (same physical/rating parameters --
   they're still identical gates, just not synchronized), each with its
   own `staged_trigger_rule()` reading that gate's own trigger levels.
   Cross-check the total discharge against the source table across the
   full head range before trusting it.

2. **Discharge coefficient** `Cd`, or `Cd(H, opening)` -- decides HOW
   MUCH WATER ACTUALLY FLOWS through that opening, given the physical
   hydraulics. A constant is fine as a first approximation; pass a
   callable `Cd(H, opening) -> float` instead if your rating data shows
   the coefficient varying with head and/or opening (e.g. distinguishing
   free vs. submerged flow regimes based on the ratio opening/head).

Keeping these separate means you can change a gate's control logic
(which is exactly what optimizing operating rules, in a later phase,
will do) without touching its physical rating curve, and vice versa.

### Multiple identical gates/outlets, and the "N-1" design rule

`GatedSpillway` and `BottomOutlet` both support `n_gates_total`/
`n_gates_operational` (or `n_outlets_total`/`n_outlets_operational`)
for the common case of several identical, interchangeably-operated
gates or outlets at one structure. Total discharge is simply the
per-gate/per-outlet discharge times however many are set operational.

This directly covers the "N-1" design/safety rule used by many
countries' dam safety guidelines: for the design discharge check,
assume the single most capacious gate is unavailable, and verify the
remaining gates still pass the design flood. Set
`n_gates_operational = n_gates_total - 1` to model that.

If your gates/outlets are NOT identical (different widths/capacities),
this shortcut doesn't apply cleanly -- instead, define one outlet
instance per physical gate, and set `on_off=False` on whichever
specific (largest) one you want to take out of service for the check.

**`on_off=False` vs. a gate that's merely stuck closed**: `on_off=False`
means ZERO discharge from that outlet, always -- right for a gate
that's been physically removed or fully isolated, but NOT right for a
gate that's stuck closed due to a mechanical/electrical fault or
maintenance -- such a gate usually can still be passively OVERTOPPED
like a fixed weir if the reservoir rises high enough over the top of
the (immobile) leaf. Model that instead as a `FreeOverflowSpillway`
with `crest_level = sill_level + <the gate's own physical leaf
height>` (the top of the closed gate leaf) and its OWN discharge
coefficient -- deliberately NOT the same value as `free_flow_C`, which
is calibrated for a gate intentionally lifted clear of the flow, a
different hydraulic geometry than water passing over a fixed closed
leaf.

Use a SEPARATE, explicitly-named constant for the gate's physical leaf
height (e.g. `GATE_1_LEAF_HEIGHT` in `Data/Template/case_config.py`'s
`gate_1_stuck_example`) -- deliberately NOT the operating rule's own
`a_max`, even though the two often coincide (a common, valid design
for a vertical-lift gate whose full travel clears it completely from
the flow path). `a_max` is how far the PLC will ever COMMAND the gate
to lift -- an operational limit; leaf height is a physical dimension of
the gate itself. Reusing one number for both is an easy way to
silently get the wrong overtop crest if your real gate's maximum
travel and its leaf height ever differ, or if one gets edited later
and the other doesn't (this happened for real on this project -- see
chat history). If you're confident they're the same for your actual
gate, set both constants to the same value explicitly; that's a fine,
common outcome, just make it a stated fact about your gate rather than
an implicit, reused number. See `Data/Template/case_config.py`'s
commented-out `gate_1_stuck_example` for the worked pattern.

## Why outlets are Python objects, not CSV rows

Stage-surface curves, hydrographs, and scalars are plain tables --
CSV is a natural fit. Discharge laws and gate operating rules are
functions with conditional logic (weir/orifice equations, level
triggers, irreversible fuse-gate state) -- trying to encode that in
CSV/YAML means inventing a formula mini-language, which is more
fragile than just using Python directly. So `case_config.py` is where
you list a dam's outlets as objects, kept alongside that case's own
data so a new scenario is a single folder duplication.

## Layer 2 -- gate operating-rule optimization

```
python Module/optimize_layer2.py <CaseName> [--pop-size N] [--n-gen N] [--seed N]
```
e.g. `python Module/optimize_layer2.py Template --pop-size 40 --n-gen 30`

Tunes the level-triggered gate operating rule -- `H_open`, `H_full`,
`a_min`, `a_max` per gate -- to delay downstream releases as long as
possible without exceeding `max_flood_level`. It never touches
physical/rating parameters (crest levels, widths, discharge
coefficients, N-1 counts); only the control logic. See "Gated
structures: operating rule vs. discharge coefficient" above for why
that split matters here, and for what `a_min` protects against.

**Problem formulation**
- Decision variables: whatever `case_config.optimizable_gates()`
  declares, per group -- NOT hard-coded to a fixed shape. The common
  case is `(H_open, H_full, a_min, a_max)` per gate (every listed gate
  needs bounds for all four, even if some are degenerate, e.g.
  `"a_max": (1.0, 1.0)` to keep a fully-open ceiling fixed) -- but a
  case can declare a completely different parameter set instead, e.g.
  a single `"shift"` that slides an entire `staged_trigger_rule()` step
  schedule earlier/later while preserving a real staggered multi-gate
  sequence (useful when the gate-to-gate stagger itself is a deliberate
  design not meant to be searched over independently). Whatever names
  a group declares are simply the keys its own
  `build_outlets(rule_overrides)` needs to know how to use.
- Objective: maximize the time before total downstream release
  (`Qout_total + Qwithdrawal` -- edit `downstream_release()` in
  `optimize_layer2.py` if withdrawal doesn't rejoin the downstream
  channel for your dam) exceeds `downstream_threshold_m3s`, a new
  scalar you set in `scalars.csv`
- Hard constraint: peak reservoir level must stay at or below
  `max_flood_level`
- Secondary constraint: for any group that declares BOTH `H_open` and
  `H_full` (or both `a_min` and `a_max`), `H_full` must exceed
  `H_open` (and `a_max` must exceed `a_min`) by at least a small margin
  (`min_gate_gap`, default 0.2) -- enforced as constraints rather than
  folded into the search bounds, so bounds can overlap freely in
  `optimizable_gates()`. Groups using a different parameterization
  (like a shared "shift") have nothing checked here -- any value
  within their own declared bounds is inherently valid.
- Optimizer: `pymoo`'s NSGA-II (`pip install pymoo`)

**Making a gate optimizable in a new case**: in that case's
`case_config.py`,
1. Build the gate's `operating_rule` from a small dict of defaults
   updated by `rule_overrides.get("<gate_name>", {})`, the same
   pattern `gate_1` and `bottom_outlet_1` use in
   `Data/Template/case_config.py` (for `level_trigger_rule`'s 4
   parameters) -- or, for a different parameterization (e.g. a shared
   `shift`), whatever pattern fits your own `build_outlets()`
2. Add `optimizable_gates()` returning `{"<group_name>": {"<param>":
   (lo,hi), ...}}` for every group you want tuned, with whatever
   parameter names that group's `build_outlets()` knows how to use
3. Add `downstream_threshold_m3s` to that case's `scalars.csv`

Only groups wired up this way are touched -- any other outlet
(ungated spillway, fuse gate, always-open bottom outlet, or a gate
you deliberately leave out) is left exactly as `case_config.py`
defines it.

Writes results for the OPTIMIZED rule only, in the exact same format
`run_case.py` uses (via the shared `Module/report.py`), to:
- `Output/<CaseName>/Optimised/results.csv`, `Plot/<CaseName>/Optimised/routing_result.png`
- `Output/<CaseName>/Optimised/best_rule.txt` -- the winning rule
  parameters, plus peak level / delay-to-threshold for both the
  optimized rule and the baseline, for a quick before/after read
  without opening either plot.

The `Baseline/` folder (this case's own hand-set rule, i.e.
`rule_overrides=None`) is owned by `run_case.py`, not this script --
run `python Module/run_case.py <CaseName>` (once, or again after
editing `case_config.py`) to (re)generate
`Output/<CaseName>/Baseline/results.csv` and
`Plot/<CaseName>/Baseline/routing_result.png`. `optimize_layer2.py`
still re-simulates the baseline internally to print/report the
before-vs-after comparison, it just doesn't write those two files
itself, so the two scripts don't fight over who owns `Baseline/`.

`argparse` is optional everywhere: run `python Module/optimize_layer2.py`
with no arguments to be prompted for the case, population size (default
40), and number of generations (default 30) interactively, the same way
`run_case.py` prompts for a case when run without one.

## Layer 3 -- Monte Carlo / uncertainty analysis

```
python Module/mc_layer3.py <CaseName> [--n-outer N] [--n-inner N] [--rule baseline|optimized] [--seed N]
```
e.g. `python Module/mc_layer3.py Template --n-outer 20 --n-inner 50`

A single `run_case.py` simulation answers "what happens for this one
specific flood, assuming these specific gate/reservoir parameters."
Layer 3 answers the more useful question: given that neither the
flood's true size nor the physical structures' exact performance is
known with certainty, what's the realistic SPREAD of possible peak
reservoir levels? It runs the same solver many times over, each time
with slightly different, still-realistic inputs, and reports the
resulting distribution rather than one number.

**Two nested loops, deliberately kept separate:**
- **Outer loop** -- "which flood magnitude are we facing": by default
  (`mc_outer_hydrograph_source="scaled"`), draws a peak scale factor
  and a volume scale factor, applied to this case's own
  `inflow_hydrograph.csv`. Independent by default (see "Peak-volume
  dependence" below) -- optionally correlated via a copula if a case
  opts in. A case can opt into a fundamentally different, EMPIRICAL
  alternative instead (`mc_outer_hydrograph_source="ensemble"`) --
  see "Empirical hydrograph ensembles" below.
- **Inner loop** -- "how do the physical structures perform, for a
  fixed flood": for each outer draw, re-runs the simulation many times,
  varying discharge/rating coefficients, the reservoir's initial level
  (H0), and the stage-storage curve, all within realistic bounds. H0 is
  drawn from a NORMAL distribution (additive, can go either direction)
  centered on this case's own `scalars.csv` H0, sigma 0.10 m by default
  (`H0_DEFAULT_SIGMA`); the reservoir curve gets a log-normal multiplier,
  median 1.0, CV 0.05 by default (`AREA_DEFAULT_CV`) -- both overridable
  per-case via `mc_h0_sigma`/`mc_area_cv` scalars, same pattern as the
  outer loop's CV overrides above. Discharge-coefficient perturbation
  only happens if the case defines `mc_uncertain_params()` (below); with
  no case-specific data at all, H0 and the reservoir curve still vary,
  just not the coefficients.

Every (outer, inner) combination is one full `solver.run_simulation()`
call. Each records: peak reservoir level (the primary metric, compared
against `max_flood_level`/`dam_crest_level`), peak downstream release,
time to exceedance of `downstream_threshold_m3s` (if set), and
exceedance flags.

**Every case must explicitly declare its sampling sources -- no silent
default:** unlike earlier versions of this tool, Layer 3 no longer
falls back to a generic, unlabeled assumed spread when a case doesn't
provide real data. Every case's `scalars.csv` MUST set
`mc_outer_peak_source` and `mc_outer_volume_source` explicitly (see
below); `mc_layer3.py`'s `validate_mc_sources()` hard-errors, before
any simulation runs, if either is missing, misspelled, or missing the
data/scalar its chosen value requires. This is deliberate -- a case
that hasn't decided yet should fail loudly at startup, not run
silently against an unlabeled generic number.

**What "median" means, and why it's (almost) always locked to 1.0:**
this case's own `inflow_hydrograph.csv` IS the design hydrograph -- its
peak and volume already ARE the target-return-period central estimate.
So `peak_scale`/`volume_scale` are drawn log-normally with median = 1.0
(this case's own hydrograph, unchanged) for every source EXCEPT
`mc_outer_peak_source="bootstrap"` -- external curve/table files are
CV-ONLY inputs; they inform the SPREAD around the hydrograph, they
never recompute or shift the central value away from it. `"bootstrap"`
is the one deliberate exception (see below) -- an FFA stress-test
case's whole point is that the FFA tool's own median is a genuinely
different, disputed central estimate, not a refinement of the
hydrograph's. Log-normal rather than normal because scale factors are
strictly positive and the distribution should be symmetric in relative
(percentage) terms rather than absolute terms -- a factor of 0.85 and a
factor of 1.15 are treated as equally likely, not a factor of -0.15 and
+0.15 on a scale that shouldn't go negative.

**`mc_outer_peak_source`** -- one of three values:
- `"curve"` -- median: this case's own `inflow_hydrograph.csv`,
  unchanged. Spread (CV): from `alt_studies_csv`, which is REQUIRED
  alongside `curve_csv` for this source (`curve_csv` alone is a single
  value per return period with no spread information in it --
  `alt_studies_csv`, 2+ independent past studies compared against the
  adopted curve, is what actually supplies the CV). Selecting `"curve"`
  without a usable `alt_studies_csv` is a hard error, not a silent
  fallback -- see `Data/Template/case_config.py`'s `mc_outer_distribution()`
  comments.
- `"bootstrap"` -- median AND spread BOTH from a real FFA (Flood
  Frequency Analysis) tool's bootstrap confidence interval (or
  model-averaged quantiles in the same shape), via
  `mc_outer_distribution()`'s `bootstrap_ci_csv` key. The median here
  DELIBERATELY differs from this case's own hydrograph peak -- that
  disagreement is the entire reason this source exists. Almost always
  used for a dedicated STRESS-TEST case duplicated from the primary
  one, not the primary case itself -- see "FFA bootstrap stress-test
  cases" below.
- `"user"` -- no file needed at all; median stays this case's own
  hydrograph peak, spread = `mc_outer_peak_cv` (required in
  `scalars.csv` when this source is selected) as a pure
  engineering-judgment number.

**`mc_outer_volume_source`** -- one of two values:
- `"table"` -- median: this case's own `inflow_hydrograph.csv`,
  unchanged. Spread (CV): from `volume_duration_csv`, which is
  SELF-SUFFICIENT for this -- its own multiple return-period rows
  already carry the spread signal (via the characteristic V/Q duration
  across rows), unlike `curve_csv`, which needs a companion
  `alt_studies_csv`.
- `"user"` -- no file needed; spread = `mc_outer_volume_cv` (required
  when this source is selected).

**Overriding a derived CV while keeping the real median:** even when
`mc_outer_peak_source="curve"`/`"bootstrap"` or
`mc_outer_volume_source="table"`, an explicit `mc_outer_peak_cv`/
`mc_outer_volume_cv` in `scalars.csv` still overrides the data-derived
CV (median stays whatever that source dictates) -- i.e. "use the real
curve/table/bootstrap data, but I want to set the SPREAD myself" is a
supported combination, not all-or-nothing.

**Peak-volume dependence:** a single uniform scale factor on the whole
hydrograph moves peak and volume together, silently assuming they're
perfectly correlated. That's rarely true -- real design-hydrograph
studies often use a shorter, more intense storm duration for the most
extreme events and a longer one for more frequent floods, which changes
the volume-to-peak ratio, not just the magnitude. For any reservoir
that provides meaningful attenuation, volume (not just peak) is often a
first-order control on the routed peak level, so `ScaledHydrograph` (in
`mc_layer3.py`) supports independent `peak_scale` and `volume_scale`,
achieved by stretching the hydrograph's time axis around its own
time-to-peak rather than just rescaling its magnitude.

Full independence is itself a real, documented limitation, though --
peak and volume come from the SAME storm event, so flood-frequency
literature treats them as physically correlated, not independent (see
Requena, Mediero & Garrote, 2013, *HESS* 17:3023-3038, the direct
precedent for this project's own peak-volume -> synthetic hydrograph ->
reservoir routing -> overtopping-risk chain). `mc_peak_volume_dependence`
in `scalars.csv` is OPTIONAL and defaults to `"independent"` (the
behavior above, unchanged, no row needed) -- but a case can opt into
copula-based joint sampling instead:

- `"independent"` -- default; no row needed at all.
- `"gaussian"` -- correlated via a Gaussian copula.
- `"gumbel"` -- correlated via a Gumbel-Hougaard copula, which has
  UPPER TAIL DEPENDENCE (extreme peak and extreme volume co-occur MORE
  than under `"gaussian"` at the same `tau`) -- more often the
  literature-preferred fit for this specific pair, and the more
  conservative choice for a dam-safety application, since it doesn't
  understate exactly the tail that drives crest-exceedance risk.

Either non-default choice requires `mc_peak_volume_tau` -- Kendall's
tau, `0 <= tau < 1` -- the quantity you'd estimate from a case's own
paired annual-maximum peak/volume record (the same record an FFA tool
would use), via a rank-correlation statistic across historical
peak/volume pairs. **This is a genuinely different piece of evidence
from either marginal's own CV** -- an FFA bootstrap fit to the peak
series ALONE, or a duration/volume table, tells you nothing about how
peak and volume co-vary; that requires the two measured together,
event by event. Literature-reported values commonly fall around
0.4-0.7 for peak-volume pairs across various basins -- a starting
REFERENCE range only, not a default to assume for any specific case
without a paired record to check it against.

Marginal medians/CVs (`mc_outer_peak_source`/`mc_outer_volume_source`,
above) are completely UNCHANGED by this setting either way -- the
copula only links two marginals that already exist; it doesn't
substitute for having a volume CV in the first place, and it doesn't
touch how the median or CV of either marginal was derived. See
`Module/mc_layer3.py`'s `_sample_gaussian_copula_normals()`/
`_sample_gumbel_copula_uniforms()` docstrings for the sampling
algorithms (Marshall-Olkin/Chambers-Mallows-Stuck for Gumbel).

**Empirical hydrograph ensembles:** everything above (`mc_outer_peak_source`,
`mc_outer_volume_source`, `mc_peak_volume_dependence`) is a PARAMETRIC
approach -- a marginal CV for each of peak and volume, and, if
correlated draws are wanted, a copula family plus a `tau` that's often
an unfitted literature-range guess absent a paired record. A case can
opt into a fundamentally different, EMPIRICAL alternative instead:

```
mc_outer_hydrograph_source,ensemble,-
```

OPTIONAL, defaults to `"scaled"` (everything above, unchanged) if the
row is omitted -- every case built before this option existed keeps
working exactly as it always has. When set to `"ensemble"`,
`mc_outer_peak_source`/`mc_outer_volume_source`/`mc_peak_volume_dependence`
are NOT read or required at all -- each outer scenario instead uses a
WHOLE, already-coherent hydrograph pulled directly from a pre-generated
ensemble file, via `mc_outer_distribution()`'s new `ensemble_csv` key
(a wide CSV: `Time_hours,rep_1,rep_2,...,rep_N`, one row per timestep,
one column per replicate).

The motivating case: the companion `Design_Flood_IIUNAM` project's
`Module/task7_hydrograph_ensemble.py`, which generates hundreds of
individually-coherent hydrographs via year-block bootstrap of an
annual-maximum record (the same resampled year driving every duration
at once, so peak, volume, AND shape stay correlated correctly and
EMPIRICALLY in every replicate -- not from an assumed copula parameter).
Its own wide-format hourly CSV output is exactly this section's
`ensemble_csv` format, ready to point at directly with no conversion.
Each outer scenario draws a DIFFERENT replicate, WITHOUT replacement,
so every scenario stays genuinely independent -- this project's own
SE/convergence calculations (`_cluster_bootstrap_stats()`, Eq.176/177)
assume `n_outer` independent scenarios, and resampling with replacement
would silently reintroduce duplicates and undermine that (the same
class of bug this project already found and fixed once, for the
nested outer/inner draw structure itself). If a run requests more
`n_outer` than the ensemble has replicates, Layer 3 CAPS `n_outer` at
the available count and prints a clear warning, rather than either
erroring outright or silently resampling with replacement -- generate
more replicates from `Design_Flood_IIUNAM`'s Task 7 if you need more
than that.

`peak_scale`/`volume_scale` still appear in `mc_results.csv`/
`mc_summary.txt` in this mode, for compatibility with existing
plotting/analysis tooling -- but they're EMPIRICAL, INFORMATIONAL
summary statistics of the selected replicates relative to this case's
own `inflow_hydrograph.csv`, not sampling inputs; the actual forcing
hydrograph for each outer scenario is the ensemble replicate itself,
used directly, never reconstructed from those ratios. See
`Data/Template/case_config.py`'s
`mc_outer_hydrograph_source="ensemble"` comments for a worked
`mc_outer_distribution()` example, and `Module/mc_layer3.py`'s
`load_hydrograph_ensemble()`/`select_ensemble_hydrographs()`/
`EnsembleHydrograph` docstrings for the implementation.

**Opting a case into real CV data, on top of the mandatory source
selection above:** `mc_uncertain_params()` and `mc_outer_distribution()`
are two optional hooks in a case's own `case_config.py`, mirroring
Layer 2's `optimizable_gates()` pattern, that supply the REAL data
behind a `"curve"`/`"bootstrap"`/`"table"` source selection:

- **`mc_uncertain_params()`** -- declares which `physical_overrides`
  keys the inner loop may sample (whatever this case's own
  `build_outlets(physical_overrides=...)` already understands -- see
  its own docstring for the exact keys), each with a distribution
  (`lognormal_cv` or `normal`). Without it, the inner loop still varies
  H0 and the reservoir curve, just not discharge coefficients.
- **`mc_outer_distribution()`** -- required whenever `mc_outer_peak_source`
is `"curve"`/`"bootstrap"`, or `mc_outer_volume_source` is `"table"`
(not needed at all if both sources are `"user"`). Points at whichever
of these CSVs the selected sources need:
  - `curve_csv` (`return_period_years,peak_Q_m3s[,source]`) -- your
    adopted design flood-frequency curve. Used only as the comparison
    reference for `alt_studies_csv` (interpolated at each alt study's
    own return period) -- it does NOT set the outer loop's median;
    that stays this case's own hydrograph (see "What 'median' means"
    above).
  - `alt_studies_csv` (`return_period_years,study,peak_Q_m3s`) --
    independent past estimates of the same flood(s). REQUIRED
    alongside `curve_csv` when `mc_outer_peak_source="curve"` -- the
    spread between them and your adopted curve is what sets the peak
    scale factor's CV. Needs 2+ rows; fewer than that is a hard error,
    not a silent fallback.
  - `bootstrap_ci_csv` (`T,lower,median,upper`) -- an FFA tool's
    bootstrap confidence interval (or model-averaged quantiles in the
    same shape). Required when `mc_outer_peak_source="bootstrap"`:
    supplies BOTH the peak scale factor's median (this file's own
    `median` column at `target_return_period`) AND its CV (from the
    `(lower, upper)` width, converted to an equivalent lognormal CV at
    an optional `bootstrap_ci_confidence`, default 0.95). See "FFA
    bootstrap stress-test cases" below.
  - `volume_duration_csv`
    (`return_period_years,storm_duration_hours,peak_Q_m3s,volume_Mm3`)
    -- your adopted design hydrographs' PAIRED peak and volume by
    return period. Required when `mc_outer_volume_source="table"`:
    self-sufficient for the volume scale factor's CV (median stays
    this case's own hydrograph volume, unchanged).
  - `target_return_period` -- required only when
    `mc_outer_peak_source="bootstrap"` (curve/table sources are CV-only
    and no longer need a specific T to anchor a median lookup at).

  See `Data/Template/case_config.py`'s comments for the exact function
  signature and file formats to copy.

**FFA bootstrap stress-test cases:** when a real Flood Frequency
Analysis tool's bootstrap fit to a case's own annual-maximum record
disagrees sharply with its adopted design curve -- typically because a
short record is being extrapolated far beyond its length, and the
FFA tool's own model-averaged-disagreement diagnostics confirm that's
sample-skew instability rather than new hazard evidence -- the
recommended approach is NOT to replace the primary case's adopted
curve. Duplicate the case folder into a second, clearly-labeled
STRESS-TEST case whose `mc_outer_distribution()` sets
`mc_outer_peak_source="bootstrap"` + `bootstrap_ci_csv`, and run both
side by side, reporting both, rather than picking a winner or blending
them into one number. The stress-test case deliberately leaves
`mc_outer_volume_source` exactly as the primary case has it (usually
`"table"`, inheriting the same `volume_duration_csv` unchanged) --
the FFA bootstrap is peak-only evidence, so volume uncertainty stays
exactly as uncertain (or un-uncertain) as it already was, rather than
being artificially widened or narrowed on the strength of a change to
a completely different quantity.

**Important caveat on `alt_studies_csv`-derived uncertainty:** the
resulting CV is a measure of how much INDEPENDENT PAST STUDIES of the
same structure have disagreed with each other over time -- a
reasonable, real-data starting point, but NOT a formal statistical
confidence interval for the flood-frequency curve itself. A rigorous
confidence interval for a return period far beyond the length of the
measured record (which is the usual situation for rare design floods)
comes from the sampling uncertainty of the fitted distribution's own
parameters -- e.g. a bootstrap: resample the annual-maximum series with
replacement, refit the distribution, and re-extrapolate, many times
over, then look at the spread of results. That is meaningfully
different from (and often noticeably WIDER than) inter-study
disagreement, particularly for return periods that require extrapolating
many multiples beyond the observed record length. If the raw
annual-maximum series behind your adopted curve is available,
implementing that bootstrap in place of (or alongside) the
`alt_studies_csv` comparison would be a real improvement -- flag this
explicitly in any downstream reporting until it's done, since the
current CV likely UNDERSTATES the true extrapolation uncertainty.

**Gate-availability (failure-to-open) risk:** a third opt-in hook, in
the same spirit as `mc_uncertain_params()`, but modeling a discrete
failure event rather than a continuous physical/rating perturbation --
whether each of a case's named gates actually opens at all this draw,
not just how well it performs given that it does.

- **`mc_gate_availability(source, p_fail=None)`** -- can be populated
  two ways, both producing the SAME dict shape; the inner loop doesn't
  know or care which one a case used. `mc_layer3.py` calls this hook
  with BOTH arguments every time, driven by `scalars.csv`'s
  `mc_gate_reliability_source` (`"ci"`\|`"flat"`, REQUIRED whenever
  this hook is defined -- `validate_mc_sources()` hard-errors
  otherwise) and `mc_gate_p_fail` (required only when
  `mc_gate_reliability_source="flat"`):

  **Simplest form** -- a single, hand-written judgment-call rate,
  applied EQUALLY to every named gate (no common-cause mechanism):

  ```python
  def mc_gate_availability(source: str, p_fail: float | None = None):
      gates = ["gate_1", "gate_2", "gate_3", "gate_4", "gate_5", "gate_6"]
      if source == "flat":
          return {"gates": gates, "p_fail": p_fail}
      ...  # source == "ci", see below
  ```

  **CI-evidence form** -- derived from real Condition Index evidence
  (an inspection, or a documented engineering judgement with an
  uncertainty band) via `Module/reliability/gate_reliability.py`'s
  `build_gate_availability()`, instead of writing the dict by hand.
  Reads per-gate, per-subsystem CI estimates from a `components.csv`,
  combines them through a series fault tree into one failure-on-demand
  probability per gate, and optionally reads a
  `common_cause_events.csv` for named, shared-vulnerability branches
  (a shared control link, a shared component model/batch, a shared
  power source) that can fail multiple gates at once:

  ```python
  from reliability.gate_reliability import build_gate_availability

  def mc_gate_availability(source: str, p_fail: float | None = None):
      gates = [...]
      if source == "flat":
          return {"gates": gates, "p_fail": p_fail}
      case_dir = os.path.dirname(os.path.abspath(__file__))
      return build_gate_availability(
          os.path.join(case_dir, "reliability", "components.csv"),
          os.path.join(case_dir, "reliability", "common_cause_events.csv"),
      )
  ```

  See `Module/reliability/README.md` for the full methodology (what CI
  means, how a `(CI_mean, CI_sd)` estimate becomes a probability, and
  what's deliberately not modeled yet -- dormancy/PSSD, deterioration,
  Bayesian updating), and `Data/Template/reliability/*.csv` for the
  illustrative input format.

  **Switching between the two forms** is now just editing
  `scalars.csv`'s `mc_gate_reliability_source` -- no environment
  variable, no code change, and unlike the environment variable, the
  tier actually used for a run is now RECORDED in that run's own
  `scalars.csv` and echoed into `mc_summary.txt`:

  ```
  mc_gate_reliability_source,ci,-      # CI-evidence form
  mc_gate_reliability_source,flat,-    # flat rate; mc_gate_p_fail sets the rate
  ```

  A case's `mc_gate_availability()` still using the old no-argument
  signature (from the previous `USE_CI_RELIABILITY`/`GATE_P_FAIL`
  environment-variable pattern) will raise a clear migration error at
  startup -- update it to accept `(source, p_fail=None)` as shown
  above.

  **Sampling mechanism**, common to both forms -- each inner draw runs
  a two-stage process:
  1. Each declared common-cause group (CI-evidence form only) draws
     its own independent Bernoulli(`p_ccf`) trial; if it fires, every
     gate named in that group is forced to the failed state for this
     draw, regardless of that gate's own individual outcome below.
  2. Every gate is then ALSO drawn independently at its own rate
     (`p_fail_by_gate[gate]`, which may differ per gate under the
     CI-evidence form, or be one flat rate under the simplest form).

  A gate's final state for the draw is failed if EITHER stage says so.
  When there are no common-cause groups, this reduces exactly to the
  original independent-per-gate binomial mechanism
  (`P(K=k) = C(n,k)*p^k*(1-p)^(n-k)` for the number of gates `K` that
  fail to open, `n = len(gates)`).

  Only makes sense for a case whose gates are modeled as
  INDIVIDUALLY NAMED outlet instances (see "Multiple identical
  gates/outlets" above) -- a case using the bundled
  `n_gates_total`/`n_gates_operational` form has no individual gate
  identities for this hook to refer to. The case's `build_outlets()`
  also needs a `gates_out_of_service` parameter (a list of gate names
  to force out of service for one call, overriding whatever a case's
  own module-level "out of service" list defaults to for a normal
  run) -- see `Data/Template/case_config.py`'s commented worked
  `gates_out_of_service` example for the pattern, including how a
  gate that's merely stuck closed (rather than fully removed) can
  still be modeled as passively overtoppable via `on_off`/a substitute
  `FreeOverflowSpillway`, not simply zero discharge.

  Without this hook, Layer 3 assumes every gate stays operational for
  every draw, and prints a note saying so -- same "not every case
  needs every layer" fallback as the other two hooks.

**Output**, mirroring Layer 1/2's own conventions:
- `Output/<CaseName>/MonteCarlo/mc_results.csv` -- one row per draw,
  every sampled parameter plus its outcome metrics (including
  `gates_failed`, `n_gates_failed`, and `ccf_fired` -- present,
  possibly always empty/zero, regardless of whether the case declares
  `mc_gate_availability()`, so the column layout stays consistent
  across cases)
- `Output/<CaseName>/MonteCarlo/mc_summary.txt` -- headline percentiles
  (P5/P50/P95/mean, each with a standard error), exceedance
  probabilities (with their own standard error, or, if zero events
  were observed in the run, an explicit "bounded above by roughly
  3/n_outer, NOT zero" caveat instead of a bare `0.0000` -- a rare
  event producing zero observed occurrences in a finite run doesn't
  mean the true probability is zero). ALL of these standard errors --
  percentiles, mean, AND exceedance probabilities -- come from a
  CLUSTER (block) bootstrap that resamples whole outer scenarios, not
  individual draws (`_cluster_bootstrap_stats()` in `mc_layer3.py`),
  because the 15,000-ish draws in a typical run are NOT that many
  independent observations -- they're `n_outer` independent flood
  scenarios, each contributing `n_inner` CORRELATED draws (correlated
  because they share the same flood magnitude/volume). A flat, per-
  draw bootstrap or closed-form formula that ignores this structure
  understates the true standard error, sometimes severely -- confirmed
  empirically via a cross-seed comparison (two otherwise-identical
  runs differing only in `--seed` showed a ~11.5x discrepancy between
  a flat bootstrap's stated P95 SE and the actual cross-seed spread,
  resolved to ~1.3x once corrected). Also reported: exactly which
  sampling source was declared and used
  (`mc_outer_peak_source`/`mc_outer_volume_source`/
  `mc_gate_reliability_source`/`mc_peak_volume_dependence`, each an
  explicit per-case choice -- see "Every case must explicitly declare
  its sampling sources" and "Peak-volume dependence" above), the run's
  `--seed`, and -- if `mc_gate_availability()` is declared -- each
  named common-cause branch's declared probability and empirical fired
  fraction over the run, the empirical-vs-closed-form
  `P(>=1 gate failed)` sanity check, and, for the CI-evidence form
  specifically, an importance ranking (which subsystem type or
  common-cause branch is carrying the most probability weight -- see
  `Module/reliability/importance.py`'s docstring for what this
  ranking does and doesn't mean).
- Convergence check, alongside the above: minimum-realizations
  thresholds for the reported mean and P95, per the USACE RMC-
  TotalRisk Technical Reference Manual's Equations 176/177, compared
  against `n_outer` -- the TRUE EFFECTIVE INDEPENDENT sample size for
  this nested design, NOT the raw total draw count -- as a margin
  (e.g. "P95 needs >=1825 realizations -- this run: n_outer=100,
  0.055x margin"). Comparing against the raw draw count instead (as an
  earlier version of this check did) dramatically overstates how
  converged a run actually is, for exactly the same reason the SE
  calculations above needed the same nested-design correction. A
  margin below 1x is a real, actionable flag that random sampling noise
  alone could shift the reported value beyond the stated tolerance,
  not just a theoretical caveat -- and, for P95 specifically, is common
  even at draw counts that look large in raw terms; more `n_outer`
  (not more `n_inner`) is what closes this margin.
- `Plot/<CaseName>/MonteCarlo/mc_distribution.png` -- a histogram of
  peak levels pooled across all draws, alongside an exceedance curve
  (worst-to-best ranking) with each outer draw's own curve shown as a
  thin line underneath the pooled curve, so climate/frequency-curve
  uncertainty is visually distinct from physical/rating uncertainty

`--rule optimized` stress-tests an already-optimized Layer 2 rule
instead of this case's baseline hand-set one (parsed from
`Output/<CaseName>/Optimised/best_rule.txt`, if it exists); default is
`baseline`. Same interactive-prompt fallback as `run_case.py`/
`optimize_layer2.py` if run with no arguments.

### Layer 3 outer-loop-only -- isolating climate/flood-frequency spread

```bash
uv run python Module/mc_layer3_outer_only.py <CaseName> --n-outer 300
```

A separate script, `Module/mc_layer3_outer_only.py`, runs each outer
scenario EXACTLY ONCE, with every INNER-loop quantity held at its
deterministic, unperturbed BASELINE value instead of sampled:
`H0` exactly this case's own `scalars.csv` value (no `h0_sigma` shift),
`area_scale=1.0` (no stage-storage curve perturbation), every physical/
rating multiplier declared in `mc_uncertain_params()` at its nominal
value (`1.0` for `lognormal_cv` params -- the median of a CV spread
around nominal IS 1.0; the declared `mean` for `normal` params, since
those are typically SHIFTS rather than multipliers), and ALL gates
operational -- no failures modeled at all. That last point is a
deliberate simplification of this script specifically (not a claim
that gate failure is impossible), chosen so the isolated comparison
below is clean: gate-reliability uncertainty is unambiguously an INNER-
loop quantity, so it's held at its best case here, same as everything
else in that loop.

**What this answers, that a full Layer 3 run can't in isolation**: "how
much does peak level vary due to WHICH FLOOD we get, holding every
physical/reliability uncertainty at its nominal value?" Running this
script and a full Layer 3 run with the SAME `--n-outer`/`--seed`/case,
then comparing the two runs' spreads (e.g. their P95-P5 range, or their
`mc_summary.txt`/`outer_only_summary.txt` standard errors side by
side), gives a direct, empirical answer to which loop dominates a given
case's reported uncertainty -- rather than inferring it indirectly from
a single combined run, the way earlier analysis in this project's
history had to.

Reuses `Module/mc_layer3.py`'s own `build_case()`/`build_outer_draws()`
-- the SAME outer-loop derivation logic (all of `mc_outer_peak_source`/
`mc_outer_volume_source`/`mc_peak_volume_dependence`/
`mc_outer_hydrograph_source`, including ensemble-mode capping) as the
full Layer 3 run, so this script's outer scenarios are identical to
what a full run with the same `scalars.csv` would draw; this script
deliberately does not duplicate that logic a second time. Its
`_cluster_bootstrap_stats()`/`_min_realizations_*()`/`_format_margin()`
reporting machinery is reused too -- with exactly one deterministic
value per outer scenario, that machinery naturally degenerates into a
standard, textbook i.i.d. bootstrap/Eq.176-177 check over `n_outer`
independent scenarios (there's no inner-loop nesting here to correct
for, unlike full Layer 3's `n_outer`-vs-`n_total` distinction).

Writes: `Output/<CaseName>/OuterOnly/outer_only_results.csv`,
`Output/<CaseName>/OuterOnly/outer_only_summary.txt`,
`Plot/<CaseName>/OuterOnly/outer_only_distribution.png` -- a separate
`OuterOnly/` subfolder, so running this never overwrites or gets
confused with a full Layer 3 run's own `MonteCarlo/` output.

## Extensibility built in for later phases

- `outlets.py` also includes `SwitchedOutlet` and `ConstantOutlet` --
  composition helpers for building up more complex real-world behavior
  from simple pieces (e.g. a gate whose discharge law itself changes
  character once it opens, rather than just scaling by an opening
  fraction).
- Operating rules take `(H, t)`, not just `H`, so a future forecast-aware
  (FIRO-style) rule can be swapped in without touching outlet classes
  or the solver.
- Fuse gates only update their tripped state via `commit()`, called
  once per accepted RK4 step -- not from intermediate RK4 stage
  evaluations -- so the trigger can't fire spuriously from a trial
  sub-step overshoot.
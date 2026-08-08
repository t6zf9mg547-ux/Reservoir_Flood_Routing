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
    dam_crest_level, max_flood_level (see "Safety thresholds" below)
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
with `crest_level = sill_level + a_max` (the top of the closed gate
leaf) and its OWN discharge coefficient -- deliberately NOT the same
value as `free_flow_C`, which is calibrated for a gate intentionally
lifted clear of the flow, a different hydraulic geometry than water
passing over a fixed closed leaf. See `Data/Template/case_config.py`'s
commented-out `gate_1_stuck_example` for the pattern.

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
- **Outer loop** -- "which flood magnitude are we facing": draws a
  peak scale factor and, independently, a volume scale factor (see
  "Why peak and volume are sampled independently" below), applied to
  this case's own `inflow_hydrograph.csv`.
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

**Default sampling, before any case-specific data is provided:** both
`peak_scale` and `volume_scale` are drawn independently from a
LOG-NORMAL distribution with median 1.0 (i.e. "this case's own
`inflow_hydrograph.csv` is assumed correct as-is") and a coefficient of
variation (CV) of 0.15 (`OUTER_DEFAULT_CV` and `VOLUME_DEFAULT_CV` in
`mc_layer3.py`) -- a flat, generic, ASSUMED spread, not derived from
any data. Log-normal rather than normal because scale factors are
strictly positive and the distribution should be symmetric in relative
(percentage) terms rather than absolute terms -- a factor of 0.85 and a
factor of 1.15 are treated as equally likely, not a factor of -0.15 and
+0.15 on a scale that shouldn't go negative. A case can override just
the CV, without providing full curve data, via `mc_outer_peak_cv` and
`mc_outer_volume_cv` scalars in its own `scalars.csv`. Providing a real
`mc_outer_distribution()` (see below) replaces both the median AND the
CV with values grounded in that case's own data instead.

**Why peak and volume are sampled independently:** a single uniform
scale factor on the whole hydrograph moves peak and volume together,
silently assuming they're perfectly correlated. That's rarely true --
real design-hydrograph studies often use a shorter, more intense storm
duration for the most extreme events and a longer one for more
frequent floods, which changes the volume-to-peak ratio, not just the
magnitude. For any reservoir that provides meaningful attenuation,
volume (not just peak) is often a first-order control on the routed
peak level, so `ScaledHydrograph` (in `mc_layer3.py`) supports
independent `peak_scale` and `volume_scale`, achieved by stretching
the hydrograph's time axis around its own time-to-peak rather than
just rescaling its magnitude.

**Opting a case into real data, instead of the generic placeholder
spread:** by default, Layer 3 samples both loops from a generic,
flat, assumed spread around this case's own data, and prints a note
saying so -- it still runs correctly with zero extra setup. A case can
opt into something better-grounded via two optional hooks in its own
`case_config.py`, mirroring Layer 2's `optimizable_gates()` pattern:

- **`mc_uncertain_params()`** -- declares which `physical_overrides`
  keys the inner loop may sample (whatever this case's own
  `build_outlets(physical_overrides=...)` already understands -- see
  its own docstring for the exact keys), each with a distribution
  (`lognormal_cv` or `normal`). Without it, the inner loop still varies
  H0 and the reservoir curve, just not discharge coefficients.
- **`mc_outer_distribution()`** -- points at up to three CSVs so the
  outer loop can be grounded in real flood-frequency data instead of a
  flat assumed spread:
  - `curve_csv` (`return_period_years,peak_Q_m3s[,source]`) -- your
    adopted design flood-frequency curve. Sets the peak scale factor's
    median.
  - `alt_studies_csv` (`return_period_years,study,peak_Q_m3s`) --
    independent past estimates of the same flood(s), if more than one
    hydrology study exists for this structure. The spread between them
    and your adopted curve sets the peak scale factor's uncertainty
    (CV) from real data instead of an assumed number.
  - `volume_duration_csv`
    (`return_period_years,storm_duration_hours,peak_Q_m3s,volume_Mm3`)
    -- your adopted design hydrographs' PAIRED peak and volume by
    return period. Sets the volume scale factor (median AND
    uncertainty), independently of the peak scale factor.

  See `Data/Template/case_config.py`'s comments for the exact function
  signature and file formats to copy.

  **Overriding the derived spread while keeping the real median:** a
  case's `scalars.csv` can set `mc_outer_peak_cv`/`mc_outer_volume_cv`
  explicitly. If `mc_outer_distribution()` is ALSO set and derivation
  from `alt_studies_csv`/`volume_duration_csv` succeeds, the explicit
  scalar wins over the derived CV -- i.e. you can use the real
  curve/table data for the MEDIAN (the actual design-flood anchor)
  while setting the SPREAD yourself as an engineering-judgment number,
  rather than it being all-or-nothing. Leave the scalar unset to keep
  the old behavior (derived CV wins when available, generic default
  otherwise).

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

- **`mc_gate_availability()`** -- can be populated two ways, both
  producing the SAME dict shape; the inner loop doesn't know or care
  which one a case used:

  **Simplest form** -- a single, hand-written judgment-call rate,
  applied EQUALLY to every named gate (no common-cause mechanism):

  ```python
  def mc_gate_availability():
      return {
          "gates": ["gate_1", "gate_2", "gate_3", "gate_4", "gate_5", "gate_6"],
          "p_fail": 0.05,   # SAME probability applied to EVERY gate
      }
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

  def mc_gate_availability():
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

  **Switching between the two forms at run time**, without editing
  `case_config.py` each time -- useful for a quick side-by-side
  comparison. Combine both forms into one function, gated by a
  module-level flag that reads an environment variable (defaulting to
  the CI-evidence form if unset):

  ```python
  USE_CI_RELIABILITY = os.environ.get("USE_CI_RELIABILITY", "true") \
      .strip().lower() not in ("0", "false", "no")
  TIER0_P_FAIL = float(os.environ.get("GATE_P_FAIL", 0.05))

  def mc_gate_availability():
      if not USE_CI_RELIABILITY:
          return {"gates": [...], "p_fail": TIER0_P_FAIL}
      case_dir = os.path.dirname(os.path.abspath(__file__))
      return build_gate_availability(
          os.path.join(case_dir, "reliability", "components.csv"),
          os.path.join(case_dir, "reliability", "common_cause_events.csv"),
      )
  ```

  Then, from the command line, no code edits needed to switch:
  ```
  python Module/mc_layer3.py <CaseName> ...                                    # CI-evidence form (default)
  USE_CI_RELIABILITY=false python Module/mc_layer3.py <CaseName> ...           # flat rate, default p_fail
  USE_CI_RELIABILITY=false GATE_P_FAIL=0.10 python Module/mc_layer3.py <CaseName> ...  # flat rate, custom p_fail
  ```

  `Module/mc_layer3.py` itself never reads either environment variable
  and has no branching logic for this at all -- the switch lives
  entirely inside `case_config.py`'s own `mc_gate_availability()`,
  which `mc_layer3.py` just calls and consumes the result of, same as
  always. This keeps the same "case-specific choices live in
  `case_config.py`, the engine stays generic" split used everywhere
  else in this project -- see `Data/Template/case_config.py`'s
  commented worked version of this pattern.

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
  (P5/P50/P95/mean, each with a standard error -- bootstrap for the
  percentiles, closed-form SEM for the mean), exceedance probabilities
  (with their closed-form binomial standard error, or, if zero events
  were observed in the run, an explicit "bounded above by roughly
  3/N, NOT zero" caveat instead of a bare `0.0000` -- a rare event
  producing zero observed occurrences in a finite run doesn't mean the
  true probability is zero), exactly which sampling settings were used
  (generic placeholder vs. case-provided real data, for all three
  hooks), the run's `--seed`, and -- if `mc_gate_availability()` is
  declared -- each named common-cause branch's declared probability
  and empirical fired fraction over the run, the empirical-vs-closed-
  form `P(>=1 gate failed)` sanity check, and, for the CI-evidence
  form specifically, an importance ranking (which subsystem type or
  common-cause branch is carrying the most probability weight -- see
  `Module/reliability/importance.py`'s docstring for what this
  ranking does and doesn't mean).
- Convergence check, alongside the above: minimum-realizations
  thresholds for the reported mean and P95, per the USACE RMC-
  TotalRisk Technical Reference Manual's Equations 176/177, compared
  against this run's actual draw count as a margin (e.g. "P95 needs
  >=1825 realizations -- this run: 15000, 8.2x margin"). A margin
  below 1x is a real, actionable flag that random sampling noise
  alone could shift the reported value beyond the stated tolerance,
  not just a theoretical caveat.
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
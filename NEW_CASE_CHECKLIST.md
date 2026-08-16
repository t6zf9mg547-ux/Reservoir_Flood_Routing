# New Case Checklist

A step-by-step guide for onboarding a new dam/reservoir case into this
project. Work through it in order — Layer 2 and Layer 3 both depend on
Layer 1 being right first, and Layer 1 depends on the physical features
below being gathered accurately before any code gets written.

Opening prompt for a new Claude Code session, once the reference
material below is in the repo: *"New case for `<DamName>`. Read
`README.md`, `Data/Template/case_config.py`, and the reference material
in `Data/<DamName>/`. Propose a `case_config.py` design before writing
anything"* — per `CLAUDE.md`'s `Working style` rule, this should produce
a design proposal to review, not a finished file on the first pass.

## 0. Gather the project's main features and constraints FIRST

Do this before creating any file. Every item below should be a real,
sourced number or a stated placeholder — not a guess standing in for
one. If a number isn't known yet, write "TBD — see [source]" rather
than a plausible-looking placeholder; a placeholder that looks real is
easy to forget to revisit later.

**Reservoir**
- [ ] Stage-storage/surface curve source (survey? design report table?)
      and its format — becomes `reservoir_curve.csv`
- [ ] Normal operating level (`H0`)
- [ ] Dam crest level (`dam_crest_level`) — true overtopping threshold
- [ ] Max flood/design level (`max_flood_level`) — the case's own
      no-overtopping design target, distinct from the crest itself

**Spillway / gates**
- [ ] Number of gates, and whether they're identical or vary
      (width/capacity) per gate
- [ ] Sill level(s)
- [ ] Gate width per gate
- [ ] **Gate physical leaf height, confirmed separately from the
      gate's maximum OPERATIONAL opening** — these are two different
      things that are easy to conflate into one reused number (real
      bug found and fixed across every case in this project — see
      `README.md`'s `on_off=False` vs. stuck-gate section). Ask: what
      is the gate's actual manufactured/as-built leaf height? Is the
      maximum commanded opening confirmed equal to it, or could a
      mechanical restriction make them differ?
- [ ] Operating philosophy — which of these (or something else) does
      the real dam actually do:
      - a level-triggered ramp (`level_trigger_rule`) — smooth open
        between two levels
      - a fixed step table from an operations manual
        (`staged_trigger_rule`) — discrete opening steps at named
        trigger levels
      - a reactive/pass-through rule (`sequential_fill_rule`) —
        release tracks inflow, gates come into service one at a time
      Get this from the real operations manual or SCADA logic, not an
      assumption — the discharge behavior differs meaningfully between
      these (see `README.md`'s outlet/operating-rule section)
- [ ] Discharge coefficient(s) — gated/orifice regime and free-flow
      regime, and their SOURCE (a rating table fit? manufacturer data?
      a generic default?) — state explicitly if a value is a
      rule-of-thumb rather than case-calibrated (e.g. `free_flow_ratio`
      is very often not independently calibrated — say so if true here
      too)
- [ ] `a_min` — smallest allowed opening once triggered, and whether
      it's a real spec or a generic default

**Fuse gates / auxiliary spillway (if any)**
- [ ] Number of modules, module width, trigger level(s) per group
- [ ] Pre-trip crest level and coefficient (before tipping)
- [ ] Post-trip sill level and coefficient (after tipping) — confirm
      these use a genuinely different sill/coefficient than pre-trip,
      not the same numbers reused

**Bottom outlet(s) (if any)**
- [ ] Invert level, area, operating rule

**N-1 / reliability considerations**
- [ ] Is a stuck-gate/out-of-service scenario relevant for this dam?
      If so, the gate's physical leaf height (above) is what its
      overtop crest will be built from
- [ ] Is Condition Index data available for gate reliability (Layer 3),
      or will this case use a flat judgment-call rate instead?

**Downstream**
- [ ] Downstream threshold(s) for safety/exceedance reporting
      (`downstream_threshold_m3s`)

**Design flood**
- [ ] Source of the design inflow hydrograph (`inflow_hydrograph.csv`)
      — which return period, which method produced it
- [ ] Is there a real flood-frequency curve, alternate studies, or a
      duration/volume table available (for Layer 3's `curve`/`table`
      sources), or will this case start on `user`-source judgment CVs?
- [ ] Is a real FFA bootstrap fit or hydrograph ensemble available for
      this dam (Layer 3's `bootstrap`/`ensemble` sources), or is that
      a later addition?

**Put all of this — reports, rating tables, survey data — directly in
`Data/<DamName>/` in the real repo, not in a separate knowledge
space.** Claude Code reads the actual repo; keeping reference material
anywhere else risks it silently drifting from what's actually there.

## 1. Layer 1 — routing kernel

- [ ] `cp -r Data/Template Data/<DamName>` (or build fresh from the
      gathered features, but Template's structure is the fastest start)
- [ ] Populate `reservoir_curve.csv`, `inflow_hydrograph.csv`,
      `withdrawal_schedule.csv` (or leave absent — falls back to a flat
      0 automatically), `scalars.csv` (`t_max`, `dt`, `H0`, and
      optionally `dam_crest_level`/`max_flood_level`/
      `downstream_threshold_m3s`)
- [ ] Write `build_outlets()` from the gathered features above — propose
      the design first (per the opening prompt), particularly the
      operating-rule choice and the leaf-height/`a_max` split
- [ ] Run it for real: `uv run python Module/run_case.py <DamName>` —
      don't consider this step done on a syntax check alone; look at
      the actual `Output/<DamName>/Baseline/results.csv` and
      `Plot/<DamName>/Baseline/routing_result.png`
- [ ] Sanity-check: does peak level respond to peak inflow the way you
      expect? Does the routed release make physical sense against the
      gates' known capacity at that level?
- [ ] If a reference result exists to check against, update
      `ORIGINAL_CSV` in `validate_against_original.py` and run it

## 2. Layer 2 — gate operating-rule optimization

- [ ] Decide what's actually tunable for this dam's real operating
      rule — a single shared parameter (e.g. a whole-schedule shift, if
      the gate sequence itself is a fixed, deliberate design not to be
      searched independently) vs. per-gate bounds
      (`H_open`/`H_full`/`a_min`/`a_max`) if there's a real reason to
      tune them independently
- [ ] Write `optimizable_gates()` with bounds justified by where the
      case's own calibration data actually covers — don't let the
      optimizer extrapolate discharge coefficients into head ranges
      they were never checked against
- [ ] Run it for real: `uv run python Module/optimize_layer2.py <DamName>`
      and check `Output/<DamName>/Optimised/best_rule.txt` against
      physical plausibility, not just that it ran

## 3. Layer 3 — Monte Carlo / uncertainty analysis

`scalars.csv` MUST declare these explicitly or the run hard-errors —
see `README.md`'s "Layer 3" section for the full contract:

- [ ] `mc_outer_peak_source` (`curve`\|`bootstrap`\|`user`) — pick based
      on what real data was found in step 0; `user` with an honest
      judgment-call CV is legitimate if nothing better exists yet
- [ ] `mc_outer_volume_source` (`table`\|`user`)
- [ ] `mc_gate_reliability_source` (`ci`\|`flat`) — only needed if
      `mc_gate_availability()` is defined
- [ ] If `mc_outer_peak_source=curve`: also needs `alt_studies_csv`
      (curve alone has no spread information)
- [ ] If `mc_outer_peak_source=bootstrap`: needs `bootstrap_ci_csv` +
      `target_return_period` in `mc_outer_distribution()`
- [ ] Optional: `mc_peak_volume_dependence`/`mc_peak_volume_tau` — only
      set a real (non-`independent`) value if there's a paired
      peak/volume record to estimate `tau` from, or label it clearly as
      an unfitted literature-range judgment call if not
- [ ] Optional: `mc_outer_hydrograph_source=ensemble` instead of the
      above entirely, if a real hydrograph ensemble exists for this dam
      (e.g. from `Design_Flood_IIUNAM`'s Task 7)
- [ ] Write `mc_uncertain_params()` (inner-loop physical/rating
      uncertainty) and `mc_gate_availability()` (gate reliability) if
      this case needs them — both are opt-in; Layer 3 runs fine without
      either, just with less uncertainty modeled
- [ ] Smoke test first: `uv run python Module/mc_layer3.py <DamName>
      --n-outer 5 --n-inner 20` — confirm it runs and
      `mc_summary.txt` reports the sources you actually intended, not
      just that it didn't crash
- [ ] Full run once the smoke test looks right, checking the Eq.176/177
      convergence margins before trusting P95/tail numbers

## Before calling the case done

- [ ] Every physical parameter above traces to a real source, not an
      unlabeled guess
- [ ] `build_outlets()` has actually been run, not just read
- [ ] Each layer's smoke test/small run has been executed and its
      output inspected, not just its exit code checked
- [ ] `git add`/commit proposed with a concise message — wait for
      confirmation before committing or pushing, per `CLAUDE.md`

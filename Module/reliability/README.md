# Module/reliability/

Computes gate failure-on-demand (P_FOD) probabilities for Layer 3's
`mc_gate_availability()` hook, from Condition Index (CI) evidence
rather than a single guessed rate. Implements a simplified version of
the CI-based / fault-tree / common-cause methodology described in the
project's spillway-gate-reliability technical note (synthesizing
Kalantarnia et al. 2013, Kalantarnia thesis 2013, Lewin et al. 2003,
and Patra et al. 2019) -- CI-only for this first version: no failure-
rate/Weibull deterioration model, no dormancy/PSSD, no Bayesian
updating. Those are documented as deferred, not silently dropped; see
"Not yet implemented" below.

## Two tiers, one contract

`mc_gate_availability()` always returns the same shape, regardless of
which tier produced it:

```python
{"gates": [...],
 "p_fail_by_gate": {"gate_1": 0.03, "gate_2": 0.09, ...},
 "ccf_groups": [{"label": "...", "gates": [...], "p_ccf": 0.02}, ...]}
```

- **Tier 0 -- poorly documented case, no CI evidence.** Write the dict
  by hand: one flat probability copied across every gate,
  `ccf_groups` empty. No different, in spirit, from the very first
  version of this hook.
- **Tier 1 -- CI evidence available.** Call
  `gate_reliability.build_gate_availability()` (below) instead of
  writing the dict by hand; it derives the same shape from
  `components.csv`/`common_cause_events.csv`.

`Module/mc_layer3.py`'s inner loop consumes this dict exactly the same
way either way -- it has no idea, and doesn't need to, which tier
produced it.

## Files

- `condition_index.py` -- `pf_from_ci(ci_mean, ci_sd, ci_failure)`:
  the one probability primitive everything else is built on. CI is
  treated as normally distributed (mean, sd), not a fixed point value
  -- see its docstring for why.
- `fault_tree.py` -- `series_pf()`/`parallel_pf()`, and
  `gate_pfod()`: combines a gate's named subsystem failure
  probabilities into one P_FOD, in series by default (any essential
  subsystem failing prevents the gate opening).
- `common_cause.py` -- `load_ccf_groups()`: reads
  `common_cause_events.csv` into named, shared-vulnerability groups,
  each with its own CI-derived probability.
- `gate_reliability.py` -- `build_gate_availability()`: the main
  entry point a case's `mc_gate_availability()` should call.
- `importance.py` -- `rank_importance()`: which subsystem type or
  common-cause branch is carrying the most probability weight, for
  prioritizing where a limited inspection/repair budget goes next. A
  simple first-pass ranking (see its docstring), not a rigorous
  sensitivity measure.

## Input CSV schemas

**`components.csv`** -- one row per (gate, subsystem):

| Column | Meaning |
|---|---|
| `component_id` | Free-text identifier for this row |
| `gate_id` | Must match the gate names the case's own `build_outlets()` uses |
| `subsystem` | Free-text category (e.g. `hydraulic`, `structural`, `control_electrical`) -- every row sharing a `(gate_id, subsystem)` combination is one entry; different subsystems for the same gate are combined in series |
| `CI_mean`, `CI_sd` | The CI estimate (0-100 scale), supplied by whoever wrote the row -- inspection evidence or engineering judgement. **This module never invents these values.** |
| `CI_failure` | The CI value at or below which this row is considered failed. Typically one constant across a case's whole file (a single documented assumption), but each row may set its own. |
| `notes` | Free text -- what the estimate is based on |

**`common_cause_events.csv`** -- one row per named shared vulnerability:

| Column | Meaning |
|---|---|
| `event_id`, `label` | Identifier and display name |
| `affected_gates` | `ALL`, or a semicolon-separated list of gate names |
| `CI_mean`, `CI_sd`, `CI_failure` | Same meaning as above, describing the condition of the SHARED item (e.g. a comms link, a shared component batch/model) |
| `notes` | Free text |

A case with no `common_cause_events.csv` (or an empty one) simply has
no common-cause mechanism modeled -- not an error.

## Not yet implemented (deferred, not dropped)

- Failure-rate / exponential K-factor model (Section 5.1 of the
  technical note) -- needs failure-rate data most cases won't have.
- Weibull deterioration projection (Section 5.2) -- same reason.
- Dormancy / PSSD (Section 7) -- needs test-history data (successes
  and failures of actual gate-lift tests), not condition observation.
- Bayesian updating from new inspections/tests (Section 15) -- natural
  next step once a case has more than one dated inspection on file.
- Dynamically-required gate count / static insufficient-capacity check
  (Section 11) -- NOT needed here: `mc_layer3.py`'s full Monte Carlo
  flood routing already determines the actual consequence of whatever
  gates fail on a given draw (peak level, exceedance), which is more
  accurate than a separate static binomial check would be.

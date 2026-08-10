# Run Commands

Quick reference for every standard command in this project. All commands
are run through `uv run` so they use the project's managed virtual
environment (`uv sync` once, beforehand, to set that up -- see README.md).

Replace `<CaseName>` with a real case folder name under `Data/`
(e.g. a case you've duplicated from `Data/Template/`).

## Layer 1 -- single simulation

Prompts interactively for the case:
```bash
uv run python Module/run_case.py
```

Skip the prompt (useful for scripted/batch runs):
```bash
uv run python Module/run_case.py <CaseName>
```

Compare a run against a reference results file (edit `ORIGINAL_CSV` at
the top of that script to point at your reference first):
```bash
uv run python Module/validate_against_original.py <CaseName>
```

Writes: `Output/<CaseName>/Baseline/results.csv`,
`Plot/<CaseName>/Baseline/routing_result.png`.

## Layer 2 -- gate operating-rule optimization

```bash
uv run python Module/optimize_layer2.py <CaseName> [--pop-size N] [--n-gen N] [--seed N]
```
e.g.
```bash
uv run python Module/optimize_layer2.py <CaseName> --pop-size 40 --n-gen 30
```

Run with no arguments to be prompted for the case, population size
(default 40), and number of generations (default 30) interactively.

Writes: `Output/<CaseName>/Optimised/results.csv`,
`Plot/<CaseName>/Optimised/routing_result.png`,
`Output/<CaseName>/Optimised/best_rule.txt`.

## Layer 3 -- Monte Carlo / uncertainty analysis

```bash
uv run python Module/mc_layer3.py <CaseName> [--n-outer N] [--n-inner N] [--rule baseline|optimized] [--seed N]
```

Standard run (project defaults: 100 outer x 150 inner = 15,000 draws):
```bash
uv run python Module/mc_layer3.py <CaseName>
```

Custom draw counts:
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

Quick smoke test before committing to a full run:
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 5 --n-inner 20
```

Against the Layer 2-optimized rule instead of the case's hand-set baseline:
```bash
uv run python Module/mc_layer3.py <CaseName> --rule optimized
```

Reproducible run (fixed seed, e.g. for comparing before/after a change):
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150 --seed 1
```

Run with no arguments to be prompted interactively, same as `run_case.py`.

Writes: `Output/<CaseName>/MonteCarlo/mc_results.csv`,
`Output/<CaseName>/MonteCarlo/mc_summary.txt`,
`Plot/<CaseName>/MonteCarlo/mc_distribution.png`.

### Before running Layer 3 for the first time on a case

`scalars.csv` MUST set `mc_outer_peak_source`, `mc_outer_volume_source`,
and (if `mc_gate_availability()` is defined) `mc_gate_reliability_source`
explicitly -- the run refuses to start otherwise, with a message naming
exactly what's missing. `mc_peak_volume_dependence`/`mc_peak_volume_tau`
are OPTIONAL (default to fully independent sampling if omitted -- see
below). See the README's "Layer 3" section and `Data/Template/scalars.csv`'s
`options` column for the valid values and what each one needs.

### Peak-volume dependence switch

Optional for every case -- omit both rows entirely for the default,
fully independent `peak_scale`/`volume_scale` sampling (unchanged from
every prior version of this tool). To opt into copula-based correlated
sampling instead, both rows go in that case's own `scalars.csv`, same
"no command-line flag or environment variable, edit the CSV" pattern
as the gate-reliability switch above:

Gumbel-Hougaard (has upper tail dependence -- usually the better fit
for peak-volume specifically, and the more conservative choice for a
dam-safety application):
```
mc_peak_volume_dependence,gumbel,-
mc_peak_volume_tau,0.5,-
```
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

Gaussian (symmetric, no tail dependence):
```
mc_peak_volume_dependence,gaussian,-
mc_peak_volume_tau,0.5,-
```
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

`mc_peak_volume_tau` is Kendall's tau (`0 <= tau < 1`), ideally
estimated from this case's own paired annual-maximum peak/volume record
if one exists -- see the README's "Peak-volume dependence" section for
why this is a genuinely different piece of evidence from either
marginal's own CV, and why an un-calibrated literature-range value
(commonly ~0.4-0.7 for this pair) should be flagged as a judgment call,
not presented as a fitted result. A useful way to present that
uncertainty honestly, absent a paired record: run the same case at a
few `tau` values spanning the reference range as separate, clearly-
labeled variants, rather than picking one number -- the same "run both/
several, report all, don't pick a winner" pattern already used for the
FFA bootstrap stress-test case.

### Empirical hydrograph ensemble mode

Optional for every case -- omit the row entirely for the default,
`"scaled"` mode (this case's own `inflow_hydrograph.csv` scaled by
`peak_scale`/`volume_scale`, everything above). To opt a case into
drawing whole hydrographs directly from a pre-generated ensemble file
instead:

```
mc_outer_hydrograph_source,ensemble,-
```

When set, `mc_outer_peak_source`/`mc_outer_volume_source`/
`mc_peak_volume_dependence` are NOT read or required at all for that
case -- see `Data/Template/case_config.py`'s
`mc_outer_hydrograph_source="ensemble"` comments for the required
`mc_outer_distribution()` `ensemble_csv` key and the expected wide-CSV
format, and the README's "Empirical hydrograph ensembles" section for
the full reasoning.

```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 300 --n-inner 100
```

If `--n-outer` exceeds the ensemble's own replicate count, the run
still proceeds -- Layer 3 caps `n_outer` at the available count and
prints a clear warning (check the console output/`mc_summary.txt`'s
"Draws:" line for the actual count used, which may be lower than what
you requested). Generate more replicates from the ensemble tool itself
if you need more outer draws than that.

### Gate-reliability tier switch

Only relevant for a case whose `case_config.py` defines
`mc_gate_availability()` (see `Data/Template/case_config.py`'s
commented example, and the README's "Gate-availability (failure-to-
open) risk" section, for how to add it to a case that doesn't have it
yet). The switch lives entirely in that case's own `scalars.csv` --
edit `mc_gate_reliability_source`, no command-line flag or environment
variable needed:

CI-evidence-based reliability:
```
mc_gate_reliability_source,ci,-
```
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

Flat rate instead:
```
mc_gate_reliability_source,flat,-
mc_gate_p_fail,0.05,-
```
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

Flat rate, custom `p_fail` -- just change the value in `scalars.csv`:
```
mc_gate_reliability_source,flat,-
mc_gate_p_fail,0.10,-
```
```bash
uv run python Module/mc_layer3.py <CaseName> --n-outer 100 --n-inner 150
```

`Module/mc_layer3.py` itself never reads an environment variable for
this -- the switch is entirely a `scalars.csv` value, read via
`validate_mc_sources()` and passed to that case's own
`mc_gate_availability(source, p_fail)`.

## Setup (once, before any of the above)

```bash
uv sync
```

## Git (this repo specifically)

`master` (local) tracks `origin/Master` (remote) -- a case-sensitivity
mismatch this repo has hit before. `push.default` should already be set
to `upstream` (see chat history) so plain `git push`/`git pull` work
correctly; if a fresh clone ever shows the old "branch names don't
match" error again:
```bash
git config push.default upstream
```
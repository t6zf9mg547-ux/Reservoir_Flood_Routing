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
exactly what's missing. See the README's "Layer 3" section and
`Data/Template/scalars.csv`'s `options` column for the valid values and
what each one needs.

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
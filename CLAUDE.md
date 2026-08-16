# CLAUDE.md

Guidance for Claude Code when working in this project.

## What this project is

A Python flood-routing / dam-safety toolkit, built in layers:

- **Layer 1** (`Module/run_case.py`) — single-scenario routing kernel:
  given an inflow hydrograph, a stage-storage curve, and an operating
  rule, routes the flood and reports the peak reservoir level.
- **Layer 2** (`Module/optimize_layer2.py`) — NSGA-II optimization of a
  case's gate operating rule (via pymoo).
- **Layer 3** (`Module/mc_layer3.py`, `Module/mc_layer3_outer_only.py`)
  — two-loop Monte Carlo uncertainty analysis. Outer loop: climate/
  flood-frequency uncertainty (peak/volume scale factors, optionally
  correlated via a copula, or a real hydrograph ensemble). Inner loop:
  physical/rating-coefficient and gate-reliability uncertainty.

Full detail on all three lives in `README.md` — this file is
orientation and working conventions, not a substitute for reading it.
`RUN_COMMANDS.md` has the exact command for every standard operation.

## Package management

- **uv only.** Never use `pip`, `pip install`, `conda`, or
  `python -m venv` directly.
  - Add a dependency: `uv add <package>`
  - Remove a dependency: `uv remove <package>`
  - Install/sync environment: `uv sync`
  - Run any script: `uv run python Module/your_script.py`
- Never edit `.venv` by hand and never hand-edit dependency versions in
  `pyproject.toml` — use `uv add`/`uv remove` so the lockfile stays
  consistent.
- `[tool.uv] package = false` is intentional (scripts project, not an
  installable package) — don't remove it.

## Project structure

- `Module/` — all Python source. Generic engine code (`mc_layer3.py`,
  `outlets.py`, `reservoir.py`, `hydrograph.py`, `solver.py`,
  `reliability/`) contains ZERO case-specific numbers or references —
  see "The generic-engine / case-specific-data split" below.
- `Data/<CaseName>/` — one folder per dam/case: `case_config.py`,
  `scalars.csv`, `reservoir_curve.csv`, `inflow_hydrograph.csv`,
  reference reports/rating tables. Not tracked in git.
- `Output/<CaseName>/<Layer>/`, `Plot/<CaseName>/<Layer>/` — generated
  results/figures, one subfolder per layer (`Baseline`, `Optimised`,
  `MonteCarlo`, `OuterOnly`). Not tracked in git.
- `NEW_CASE_CHECKLIST.md` — step-by-step guide for onboarding a new
  dam/case. Start here, and start with its "gather the project's main
  features and constraints" section, before writing any code.

## The generic-engine / case-specific-data split

This is the core architectural rule of the project: `Module/` code
must work for ANY case without modification. Case-specific physics
(gate dimensions, discharge coefficients, operating rules, calibration
values) live entirely in that case's own `Data/<CaseName>/case_config.py`
— never hardcoded into `Module/`. If you find yourself wanting to
special-case something in `Module/` for a specific dam, that's a sign
the abstraction is wrong, not that the special case is justified.

## `case_config.py` conventions

- `build_outlets(rule_overrides=None, physical_overrides=None, ...)` —
  returns the list of Outlet objects. Read `Data/Template/case_config.py`'s
  extensively-commented version before writing a new one; propose the
  design before implementing (see "Working style" below).
- **Keep a gate's physical leaf height and its maximum OPERATIONAL
  opening as two separately-named values, even when they're set
  equal.** These are conceptually different things (a physical
  dimension vs. an operating-rule limit) that are easy to silently
  conflate into one reused number — this was a real bug, found and
  fixed across every case in this project (see `README.md`'s
  `on_off=False` vs. stuck-gate section for the full reasoning). Use a
  named constant like `GATE_HEIGHT` for the physical dimension,
  distinct from `a_max`.
- `mc_uncertain_params()`, `mc_outer_distribution()`,
  `mc_gate_availability(source, p_fail=None)` — all OPT-IN hooks for
  Layer 3. A case with none of them still runs Layer 3 correctly, just
  with less uncertainty modeled.
- `scalars(case_dir)` — reads `scalars.csv`, must return
  `dict[str, float | str]` (mixed types — Layer 3's `mc_*_source` rows
  are text, everything else is numeric). Copy this function verbatim
  from an existing case rather than rewriting it.

## Layer 3's explicit-source contract

Every case running Layer 3 MUST declare, in `scalars.csv`:
`mc_outer_peak_source`, `mc_outer_volume_source`, and (only if
`mc_gate_availability()` is defined) `mc_gate_reliability_source` —
`validate_mc_sources()` hard-errors before any simulation runs if
these are missing or don't match what's actually available. There is
no silent default. See `README.md`'s "Layer 3" section for the full
contract, including the optional `mc_peak_volume_dependence` and
`mc_outer_hydrograph_source=ensemble` extensions.

## Conventions

- Python `>=3.13`.
- `tkinter` is stdlib — never add it as a dependency.
- Keep `Data/`, `Output/`, `Plot/` out of git (already in
  `.gitignore`); don't change that without asking.
- No changelog for this project (explicitly not wanted) -- don't
  create or reference a `CHANGELOG.md`.

## General rules

- Don't add error handling, abstractions, or config options for cases
  that can't happen — keep scripts direct and readable.
- Don't introduce new top-level folders or restructure the layout
  without asking first.
- Prefer editing existing scripts in `Module/` over creating new ones
  for small changes. When a new script genuinely is warranted (e.g.
  `mc_layer3_outer_only.py`), factor out logic it shares with an
  existing script into a common, importable function rather than
  duplicating it.

## Working style

- **Confirm the design before implementing, for anything non-trivial.**
  When there's more than one reasonable way to build a feature, or a
  design choice has real consequences (what a default should be, how
  an edge case should behave, what gets validated vs. silently
  allowed), propose the approach and wait for confirmation before
  writing code. Don't let an autonomous multi-step plan run through
  several files on an unconfirmed assumption.
- **Verify by actually running the code, not just by reading it or
  syntax-checking it.** Run the real script against real (or realistic
  test) data before calling a change complete. A file that only passes
  `python -m py_compile`/linting has not been verified. If the real
  script can't be run end-to-end (missing data, missing external
  service), say so explicitly and test as much of the real logic as
  possible instead of skipping verification silently.
- **When two scripts need the same logic, factor it out into one
  shared, importable function or module — don't duplicate it.**
  Duplicated logic drifts the first time one copy gets fixed and the
  other doesn't.
- **When asked to review or audit code, actually read every relevant
  file.** Don't just grep for the one issue already mentioned — the
  most useful findings are usually the ones nobody knew to look for
  yet.
- **Don't rely on memory of a file's contents from earlier in the
  session — re-read it before making claims about it or building on
  it.** This matters most in long sessions, where a file may have
  changed since it was last read, or memory of its exact contents may
  be imprecise.
- **Git is not a silent "finishing step."** You have real commit/push
  access in this environment — use it deliberately, not automatically.
  Stage changes and propose a commit message, but wait for explicit
  confirmation before committing, and confirm again before pushing.
  Keep commit messages concise unless asked for more detail. This
  repo's local branch is `master`, tracking remote `Master`
  (case-sensitivity mismatch already resolved via
  `git config push.default upstream` — don't need to touch this again
  unless a fresh clone shows the old error).

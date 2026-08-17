---
name: new-case
description: Guides onboarding a new dam or reservoir case into this project, covering feature gathering and setup across all three layers (routing kernel, gate operating-rule optimization, Monte Carlo uncertainty analysis). Use when the user asks to add a new dam, case, or reservoir, wants to create a case_config.py for a new site, or asks how to set up a new case in this project.
---

Read `NEW_CASE_CHECKLIST.md` at the project root in full before doing
anything else -- it is the single source of truth for this workflow,
kept as a plain file (not duplicated here) so it stays readable
directly on GitHub and never drifts out of sync with this trigger.

Follow it in order:

1. Section 0 (gather the project's main features and constraints)
   FIRST, before writing any code. Every item should be a real, sourced
   number or an explicitly labeled "TBD" -- never a guess standing in
   for a real value without saying so.
2. Layer 1, then Layer 2, then Layer 3, in that order -- each layer's
   checklist ends with actually running the real script against real
   data and inspecting the output, not just writing the file.

Per this project's `CLAUDE.md`:
- Propose the `case_config.py` design and wait for confirmation before
  writing it.
- Verify every claim of "done" by actually running the relevant script,
  not by reading or syntax-checking alone.
- Stage git changes and propose a commit message; wait for explicit
  confirmation before committing or pushing.

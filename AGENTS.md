# Project delivery instructions

## The outcome that matters

Make the conventional single-driver workflow useful: an operating band and optional constraints should produce a suitable experimental driver/horn shortlist, horn dimensions, acoustic STEP geometry, response/impedance predictions and a clear reproducible report. MEH belongs in https://github.com/timini/meh-design-studio.

Prioritize a working end-to-end example, meaningful acoustic correctness, usable outputs and clear limitations. Deliver coherent improvements instead of expanding the validation infrastructure indefinitely. Keep physical assembly qualification separate from delivery of a usable experimental design tool; disclose missing evidence without pretending it has passed.

## One review round

1. Complete the proposed change and run the relevant checks.
2. Obtain **one round of review feedback** for the change. Review suggestions are advisory and must be assessed against the actual product goal.
3. Fix substantive findings: incorrect supported behavior or predictions, broken required checks, security/data-loss risks, or materially misleading claims. Defer cosmetic changes, speculative edge cases, optional refactors and extra audit machinery that do not affect the intended use.
4. Check the fixes with targeted tests, let required CI finish, then merge when authorized. **Do not request another full review after every fix or demand a fresh automated review of every final commit.**

Additional review requires an explicit user request or a specific newly introduced critical risk that cannot be resolved through focused inspection/testing. Keep any such review confined to that risk; it must not restart the whole review cycle. A new unrelated feature may receive its own single round.

If an automated reviewer fails, retry at most once, then inspect the change directly, record the service failure and proceed when the relevant checks pass and merge is authorized. Do not mistake a service error for a code defect or a successful review. Do not claim the reviewer uses the latest model unless its identity is actually available.

## Proportionate validation and follow-through

- Agree a concrete acceptance result from the task context and work toward it. Avoid moving the completion criteria as new optional ideas arise.
- Prefer checks of user-visible behavior and acoustic results over tests that merely mirror implementation details.
- Reuse existing trustworthy evidence. Repeat expensive solves only when changed solver code, geometry, inputs or unresolved failures justify them. A comparison/reporting change can be checked against preserved raw results without rerunning unchanged physics.
- Once relevant checks pass, stop broadening or repeating tests without a specific reason. Do not weaken numerical tolerances or hide failed results to obtain a pass.
- Keep research focused on a concrete missing input. When repeated searches reveal the same missing calibration/source/geometry data, record that dependency and the next actionable acquisition step instead of collecting more unusable curves indefinitely.
- Preserve user work and historical results. Use an isolated checkout when the working tree contains unrelated edits.
- Finish the delivery: push, merge when authorized and checks pass, then summarize what a user can do now and the few remaining product blockers. Do not substitute review activity or documentation volume for progress toward designing a horn.

See `CLAUDE.md` for repository commands and `docs/SINGLE_HORN_ROADMAP.md` for the current technical scope. These delivery and review rules apply to all coding agents working in this repository.

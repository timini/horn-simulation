# Documentation

Start with the [README](../README.md) for installation and a first run.

## Using the pipeline

| Page | What it covers |
| --- | --- |
| [USAGE](USAGE.md) | Every command, launcher option and pipeline parameter, with defaults. Run folders, resume and lookup. |
| [OUTPUTS](OUTPUTS.md) | Every output file and column, and how to read f3 and ripple. |
| [SEARCH_AND_SCORING](SEARCH_AND_SCORING.md) | How auto mode chooses, rejects and ranks designs. The score formula. |
| [DRIVERS](DRIVERS.md) | Driver databases, record format, adding and validating drivers, front and rear chambers, evidence gaps. |
| [WORKED_EXAMPLE](WORKED_EXAMPLE.md) | An 800 to 1600 Hz run from command to result. |
| [Examples](../examples/README.md) | Preserved example outputs you can open without running anything. |
| [Example site](https://timini.github.io/horn-simulation/) | The examples on GitHub Pages. |

## How it works

| Page | What it covers |
| --- | --- |
| [PHYSICS](PHYSICS.md) | What the solver, radiation models, driver coupling and screen compute, and their limits. |
| [ACOUSTIC_CONTRACT](ACOUSTIC_CONTRACT.md) | Sign, unit and phasor conventions in detail. |
| [MODAL_APERTURE](MODAL_APERTURE.md) | The default `modal_baffled` mouth model. |

## Evidence and status

| Page | What it covers |
| --- | --- |
| [STATUS](STATUS.md) | What works, what is experimental, open issues and remaining work. |
| [VALIDATION](VALIDATION.md) | Index of every validation record, script and CI job, with results. |
| [DESCOPED](DESCOPED.md) | Items from the original brief that were not built, and why. |

Validation write-ups, linked from [VALIDATION](VALIDATION.md):

- [INDEPENDENT_SOLVER_VALIDATION](INDEPENDENT_SOLVER_VALIDATION.md): interior FEM, motor and cavity checks against independent code.
- [EXTERIOR_RADIATION_VALIDATION](EXTERIOR_RADIATION_VALIDATION.md): mouth radiation against the MMM toolbox, and finite-baffle pilots.
- [POST_HIXSON_VALIDATION](POST_HIXSON_VALIDATION.md): measured exponential horn input impedance.
- [LOSS_PHYSICS_VALIDATION](LOSS_PHYSICS_VALIDATION.md): wall losses and finite-flange radiation against measured and published pipes.
- [CROSS_SOFTWARE_VALIDATION](CROSS_SOFTWARE_VALIDATION.md): openwind, bempp-cl and Hornresp.
- [CANDIDATE_RESOLUTION](CANDIDATE_RESOLUTION.md) and [MODAL_CANDIDATE_RESOLUTION](MODAL_CANDIDATE_RESOLUTION.md): mesh, loft and frequency convergence of generated candidates.
- [VALIDATION_DATA](VALIDATION_DATA.md): measurement sources, the importer, and the Hajicek and Neely horn standard.
- [PHYSICAL_REFERENCE_AUDIT](PHYSICAL_REFERENCE_AUDIT.md): candidate physical reference assemblies.

## Contributing

| Page | What it covers |
| --- | --- |
| [DEVELOPMENT](DEVELOPMENT.md) | Repository layout, environments, Docker images, tests, CI and conventions. |
| [AGENTS.md](../AGENTS.md) | Delivery rules for coding agents and contributors. |
| [history/](history/README.md) | Dated plans, reports and research kept for provenance. |

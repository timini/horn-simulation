# Descoped items

The original brief ([design document](<history/A Comprehensive Design Document for an Open-Source, Cloud-Deployable Horn Loudspeaker Simulation and Optimization Pipeline.md>) and its [revision](<history/Revised Horn Simulation Pipeline Architecture.md>), 2025) asked for more than was built. This page records each item that was not delivered, or was delivered differently, and why. Items the brief asked for that are delivered (SPL, electrical impedance, excursion and group delay curves; weighted scoring; parametric profiles; STEP export; Docker containers; parallel simulation) are described in the rest of the documentation.

## Not built

| Brief item | Decision | Reason |
| --- | --- | --- |
| Cloud deployment on Google Cloud: Cloud Workflows, Cloud Batch on Compute Engine, Cloud Run, Cloud Storage, Cloud SQL | Not built. Runs locally with Docker and Nextflow. | The tool is used by one person on one machine. A search of the default size finishes in minutes locally. Nextflow can target cloud executors later without changing the pipeline, if a real need appears. |
| PostgreSQL driver database | Not built. One JSON file per driver in git. | A server adds setup and gives a single local user nothing. Git history is the audit trail, and the files are easy to review in a pull request. |
| PDF datasheet parser (pdfplumber, regex) with a human-review form (Flask or Streamlit) | Not built. `horn-add-driver` takes datasheet values with explicit units. | Datasheet layouts vary too much for reliable extraction, and every value still needs a person to check it. Typing the values with units, then running the loader and validator, is that check. |
| Horn walls, wall thickness, enclosure and driver chamber solids in the CAD | Not built. STEP files are the air volume only. | The acoustic model needs only the air. Walls, flanges and adapters depend on the build method and the driver, and need mechanical design that the pipeline cannot verify. |
| STL export | Not built. | STEP is the exchange format CAD tools import. Any CAD tool can export STL from it; the report has its own render. |
| Rectangular and other non-circular throats and mouths | Not built. Axisymmetric circular horns only. | The radiation model, meshing checks and validation evidence all assume an axisymmetric volume. Rectangular horns need a different exterior model; see issues [#49](https://github.com/timini/horn-simulation/issues/49) and [#51](https://github.com/timini/horn-simulation/issues/51). |
| Hornresp "rapid path" automation | Not built. Replaced by an in-house Webster transfer-matrix screen. A Hornresp input kit is provided for a manual comparison. | Hornresp is Windows-only freeware with no scripting interface. The revised brief already dropped it. The screen is compared with openwind instead. |
| Coupled FEM-BEM exterior solve with bempp-cl | Built, then removed from the solver. | It coupled the whole horn boundary instead of the mouth, which is wrong for a horn in a baffle. `modal_baffled` replaced it for the infinite-baffle case. bempp-cl is now used only as an independent check. Free-standing horns and finite baffles are not modelled. |
| Driver cost (`price_usd`) in the score | Not built. | No reliable price data with provenance. The score has sensitivity, ripple, extension and size terms instead. |
| Formal multi-objective optimisation (pymoo, Pareto fronts) | Not built. A weighted score with fixed exchange rates and a bounded refinement. | A weighted score with stated exchange rates is easier to read and audit. The shortlist audit shows the screen keeps the best FEM designs on the tested grid. |
| Hyperbolic flare constant T as a pipeline input | Python API and `horn-geometry --hypex-t` only. The pipeline uses T = 0. | The screen, report and candidate CSV do not carry T, so allowing it in the pipeline would mix geometries. |
| Machine-learning surrogate models | Not built. | Listed as a future idea in the brief. The Webster screen is already fast enough for the grid sizes used. |
| Crossover design and multi-way system integration | Not built. | Out of scope for a single-horn tool. Multiple-entry horns and their crossovers are in [MEH Design Studio](https://github.com/timini/meh-design-studio). |
| Structural vibroacoustic coupling (wall vibration) | Not built. | Listed as a future idea. Walls are rigid. |
| Web front end | Not built. Self-contained HTML reports and a static GitHub Pages site of examples. | One user, local runs. The reports open in any browser. |

## Delivered differently

| Brief item | What was built |
| --- | --- |
| FreeCAD headless geometry generation | gmsh with the OpenCASCADE kernel, which the mesher already needs. See [history/GEOMETRY_KERNEL_EVALUATION](history/GEOMETRY_KERNEL_EVALUATION.md). |
| A central Python orchestrator | Nextflow (`main.nf`) with a Python launcher for provenance and resume. |
| `scoring.json` weights file | `--score_weights`, or `score_weights` in a `-params-file`. |
| Driver velocity applied at the throat inside the FEM solve | The FEM solves a unit-pressure transfer once per horn, and each driver is coupled through a lumped model. For a linear system this gives the same answer for every driver at the cost of one solve. The solver still has a velocity-inlet mode, used in validation. |
| Directivity on a grid of points | Experimental far-field directivity from the mouth modes, infinite baffle only, not compared with measurement. |
| Pressure at a far-field microphone from BEM | On-axis pressure at `observation_distance` from a Rayleigh integral over the mouth, infinite baffle. |

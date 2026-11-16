# Run locations, lookup and resume

The supported entry points are `just run`, `just run-auto` and `just run-fullauto`. Each fresh launch creates a timestamp/UUID folder beneath `results/`. Reports and data are under that run's `outputs/`; its manifest, logs, work cache and source snapshot are alongside them. Two fresh launches preserve separate directories.

Use a chosen directory name when useful:

```sh
python3 scripts/run_pipeline.py --run-dir results/my-500hz-horn --mode auto --target_f_low 500 --target_f_high 4000 --drivers_db data/drivers-curated
```

An existing directory is rejected for a fresh launch. Resume with the same directory:

```sh
python3 scripts/run_pipeline.py --run-dir results/my-500hz-horn -resume
```

Resume checks the original source, driver/input files, containers and Nextflow engine, and uses the original session and specification. Changed inputs require a new run. Preserve `work/` and the Nextflow cache while you need resume. A directory name is the supported run label; no separate label option is needed.

## Parameter files and images

A `-params-file` (JSON or YAML) is forwarded to Nextflow. The manifest records its sha256 under `input_sha256`, together with the driver database, STEP file or render it names; a value given on the command line takes precedence, as in Nextflow.

`just build` labels each image with `horn.source_sha256`, a hash of the files the image is built from (`scripts/image_source.py`). The launcher recomputes the hash and prints a warning when an image is unlabelled or was built from other sources, and records both values under `image_sources` in the manifest. Rebuild with `just build` after changing package code. `--require-current-images` turns the warning into an error; CI uses it. `--image-tag TAG` (or `HORN_IMAGE_TAG`) runs `horn-*:TAG` images built with `just build TAG`.

Images install third-party Python packages from `uv.lock` and pin their base images by digest. The dolfinx-based images keep the base image's numpy (1.26.4), mpi4py and gmsh, because dolfinx, petsc4py and the system gmsh are compiled against them; `bempp-cl` is pinned to 0.3.1 in the solver Dockerfile. After changing dependencies, run `uv lock` and commit the lockfile; CI runs `uv lock --check`.

The pipeline requires Nextflow 24.10.x (`manifest.nextflowVersion`). `just test-nextflow` runs nf-test 0.9.5 through `scripts/nf_test.sh`, which downloads it to `.tools/` and checks its sha256 on first use; set `NF_TEST` to use another copy.

## Find a run

```sh
just latest-run
just runs
just latest-run --status running
just latest-run --root /path/to/custom/run-parent
```

`latest-run` prints the absolute directory of the most recently **completed** run, ordered by its recorded finish time. It excludes running, failed, malformed and historical folders. An unsuccessful search exits with a diagnostic. `--status failed`, `--status running` or `--status any` explicitly changes the selection; running runs use their start time. Repeat `--root` to search multiple parents. The commands also work directly as `python3 scripts/runs.py latest` and `python3 scripts/runs.py list` without `just`.

Before the first run, `runs` returns an empty JSON list and `latest-run` reports that no completed run exists. An explicitly supplied nonexistent `--root` is an error.

`runs` returns a read-only JSON inventory, including folders without a run manifest as `unmanaged`, and malformed manifests as `invalid_manifest`. It reads only immediate child directories, does not follow directory symlinks, and does not modify results. A recorded running state may be stale after a machine crash; the locator does not infer liveness. Completed means the workflow exited successfully: it may legitimately report no feasible design or insufficient evidence. It does not mean a physically validated horn was found.

## Direct Nextflow use

Advanced direct invocation requires an explicit `--outdir`. There is no shared default output location:

```sh
nextflow run main.nf -profile docker --outdir results/manual-example/outputs
```

The caller must choose a fresh output directory and preserve the matching work/cache/session for direct resume. Direct invocation does not provide the launcher's checked provenance or a launcher manifest and will therefore not appear as a completed run in the locator. Use the launcher for normal design work.

## Historical folders and opt-in cleanup

Run `just runs` before cleanup. Older results, validation directories and folders without manifests remain visible as unmanaged; their names and filesystem dates do not establish whether they are safe to remove. Inspect reports, source/inputs, raw response files, logs and validation provenance individually. Keep the original and failed reference evidence used by published validation reports.

Cleanup is deliberately manual: archive the complete selected run outside the repository, verify the archive can be read, then remove only that specific run if its outputs and resume cache are no longer needed. The project does not delete, move or relabel historical folders automatically. Never clear all of `results/` or `work/` as a shortcut to making space: those locations may contain the only copy of a validation result or the cache needed to resume a run.

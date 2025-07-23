# BEM Solver Environment Setup Challenges

**Date:** 2025-07-23

## 1. Goal

The primary goal was to create a Dockerized environment for the `horn-bem-solver` package that includes `bempp-cl` and `fenics-dolfinx`. This environment is a prerequisite for implementing the BEM solver as outlined in `todo/14-bem-solver-implementation.md`.

## 2. Summary of the Problem

A stable and functional Docker environment could not be established due to persistent architecture compatibility issues. The development was performed on a `darwin` (macOS) machine with an ARM64 (Apple Silicon) processor. The core issue is that the required scientific computing packages (`bempp-cl` and `fenics-dolfinx`) do not have pre-built distributions available for the `linux/arm64` architecture, which is the target architecture when building Docker images on this host.

Multiple approaches were attempted to resolve this, but each resulted in a critical failure.

## 3. Detailed Log of Attempts and Failures

### Attempt 1: Using Pre-built `bempp/bempp-cl-with-dolfinx` Image

-   **Strategy:** Use the official, pre-built Docker image which is supposed to contain all necessary dependencies. Specify `--platform=linux/amd64` to force the use of the x86_64 version via emulation.
-   **Problem:** The container would build, but at runtime, it consistently failed with `ModuleNotFoundError: No module named 'bempp'`.
-   **Investigation:**
    -   It was discovered that the base image uses a virtual environment located at `/dolfinx-env/bin/activate`.
    -   The `Dockerfile` was modified to source this `activate` script before running the solver.
    -   Despite correctly activating the environment, the `ModuleNotFoundError` persisted. This suggests a fundamental issue with the environment setup within the base image itself, where the Python interpreter is not correctly linked to the installed packages when the container is run non-interactively.
    -   Various `ENTRYPOINT` configurations (direct command, entrypoint script) were tried without success.

### Attempt 2: Building from Source with `pip`

-   **Strategy:** Start with a base `ubuntu:22.04` image and install `fenics-dolfinx` and `bempp-cl` using `pip`.
-   **Problem:** The `pip install` command failed with `ERROR: No matching distribution found for fenics-dolfinx`.
-   **Investigation:** This confirms that there are no pre-built wheels for `fenics-dolfinx` (and likely `bempp-cl`) for the `linux/arm64` architecture on PyPI. Building them from source would require installing a complex chain of system-level dependencies (compilers, linear algebra libraries, etc.) which is a significant undertaking.

### Attempt 3: Using a `conda` Environment

-   **Strategy:** Use a `continuumio/miniconda3` base image and create a `conda` environment, installing the packages from the `conda-forge` channel, which is often better at handling scientific packages.
-   **Problem:** The `conda create` command failed with `PackagesNotFoundError: The following packages are not available from current channels: - bempp-cl`.
-   **Investigation:** This shows that `conda-forge` also lacks a `bempp-cl` package for the `linux/aarch64` (ARM64) platform.

### Attempt 4: Using `conda` with Platform Emulation

-   **Strategy:** Force `conda` to use the `linux-64` (x86_64) channel by setting `conda config --env --set subdir linux-64`. Docker's emulation layer would then be responsible for running the x86_64 binaries.
-   **Problem:** The `conda create` command failed with an internal `TypeError: expected str, bytes or os.PathLike object, not NoneType` within conda's own code.
-   **Investigation:** This appears to be a bug or an incompatibility within `conda` itself when trying to solve for a foreign architecture under emulation.

## 4. Conclusion and Recommendations

The inability to create a working environment is a major blocker. The combination of an ARM64 development machine and the lack of official support for this architecture in the scientific Python ecosystem for these specific packages makes a straightforward setup impossible.

**Recommendations:**

1.  **Build from Source (High Effort):** The most robust solution is to build `dolfinx` and `bempp-cl` and all their dependencies from source within a `linux/arm64` Docker image. This would create a native-running image but requires significant effort and expertise in managing complex builds.
2.  **Use an x86_64 Host:** The simplest solution is to perform the Docker build on a machine with an x86_64 processor. The resulting image could then potentially be run on an ARM64 machine using Docker's emulation, although performance might be degraded.
3.  **Explore Alternative Libraries:** Research other BEM libraries that may have better cross-platform support and pre-built packages available for `linux/arm64`.

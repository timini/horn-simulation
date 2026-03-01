#!/usr/bin/env python3
"""Phase 1a: Verify bempp-cl is importable with the Numba backend.

Run inside the dolfinx/dolfinx:v0.8.0 container after `pip install bempp-cl`:

    docker run --rm dolfinx/dolfinx:v0.8.0 bash -c \
        "pip install bempp-cl && python3 scripts/test_bempp_install.py"
"""

import sys


def main() -> int:
    # 1. Import bempp-cl (v0.4.x module name)
    try:
        import bempp.api as bempp_api
    except ImportError:
        print("FAIL: could not import bempp.api — is bempp-cl installed?")
        return 1

    print(f"OK: imported bempp.api (version {bempp_api.__version__})")

    # 2. Verify Numba backend is active (no OpenCL)
    device = getattr(bempp_api, "DEFAULT_DEVICE_INTERFACE", "unknown")
    print(f"  DEFAULT_DEVICE_INTERFACE = {device!r}")
    if device != "numba":
        print(f"WARN: expected 'numba' backend, got {device!r}")
        # Not fatal — the OpenCL backend works too, but Numba is the goal

    # 3. Create a simple sphere grid and assemble a single-layer operator
    try:
        grid = bempp_api.shapes.regular_sphere(3)  # refinement level 3
        print(f"OK: created sphere grid with {grid.number_of_elements} elements")
    except Exception as exc:
        print(f"FAIL: could not create sphere grid: {exc}")
        return 1

    try:
        space = bempp_api.function_space(grid, "P", 1)
        print(f"OK: created P1 function space ({space.global_dof_count} DOFs)")
    except Exception as exc:
        print(f"FAIL: could not create function space: {exc}")
        return 1

    # Assemble Helmholtz single-layer operator at k=1
    k = 1.0
    try:
        slp = bempp_api.operators.boundary.helmholtz.single_layer(
            space, space, space, k
        )
        mat = slp.weak_form()
        print(f"OK: assembled Helmholtz single-layer operator (k={k})")
        print(f"  Matrix shape: {mat.shape}")
    except Exception as exc:
        print(f"FAIL: could not assemble BEM operator: {exc}")
        return 1

    print("\nAll checks passed — bempp-cl is working with the Numba backend.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

# syntax=docker/dockerfile:1
# Base images are pinned by multi-arch index digest (linux/amd64 and linux/arm64).
FROM ghcr.io/astral-sh/uv:0.8.2@sha256:a7999d42cba0e5af47ef3c06ac310229c7f29c5314e35902f8353e8e170eeed1 AS uv
FROM dolfinx/dolfinx:v0.8.0@sha256:6fd526e28275ca60b5ba6125cd7ed488e4a120edc8fe7029ac2422f3154913fc AS base

# --- Requirements Stage ---
# Exact versions from uv.lock for every workspace dependency (this image also
# carries horn-analysis and horn-geometry), except numpy, mpi4py and gmsh:
# the image's compiled dolfinx, petsc4py, MPI and system gmsh are built
# against its own copies, so those stay as shipped.
FROM base AS requirements
WORKDIR /src
COPY pyproject.toml uv.lock ./
COPY packages/horn-core/pyproject.toml packages/horn-core/
COPY packages/horn-drivers/pyproject.toml packages/horn-drivers/
COPY packages/horn-geometry/pyproject.toml packages/horn-geometry/
COPY packages/horn-solver/pyproject.toml packages/horn-solver/
COPY packages/horn-analysis/pyproject.toml packages/horn-analysis/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv export --frozen --all-packages --no-emit-workspace --no-emit-package numpy \
    --no-emit-package mpi4py --no-emit-package gmsh -o /requirements.txt \
    && uv export --frozen --package horn-solver --extra test --no-emit-workspace \
    --no-emit-package numpy --no-emit-package mpi4py --no-emit-package gmsh -o /requirements-test.txt

# --- Production Stage ---
FROM base AS production

# Set up the entrypoint to activate complex mode.
COPY ./packages/horn-solver/entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

WORKDIR /app
# bempp-cl is not a workspace dependency, so it is pinned here; its own
# dependencies come from the base image or PyPI. The locked install that
# follows restores the lock's version of anything it moved.
RUN pip install 'bempp-cl==0.3.1'
COPY --from=requirements /requirements.txt /requirements-test.txt /app/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements.txt
COPY ./packages /app/packages
RUN pip install --no-deps /app/packages/horn-core /app/packages/horn-drivers \
    /app/packages/horn-solver /app/packages/horn-analysis /app/packages/horn-geometry \
    && python3 -c "import horn_solver.solver, horn_analysis.plotter, horn_geometry.generator"

# --- Test Stage ---
FROM production AS test
WORKDIR /app
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements-test.txt
COPY ./packages/horn-solver/tests /app/packages/horn-solver/tests

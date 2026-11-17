# dolfinx supplies native gmsh, including on Apple Silicon / Linux arm64.
# Base images are pinned by multi-arch index digest (linux/amd64 and linux/arm64).
FROM ghcr.io/astral-sh/uv:0.8.2@sha256:a7999d42cba0e5af47ef3c06ac310229c7f29c5314e35902f8353e8e170eeed1 AS uv
FROM dolfinx/dolfinx:v0.8.0@sha256:6fd526e28275ca60b5ba6125cd7ed488e4a120edc8fe7029ac2422f3154913fc AS base

# --- Requirements Stage ---
# Exact versions from uv.lock, except numpy and gmsh: the image's compiled
# dolfinx, petsc4py and system gmsh are built against its own copies, so
# those stay as shipped (numpy 1.26.4, gmsh 4.12.2).
FROM base AS requirements
WORKDIR /src
COPY pyproject.toml uv.lock ./
COPY packages/horn-core/pyproject.toml packages/horn-core/
COPY packages/horn-drivers/pyproject.toml packages/horn-drivers/
COPY packages/horn-geometry/pyproject.toml packages/horn-geometry/
COPY packages/horn-solver/pyproject.toml packages/horn-solver/
COPY packages/horn-analysis/pyproject.toml packages/horn-analysis/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv export --frozen --package horn-geometry --no-emit-workspace \
    --no-emit-package numpy --no-emit-package gmsh -o /requirements.txt \
    && uv export --frozen --package horn-geometry --extra test --no-emit-workspace \
    --no-emit-package numpy --no-emit-package gmsh -o /requirements-test.txt

FROM base AS production
WORKDIR /app
COPY --from=requirements /requirements.txt /requirements-test.txt /app/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements.txt
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-geometry /app/packages/horn-geometry
# Do not replace the system gmsh with an unavailable arm64 PyPI wheel.
RUN pip install --no-deps /app/packages/horn-core /app/packages/horn-geometry \
    && python3 -c "import gmsh, numpy; import horn_geometry.generator"

FROM production AS test
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements-test.txt
COPY ./packages/horn-geometry/tests /app/packages/horn-geometry/tests

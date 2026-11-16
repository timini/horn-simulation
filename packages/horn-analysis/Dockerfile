# syntax=docker/dockerfile:1
# Base images are pinned by multi-arch index digest, so linux/amd64 and
# linux/arm64 builds resolve the same release.
FROM ghcr.io/astral-sh/uv:0.8.2@sha256:a7999d42cba0e5af47ef3c06ac310229c7f29c5314e35902f8353e8e170eeed1 AS uv
FROM python:3.10.22-slim-bookworm@sha256:5be5aaecb962a41567bc8041f1e9e3432cb085b037fe057f9a562630b2e4b370 AS base

# --- Requirements Stage ---
# Export exact third-party versions from uv.lock. --frozen reads the lock
# without re-resolving; CI checks separately that the lock is current.
FROM base AS requirements
WORKDIR /src
COPY pyproject.toml uv.lock ./
COPY packages/horn-core/pyproject.toml packages/horn-core/
COPY packages/horn-drivers/pyproject.toml packages/horn-drivers/
COPY packages/horn-geometry/pyproject.toml packages/horn-geometry/
COPY packages/horn-solver/pyproject.toml packages/horn-solver/
COPY packages/horn-analysis/pyproject.toml packages/horn-analysis/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv export --frozen --package horn-analysis --no-emit-workspace -o /requirements.txt \
    && uv export --frozen --package horn-analysis --extra test --no-emit-workspace -o /requirements-test.txt

# --- Production Stage ---
FROM base AS production
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY --from=requirements /requirements.txt /requirements-test.txt /app/
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements.txt
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-drivers /app/packages/horn-drivers
COPY ./packages/horn-analysis /app/packages/horn-analysis
# Strip [tool.uv.sources] sections: uv rejects workspace refs outside a workspace.
# Workspace packages install without dependencies; the lock supplied those.
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    find /app/packages -name pyproject.toml \
    -exec sed -i '/\[tool\.uv\.sources\]/d; /workspace = true/d' {} + \
    && uv pip install --system --no-deps /app/packages/horn-core \
    /app/packages/horn-drivers /app/packages/horn-analysis \
    && python -c "import horn_analysis.plotter, horn_analysis.rank_pipeline"

# --- Test Stage ---
FROM production AS test
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*
RUN --mount=from=uv,source=/uv,target=/usr/local/bin/uv \
    uv pip install --system --no-deps --require-hashes -r /app/requirements-test.txt
COPY ./scripts/run_pipeline.py ./scripts/image_source.py /app/scripts/
COPY ./data/validation/boundary_lab_reference.json /app/data/validation/boundary_lab_reference.json
# The electrical output tests re-rank the preserved 800-1600 Hz example with a curated driver.
COPY ./data/drivers-curated /app/data/drivers-curated
COPY ./examples/800-1600-modal/ranking.json ./examples/800-1600-modal/horn-impedance.csv /app/examples/800-1600-modal/
# The chamber-policy tests screen a horn with the 6NMB420 example driver record.
COPY ./examples/6nmb420-320-5000/driver-database.json /app/examples/6nmb420-320-5000/
# The refinement-output test reads the refine process outputs from main.nf.
COPY ./main.nf /app/main.nf
WORKDIR /app
ENV PYTHONPATH=/app
COPY ./packages/horn-analysis/tests /app/packages/horn-analysis/tests

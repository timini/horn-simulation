# --- Base Image ---
# It is based on the dolfinx image because that image is guaranteed
# to have a working gmsh installation.
FROM dolfinx/dolfinx:v0.8.0 as base

# --- Builder Stage ---
# This stage installs our package along with its dependencies
FROM base as builder
RUN pip install uv
WORKDIR /app
COPY ./packages/horn-geometry/pyproject.toml .
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-geometry /app/packages/horn-geometry
RUN uv pip install --system /app/packages/horn-geometry

# --- Production Stage ---
# This is the final, clean image for the service.
FROM base as production
WORKDIR /app

# Copy the installed packages from the builder stage.
COPY --from=builder /usr/local/lib/python3.10/dist-packages /usr/local/lib/python3.10/dist-packages

# Copy the generator script, which is the entrypoint.
COPY --from=builder /app/packages/horn-geometry/src/horn_geometry/generator.py /app/generator.py

# Create and switch to a non-root user
RUN useradd --create-home --shell /bin/bash appuser
USER appuser
WORKDIR /home/appuser

# --- Test Stage ---
# This stage adds test files and dependencies to the production image.
FROM production as test
USER root
WORKDIR /app
ENV PYTHONPATH /app
RUN mkdir -p /app/.pytest_cache && chown -R appuser:appuser /app/.pytest_cache
RUN mkdir -p /app/.pytest_cache && chown -R appuser:appuser /app/.pytest_cache

# The production stage already contains the main dependencies.
# We just need to copy the source again to install the '[test]' extras.
COPY ./packages/horn-geometry/pyproject.toml .
COPY ./packages/horn-geometry /app/packages/horn-geometry
RUN pip install /app/packages/horn-geometry[test]

# Copy the tests themselves
COPY ./packages/horn-geometry/tests /app/packages/horn-geometry/tests

USER appuser
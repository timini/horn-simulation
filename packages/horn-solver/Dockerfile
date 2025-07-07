# --- Base Image ---
# Reverting to a specific, stable version of the official dolfinx image.
# This is a last resort, as other methods have failed.
FROM dolfinx/dolfinx:v0.8.0 as base

# --- Production Stage ---
FROM base as production

# The base image already contains a compatible FEniCSx environment.
# We only need to install uv.
RUN pip install uv

# Activate the complex-mode environment permanently.
# This is the crucial step that was missing.
ENV PKG_CONFIG_PATH=/usr/local/dolfinx-complex/lib/pkgconfig:$PKG_CONFIG_PATH
ENV CMAKE_PREFIX_PATH=/usr/local/dolfinx-complex/lib/cmake:$CMAKE_PREFIX_PATH
ENV PETSC_ARCH=linux-gnu-complex128-32
ENV PYTHONPATH=/usr/local/dolfinx-complex/lib/python3.10/dist-packages:$PYTHONPATH
ENV LD_LIBRARY_PATH=/usr/local/dolfinx-complex/lib:$LD_LIBRARY_PATH

WORKDIR /app
COPY pyproject.toml .
COPY ./packages /app/packages
RUN uv pip install --system /app/packages/horn-solver

# Create a non-root user to run the application
RUN useradd --create-home --shell /bin/bash appuser
USER appuser
WORKDIR /home/appuser


# --- Test Stage ---
FROM production as test
USER root
WORKDIR /app

# Copy test-specific files
COPY ./packages/horn-solver/tests /app/packages/horn-solver/tests

# Install test dependencies
RUN uv pip install --system /app/packages/horn-solver[test]

USER appuser
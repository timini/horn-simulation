# Reverting to a specific, stable version of the official dolfinx image.
# This is a last resort, as other methods have failed.
FROM dolfinx/dolfinx:v0.8.0 as base

# --- Production Stage ---
FROM base as production

# The base image already contains a compatible FEniCSx environment.
# Use pip (not uv) for bempp-cl to ensure it goes to Debian's dist-packages
# where Python's default sys.path can find it.
RUN pip install bempp-cl

# gmsh is pre-installed in the dolfinx image at /usr/local/lib (non-standard path
# without a RECORD file). Install a fresh copy to dist-packages so Python finds it.
RUN pip install --ignore-installed gmsh

# Set up the entrypoint to activate complex mode.
COPY ./packages/horn-solver/entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

WORKDIR /app
COPY ./packages/horn-solver/pyproject.toml .
COPY ./packages /app/packages
RUN pip install /app/packages/horn-solver

# --- Test Stage ---
FROM production as test
WORKDIR /app
# Copy test-specific files
COPY ./packages/horn-solver/tests /app/packages/horn-solver/tests

# Install test dependencies
RUN pip install /app/packages/horn-solver[test] /app/packages/horn-analysis

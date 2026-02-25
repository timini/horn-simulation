# Reverting to a specific, stable version of the official dolfinx image.
# This is a last resort, as other methods have failed.
FROM dolfinx/dolfinx:v0.8.0 as base

# --- Production Stage ---
FROM base as production

# The base image already contains a compatible FEniCSx environment.
# Install bempp-cl and gmsh to standard dist-packages.
RUN pip install bempp-cl
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

# Verify critical imports work at build time (no entrypoint, default PYTHONPATH)
RUN python3 -c "import bempp; print('bempp OK:', bempp.__file__)"
RUN python3 -c "import gmsh; print('gmsh OK:', gmsh.__file__)"

# Copy test-specific files
COPY ./packages/horn-solver/tests /app/packages/horn-solver/tests

# Install test dependencies
RUN pip install /app/packages/horn-solver[test] /app/packages/horn-analysis

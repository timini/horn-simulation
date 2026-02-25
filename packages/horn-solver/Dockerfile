# Reverting to a specific, stable version of the official dolfinx image.
# This is a last resort, as other methods have failed.
FROM dolfinx/dolfinx:v0.8.0 as base

# --- Production Stage ---
FROM base as production

# The base image already contains a compatible FEniCSx environment.
# bempp-cl 0.4+ renamed the module from 'bempp' to 'bempp_cl'; our code
# uses 'import bempp.api' so we need <0.4.
RUN pip install 'bempp-cl<0.4'
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

# --- Production Stage ---
FROM ubuntu:22.04 as production

# Set non-interactive frontend for package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    python-is-python3 \
    freecad \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install uv, the fast Python package installer
RUN pip install uv

# Create a virtual environment and set the path
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

WORKDIR /app
COPY pyproject.toml .
COPY ./packages /app/packages

# Install production dependencies
RUN uv pip install /app/packages/horn-geometry

# Remove tests from the production image
RUN rm -rf /app/packages/*/tests

# Set the PYTHONPATH to include the FreeCAD libraries
ENV PYTHONPATH="/usr/lib/freecad/lib"

# Set up a non-root user
RUN useradd --create-home --shell /bin/bash appuser
WORKDIR /home/appuser
USER appuser

# --- Test Stage ---
FROM production as test

# Switch to root to install test dependencies
USER root
WORKDIR /app

# Copy the tests into the test image. They were removed in the 'production' stage.
COPY ./packages/horn-geometry/tests /app/packages/horn-geometry/tests

# Install test dependencies. This will find the test extras in pyproject.toml
RUN uv pip install /app/packages/horn-geometry[test]

# Switch back to the appuser to run the tests
USER appuser
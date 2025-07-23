# --- Base Image ---
FROM python:3.10-slim as base

# --- Builder Stage ---
# This stage installs our package along with its dependencies
FROM base as builder
RUN pip install uv
WORKDIR /app
COPY ./packages/horn-analysis/pyproject.toml .
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-analysis /app/packages/horn-analysis
RUN uv pip install --system /app/packages/horn-analysis

# --- Production Stage ---
# This is the final, clean image for the service.
FROM base as production
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
WORKDIR /app

# Copy the installed packages from the builder stage.
COPY --from=builder /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

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

# The production stage already contains the main dependencies.
# We just need to copy the source again to install the '[test]' extras.
COPY ./packages/horn-analysis/pyproject.toml .
COPY ./packages/horn-analysis /app/packages/horn-analysis
RUN pip install /app/packages/horn-analysis[test]

# Copy the tests themselves
COPY ./packages/horn-analysis/tests /app/packages/horn-analysis/tests
RUN mkdir -p /home/appuser/.cache && chown -R appuser:appuser /home/appuser/.cache

USER appuser 
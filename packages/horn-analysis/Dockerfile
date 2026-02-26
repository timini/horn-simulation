# --- Base Image ---
FROM python:3.10-slim as base

# --- Builder Stage ---
# This stage installs our package along with its dependencies
FROM base as builder
RUN pip install uv
WORKDIR /app
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-drivers /app/packages/horn-drivers
COPY ./packages/horn-analysis /app/packages/horn-analysis
RUN uv pip install --system --no-deps /app/packages/horn-core \
    && uv pip install --system --no-deps /app/packages/horn-drivers \
    && uv pip install --system /app/packages/horn-analysis

# --- Production Stage ---
# This is the final, clean image for the service.
FROM base as production
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
WORKDIR /app

# Copy the installed packages from the builder stage.
COPY --from=builder /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

# --- Test Stage ---
# This stage adds test files and dependencies to the production image.
FROM production as test
WORKDIR /app
ENV PYTHONPATH /app

# The production stage already contains the main dependencies.
# We just need to copy the source again to install the '[test]' extras.
COPY ./packages/horn-core /app/packages/horn-core
COPY ./packages/horn-drivers /app/packages/horn-drivers
COPY ./packages/horn-analysis /app/packages/horn-analysis
RUN pip install --no-deps /app/packages/horn-core \
    /app/packages/horn-drivers \
    && pip install /app/packages/horn-analysis[test]

# Copy the tests themselves
COPY ./packages/horn-analysis/tests /app/packages/horn-analysis/tests
FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends octave && rm -rf /var/lib/apt/lists/*
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

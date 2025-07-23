# Dockerfile for native linux/arm64 builds
#
# This Dockerfile builds bempp-cl and its entire dependency stack from
# source. It produces a native linux/arm64 image that eliminates emulation,
# guaranteeing the best possible performance and stability on ARM64 hardware
# (e.g., Apple Silicon, AWS Graviton).
#
# The build process is time-consuming but creates a durable, reproducible,
# and high-performance environment for the project's lifetime.

# Stage 1: Build Environment
# Contains all development headers, compilers, and source code.
FROM ubuntu:22.04 AS builder

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    git \
    wget \
    unzip \
    # Python
    python3.10-dev \
    python3-pip \
    python3-venv \
    # Linear Algebra
    libblas-dev \
    liblapack-dev \
    # MPI for FEniCSx
    mpi-default-dev \
    libopenmpi-dev \
    # OpenCL for BEMPP
    opencl-headers \
    ocl-icd-opencl-dev \
    pocl-dev \
    # Other FEniCSx dependencies
    libpugixml-dev \
    # Cleanup
    && rm -rf /var/lib/apt/lists/*

# Create and activate a Python virtual environment
RUN python3.10 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN python -m pip install --upgrade pip

# Install core Python packages
RUN pip install numpy scipy numba "meshio>=4.0.16" "gmsh>=4.6.0" plotly pyopencl

# Build and install exafmm-t
RUN git clone https://github.com/exafmm/exafmm-t.git /tmp/exafmm-t
WORKDIR /tmp/exafmm-t
RUN ./configure && make -j$(nproc) && make install
RUN python setup.py install
WORKDIR /

# Build and install fenics-dolfinx
RUN git clone https://github.com/FEniCS/dolfinx.git /tmp/dolfinx
WORKDIR /tmp/dolfinx/cpp
RUN mkdir build && cd build && \
    cmake .. && \
    make -j$(nproc) && \
    make install
WORKDIR /tmp/dolfinx/python
RUN pip install -r build-requirements.txt
RUN pip install --no-build-isolation .
WORKDIR /

# Build and install bempp-cl
RUN git clone https://github.com/bempp/bempp-cl.git /tmp/bempp-cl
WORKDIR /tmp/bempp-cl
RUN pip install .
WORKDIR /

#---------------------------------------------------------------------

# Stage 2: Final Runtime Image
# Contains only the necessary runtime dependencies and the installed packages.
FROM ubuntu:22.04

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    python3.10 \
    python3.10-venv \
    libopenmpi3 \
    liblapack3 \
    libblas3 \
    ocl-icd-libopencl1 \
    pocl-opencl-icd \
    && rm -rf /var/lib/apt/lists/*

# Copy the virtual environment from the builder stage
COPY --from=builder /opt/venv /opt/venv

# Set the PATH to use the virtual environment's executables
ENV PATH="/opt/venv/bin:$PATH"

# Set a working directory for the application code
WORKDIR /app

# The container is now ready. The user can add their application code.
# For example:
# COPY ./src /app/src
# CMD ["python", "src/horn_bem_solver/solver.py"]

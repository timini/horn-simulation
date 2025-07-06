FROM ubuntu:22.04

# Set non-interactive frontend for package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-venv \
    python-is-python3 \
    freecad \
    && rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Set up the application directory
WORKDIR /app
COPY . .

# Set up FreeCAD Python path
ENV PYTHONPATH="/usr/lib/freecad/lib:${PYTHONPATH:-}"

# Final verification
RUN python -c "import pytest; print('Pytest found.')"
RUN python -c "import FreeCAD; print(f'Successfully imported FreeCAD version: {FreeCAD.Version()}')" 
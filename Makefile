# Makefile for the Horn Simulation Pipeline

.PHONY: all clean run-pipeline build-geometry-image build-solver-image test test-geometry test-solver

# --- Configuration ---
# Horn Parameters can be overridden from the command line
# e.g., make run-pipeline PROFILE=exponential
DRIVER_ID      ?= d220ti
PROFILE        ?= conical
THROAT_RADIUS  ?= 0.025
MOUTH_RADIUS   ?= 0.25
LENGTH         ?= 0.35
FREQ_MIN       ?= 100
FREQ_MAX       ?= 2000

# --- Docker Configuration ---
GEOMETRY_IMAGE := horn-geometry-app
GEOMETRY_TEST_IMAGE := $(GEOMETRY_IMAGE)-test
SOLVER_IMAGE   := horn-solver-app
SOLVER_TEST_IMAGE   := $(SOLVER_IMAGE)-test

# --- Directory and File Configuration ---
# Use a unique ID for each run to keep results separate.
RUN_ID         := $(shell uuidgen | cut -c -8)
TMP_DIR        := $(shell pwd)/tmp/$(RUN_ID)

STEP_FILE      := $(TMP_DIR)/horn_$(PROFILE)_$(LENGTH)m.stp
RESULTS_FILE   := $(TMP_DIR)/horn_$(PROFILE)_$(LENGTH)m.csv

# --- Python Environment ---
# Ensure we have a venv and the orchestrator is installed
VENV_DIR := .venv
ACTIVATE := source $(VENV_DIR)/bin/activate;
INSTALL_DEPS := $(ACTIVATE) pip install pandas

# Default target
all: run-pipeline

test: test-geometry test-solver

test-geometry: build-geometry-test-image
	@echo "\n--- Running Geometry Tests ---"
	@docker run --rm $(GEOMETRY_TEST_IMAGE) python -m pytest -s -v packages/horn-geometry/tests

test-solver: build-solver-test-image
	@echo "\n--- Running Solver Tests ---"
	@docker run --rm $(SOLVER_TEST_IMAGE) python -m pytest -s -v packages/horn-solver/tests

# Main pipeline target, depends on the final results file.
run-pipeline: $(RESULTS_FILE)
	@echo "\n--- Stage 5: Analysis ---"
	$(ACTIVATE) python scripts/run_analysis.py --results-file $(RESULTS_FILE)
	@echo "\n✅ Pipeline finished successfully!"
	@echo "   Results are in $(RESULTS_FILE)\n"

# --- File Generation Rules ---

# Rule to generate the final results file. Depends on the STEP file and the solver image.
$(RESULTS_FILE): $(STEP_FILE) build-solver-image $(VENV_DIR)/.setup_complete
	@echo "\n--- Stage 3: Solving (Meshing and Simulation) ---"
	@DRIVER_PARAMS_JSON=$$($(ACTIVATE) python scripts/get_driver_params.py --driver-id $(DRIVER_ID)); \
	docker run --rm \
		-v $(TMP_DIR):/data \
		$(SOLVER_IMAGE) \
		python -m horn_solver.solver_runner \
		--step-file /data/$(notdir $(STEP_FILE)) \
		--output-file /data/$(notdir $(RESULTS_FILE)) \
		--freq-min $(FREQ_MIN) \
		--freq-max $(FREQ_MAX) \
		--driver-params-json "$$DRIVER_PARAMS_JSON"

# Rule to generate the STEP file. Depends on the geometry Docker image.
$(STEP_FILE): build-geometry-image
	@echo "--- Creating temporary directory: $(TMP_DIR) ---"
	@mkdir -p $(TMP_DIR)
	@echo "\n--- Stage 2: Geometry Generation ---"
	@docker run --rm \
		-v $(TMP_DIR):/data \
		$(GEOMETRY_IMAGE) \
		python -m horn_geometry.generator \
		--profile $(PROFILE) \
		--throat $(THROAT_RADIUS) \
		--mouth $(MOUTH_RADIUS) \
		--length $(LENGTH) \
		--output-file /data/$(notdir $(STEP_FILE))

# --- Image Build Rules ---

# Using stamp files to avoid rebuilding images if the source hasn't changed.
build-geometry-image:
	@echo "\n--- Building Geometry Docker Image (production) ---"
	@docker build --target production -t $(GEOMETRY_IMAGE) -f packages/horn-geometry/Dockerfile .

build-geometry-test-image:
	@echo "\n--- Building Geometry Docker Image (test) ---"
	@docker build --target test -t $(GEOMETRY_TEST_IMAGE) -f packages/horn-geometry/Dockerfile .

build-solver-image:
	@echo "\n--- Building Solver Docker Image (production) ---"
	@docker build --target production -t $(SOLVER_IMAGE) -f packages/horn-solver/Dockerfile .

build-solver-test-image:
	@echo "\n--- Building Solver Docker Image (test) ---"
	@docker build --target test -t $(SOLVER_TEST_IMAGE) -f packages/horn-solver/Dockerfile .

# --- Python Environment Setup ---

# Target to set up the Python venv.
$(VENV_DIR)/.setup_complete:
	@echo "\n--- Setting up Python virtual environment ---"
	rm -rf $(VENV_DIR)
	python3 -m venv $(VENV_DIR)
	$(INSTALL_DEPS)
	@touch $@

# --- Housekeeping ---
clean:
	@echo "Cleaning up temporary run data..."
	@rm -rf tmp/*
	@echo "Done."

distclean: clean
	@echo "Cleaning up Python virtual environment..."
	@rm -rf $(VENV_DIR)
	@echo "Done." 
PACKAGES := horn-solver horn-geometry horn-analysis

.PHONY: all build test run clean

all: build

build:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):latest --target production -f ./packages/$(pkg)/Dockerfile .;)

test:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):test --target test -f ./packages/$(pkg)/Dockerfile .;)
	@for pkg in $(PACKAGES); do \
		echo "Running tests for $$pkg..."; \
		if [ "$$pkg" = "horn-solver" ]; then \
			docker run --rm -it "$$pkg:test" /bin/bash -c "echo '--- Verifying Environment ---'; echo 'PETSC_ARCH: $$PETSC_ARCH'; echo 'PYTHONPATH: $$PYTHONPATH'; python3 -c 'from petsc4py import PETSc; print(f\"PETSc.ScalarType: {PETSc.ScalarType}\")'; echo '--- Running Tests ---'; pytest /app/packages/horn-solver/tests"; \
		else \
			docker run --rm -it "$$pkg:test" pytest /app/tests; \
		fi \
	done

run:
	nextflow run main.nf -profile docker

clean:
	@echo "Cleaning up Docker images..."
	@docker rmi -f $(addsuffix :latest,$(PACKAGES)) $(addsuffix :test,$(PACKAGES)) || true
	@echo "Cleaning up Nextflow files..."
	@rm -rf .nextflow work .nextflow.log* trace.txt timeline.html report.html dag.dot dag.png 
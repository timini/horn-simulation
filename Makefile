PACKAGES := horn-solver horn-geometry horn-analysis

.PHONY: all build test run clean

all: help

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all		- Display this help message"
	@echo "  build		- Build all Docker images"
	@echo "  test		- Run all tests"
	@echo "  run		- Run the Nextflow pipeline"
	@echo "  clean		- Clean up Docker images and Nextflow files"

build:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):latest --target production -f ./packages/$(pkg)/Dockerfile .;)

test:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):test --target test -f ./packages/$(pkg)/Dockerfile .;)
	@for pkg in $(PACKAGES); do \
		echo "Running tests for $$pkg..."; \
		if [ "$$pkg" = "horn-solver" ]; then \
			docker run --rm "$pkg:test" pytest /app/packages/horn-solver/tests; \
		else \
			            else
			docker run --rm "$pkg:test" pytest /app/packages/horn-analysis/tests;
		fi \
		fi \
		done

run:
	nextflow run main.nf -profile docker

clean:
	@echo "Cleaning up Docker images..."
	@docker rmi -f $(addsuffix :latest,$(PACKAGES)) $(addsuffix :test,$(PACKAGES)) || true
	@echo "Cleaning up Nextflow files..."
	@rm -rf .nextflow work .nextflow.log* trace.txt timeline.html report.html dag.dot dag.png
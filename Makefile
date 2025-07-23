PACKAGES := horn-solver horn-geometry horn-analysis

.PHONY: all build test run clean

all: help

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all\t\t- Display this help message"
	@echo "  build\t\t- Build all Docker images"
	@echo "  test\t\t- Run all tests"
	@echo "  run\t\t- Run the Nextflow pipeline"
	@echo "  clean\t\t- Clean up Docker images and Nextflow files"

build:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):latest --target production -f ./packages/$(pkg)/Dockerfile .;)

test:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):test --target test -f ./packages/$(pkg)/Dockerfile .;)
	@for pkg in $(PACKAGES); do \
		echo "Running tests for $$pkg..."; \
		docker run --rm "$$pkg:test" pytest /app/packages/$$pkg/tests; \
	done

run:
	nextflow run main.nf -profile docker

clean:
	@echo "Cleaning up Docker images..."
	@docker rmi -f $(addsuffix :latest,$(PACKAGES)) $(addsuffix :test,$(PACKAGES)) || true
	@echo "Cleaning up Nextflow files..."
	@rm -rf .nextflow work .nextflow.log* trace.txt timeline.html report.html dag.dot dag.png
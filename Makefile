PACKAGES := horn-solver horn-geometry horn-analysis

.PHONY: all build test clean

all: build

build:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):latest --target production -f ./packages/$(pkg)/Dockerfile .;)

test:
	$(foreach pkg,$(PACKAGES),docker build -t $(pkg):test --target test -f ./packages/$(pkg)/Dockerfile .;)
	@for pkg in $(PACKAGES); do \
		echo "Running tests for $$pkg..."; \
		if [ "$$pkg" = "horn-solver" ]; then \
			docker run --rm -it "$$pkg:test" pytest /app/packages/horn-solver/tests; \
		else \
			docker run --rm -it "$$pkg:test" pytest /app/tests; \
		fi \
	done

clean:
	@docker rmi -f $(addsuffix :latest,$(PACKAGES)) $(addsuffix :test,$(PACKAGES)) || true 
.PHONY: help setup test lint atlas sweep benchmark all dag clean

help:
	@grep -E '^[a-z-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  %-12s %s\n", $$1, $$2}'

setup:      ## create the base conda environment and install the package
	conda env create -f envs/base.yaml || conda env update -f envs/base.yaml
	pip install -e ".[dev]"

test:       ## run the unit tests (no data required)
	pytest -q

lint:       ## ruff check
	ruff check src tests

atlas:      ## build the resolved atlas and stop before any TF inference
	snakemake --use-conda --cores all atlas

sweep:      ## grid search scored on the gating controls only; writes config/frozen.yaml
	snakemake --use-conda --cores all sweep

benchmark:  ## run through the positive-control benchmark and stop at the gate
	snakemake --use-conda --cores all benchmark_only

all:        ## full pipeline: gate must pass before the club prediction is produced
	snakemake --use-conda --cores all

dag:        ## render the workflow DAG
	snakemake --dag | dot -Tsvg > docs/dag.svg

clean:      ## remove results but keep downloaded raw data
	rm -rf results logs .snakemake

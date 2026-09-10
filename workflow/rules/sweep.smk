SWEEP_NETWORKS = config.get("sweep", {}).get("networks", ["aracne_pooled", "collectri"])


rule run_sweep:
    """Evaluate the parameter grid on the GATING contrasts only.

    Inputs deliberately exclude anything club-related: the sweep is scored on
    control recovery, and it cannot read the answer it is being tuned to protect.
    """
    input:
        h5ad="results/atlas/atlas_annotated.h5ad",
        networks=expand("results/networks/{network}.tsv", network=SWEEP_NETWORKS),
        markers="config/markers.yaml",
        benchmark="config/benchmark.yaml",
    output:
        table="results/sweep/sweep_results.csv",
        details="results/sweep/sweep_summary.json",
    params:
        families=SWEEP_NETWORKS,
    log:
        "logs/sweep/run_sweep.log",
    threads: 8
    conda:
        "../../envs/base.yaml"
    script:
        "../../src/basal_to_club/benchmark/sweep.py"

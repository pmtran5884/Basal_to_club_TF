rule concat_atlas:
    input:
        expand("data/interim/{dataset}_qc.h5ad", dataset=ATLAS_IDS),
    output:
        "data/interim/atlas_concat.h5ad",
    log:
        "logs/integrate/concat.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/integrate/concat.py"


rule integrate_atlas:
    """Batch-corrected latent space for neighbours/clustering/label transfer ONLY.

    Signatures are always recomputed from raw counts within donor; corrected
    expression never reaches VIPER.
    """
    input:
        "data/interim/atlas_concat.h5ad",
    output:
        h5ad="results/atlas/atlas_integrated.h5ad",
        metrics="results/atlas/integration_metrics.json",
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "logs/integrate/integrate.log",
    conda:
        "../../envs/scvi.yml"
    script:
        "../../src/basal_to_club/integrate/run_integration.py"

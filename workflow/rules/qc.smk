rule qc_dataset:
    """Per-sample adaptive QC, ambient correction, doublet removal, epithelial gate."""
    input:
        h5ad="data/interim/{dataset}_loaded.h5ad",
    output:
        h5ad="data/interim/{dataset}_qc.h5ad",
        metrics="results/qc/{dataset}_qc_metrics.json",
        fig="results/qc/{dataset}_qc.png",
    log:
        "logs/qc/{dataset}.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/qc/run_qc.py"


rule qc_summary:
    input:
        expand("results/qc/{dataset}_qc_metrics.json", dataset=ATLAS_IDS),
    output:
        "results/qc/qc_summary.csv",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/qc/summarize.py"

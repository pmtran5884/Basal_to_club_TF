rule transfer_labels:
    """scANVI/scArches surgery from the HLCA core reference."""
    input:
        query="results/atlas/atlas_integrated.h5ad",
        reference="data/interim/hlca_core_qc.h5ad",
    output:
        "results/atlas/atlas_transferred.h5ad",
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "logs/annotate/transfer.log",
    conda:
        "../../envs/scvi.yml"
    script:
        "../../src/basal_to_club/annotate/transfer.py"


rule annotate_atlas:
    """Arbitrate author label / transferred label / marker score; flag disagreements."""
    input:
        h5ad="results/atlas/atlas_transferred.h5ad",
        markers="config/markers.yaml",
    output:
        h5ad="results/atlas/atlas_annotated_coarse.h5ad",
        disagreement="results/annotate/label_disagreement.csv",
    log:
        "logs/annotate/annotate.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/annotate/run_annotate.py"

rule trajectory:
    input:
        h5ad="results/atlas/atlas_annotated.h5ad",
        validated="results/secretory/resolver_holdout_metrics.json",
    output:
        h5ad="results/trajectory/atlas_trajectory.h5ad",
        fates="results/trajectory/fate_probabilities.csv",
    log:
        "logs/trajectory.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/trajectory/run_trajectory.py"


rule signatures:
    """Three independent GES constructions per contrast. Agreement across them is
    the evidence; disagreement is reported rather than averaged away."""
    input:
        atlas="results/atlas/atlas_annotated.h5ad",
        traj="results/trajectory/atlas_trajectory.h5ad",
        meta="results/metacells/atlas_metacells.h5ad",
    output:
        "results/signatures/signatures.h5",
        "results/signatures/signature_summary.csv",
    log:
        "logs/signatures.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/signature/run_signatures.py"

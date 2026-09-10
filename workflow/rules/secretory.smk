rule resolve_secretory:
    """Re-partition the secretory compartment into club / goblet / hybrid.

    Removes submucosal gland cells first: gland mucous cells are MUC5B-high and
    gland serous cells LTF/LYZ-high, and both otherwise blur the club definition.
    """
    input:
        h5ad="results/atlas/atlas_annotated_coarse.h5ad",
        markers="config/markers.yaml",
    output:
        h5ad="results/atlas/atlas_annotated.h5ad",
        assignments="results/secretory/secretory_assignments.csv",
        fig="results/secretory/club_goblet_separation.png",
    log:
        "logs/secretory/resolve.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/secretory/run_resolver.py"


rule validate_resolver:
    """Two independent checks against GSE121600.

    1. Held-out in vivo bronchial biopsy (GSM3439925): balanced accuracy >= 0.85.
    2. ALI time course: club fraction must rise monotonically with ali_day, and the
       BEGM/PneumaCult contrast must move the goblet fraction. A resolver that gets
       the holdout labels right but the time course wrong is splitting a continuum.
    """
    input:
        h5ad="results/atlas/atlas_annotated.h5ad",
    output:
        metrics="results/secretory/resolver_holdout_metrics.json",
        fig="results/secretory/ali_timecourse_composition.png",
    log:
        "logs/secretory/validate.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/secretory/validate.py"


rule metacells:
    input:
        h5ad="results/atlas/atlas_annotated.h5ad",
    output:
        h5ad="results/metacells/atlas_metacells.h5ad",
    threads: 8
    log:
        "logs/metacells.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/trajectory/run_metacells.py"

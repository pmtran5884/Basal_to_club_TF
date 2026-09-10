rule benchmark:
    """Score the frozen positive controls. Gating and reported-only are scored
    separately; rare lineages never enter the gate or the sweep objective."""
    input:
        activity="results/viper/metaviper_activity.csv",
        by_family="results/viper/activity_by_family.csv",
        atlas="results/atlas/atlas_annotated.h5ad",
        bench="config/benchmark.yaml",
    output:
        scores="results/benchmark/control_tf_scores.csv",
        gate="results/benchmark/gate.json",
        curve="results/benchmark/lineage_recovery_curve.png",
    log:
        "logs/benchmark/score.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/benchmark/run_benchmark.py"


rule benchmark_report:
    input:
        scores="results/benchmark/control_tf_scores.csv",
        gate="results/benchmark/gate.json",
    output:
        "results/benchmark/benchmark_report.html",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/report/benchmark_report.py"


checkpoint gate:
    """HARD GATE. If the gating controls fail, the club prediction is not produced.
    A failed gate is a result about the configuration, not a reason to relax the config.
    """
    input:
        "results/benchmark/gate.json",
    output:
        "results/benchmark/.gate_passed",
    log:
        "logs/benchmark/gate.log",
    script:
        "../../src/basal_to_club/benchmark/check_gate.py"


rule club_prediction:
    input:
        gate="results/benchmark/.gate_passed",
        activity="results/viper/metaviper_activity.csv",
        by_family="results/viper/activity_by_family.csv",
    output:
        ranking="results/club/club_tf_ranking.csv",
        fig="results/club/club_tf_top.png",
    log:
        "logs/club/predict.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/benchmark/predict_club.py"


rule robustness:
    input:
        "results/club/club_tf_ranking.csv",
    output:
        "results/club/robustness.csv",
        "results/club/robustness.png",
    log:
        "logs/club/robustness.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/benchmark/robustness.py"

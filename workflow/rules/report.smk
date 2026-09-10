rule report:
    input:
        club="results/club/club_tf_ranking.csv",
        robustness="results/club/robustness.csv",
        gate="results/benchmark/gate.json",
        scores="results/benchmark/control_tf_scores.csv",
        qc="results/qc/qc_summary.csv",
        resolver="results/secretory/resolver_holdout_metrics.json",
        netqc="results/networks/network_qc.csv",
    output:
        "results/report/report.html",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/report/build_report.py"

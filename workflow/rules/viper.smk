rule viper:
    input:
        signatures="results/signatures/signatures.h5",
        network="results/networks/{network}.tsv",
        qc="results/networks/network_qc.csv",
    output:
        "results/viper/{network}_activity.csv",
    threads: 4
    log:
        "logs/viper/{network}.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/viper/run_viper.py"


rule metaviper:
    """Integrate across networks. De novo and prior families are also kept separate,
    because the benchmark requires a control TF to pass in BOTH families."""
    input:
        expand("results/viper/{network}_activity.csv", network=NETWORKS),
    output:
        combined="results/viper/metaviper_activity.csv",
        by_family="results/viper/activity_by_family.csv",
    log:
        "logs/viper/metaviper.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/viper/run_metaviper.py"

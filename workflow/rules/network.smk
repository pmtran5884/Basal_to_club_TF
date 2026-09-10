rule aracne:
    input:
        h5ad="results/metacells/atlas_metacells.h5ad",
    output:
        "results/networks/{network}.tsv",
    wildcard_constraints:
        network="aracne_.*",
    threads: 16
    resources:
        mem_mb=64000,
    log:
        "logs/network/{network}.log",
    conda:
        "../../envs/aracne.yml"
    script:
        "../../src/basal_to_club/network/run_aracne.py"


rule fetch_prior_libraries:
    """Download ChEA3 libraries once, with sha256 + date recorded. Versioned and
    offline-reproducible; the web API is deliberately not used at analysis time."""
    output:
        directory("data/external/chea3"),
        prov="results/provenance/chea3.json",
    log:
        "logs/network/fetch_chea3.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/network/fetch_prior.py"


rule build_prior_regulon:
    """Unsigned gene-set library -> signed, likelihood-weighted regulon."""
    input:
        libs="data/external/chea3",
        meta="results/metacells/atlas_metacells.h5ad",
    output:
        "results/networks/chea3_{library}.tsv",
    log:
        "logs/network/chea3_{library}.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/network/build_prior.py"


rule collectri:
    input:
        meta="results/metacells/atlas_metacells.h5ad",
    output:
        "results/networks/collectri.tsv",
    log:
        "logs/network/collectri.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/network/build_collectri.py"


rule network_qc:
    """Sanity gate on every network before it is allowed into VIPER."""
    input:
        expand("results/networks/{network}.tsv", network=NETWORKS),
    output:
        "results/networks/network_qc.csv",
    log:
        "logs/network/qc.log",
    conda:
        "../../envs/python.yml"
    script:
        "../../src/basal_to_club/network/qc.py"

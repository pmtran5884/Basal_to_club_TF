"""DREM (Dynamic Regulatory Events Miner) integration.

DREM is an Input-Output HMM over a time course: genes follow paths through a
bifurcating tree, and each bifurcation is annotated with the TFs whose targets
are over-represented on one branch. It is a different inference object from
VIPER -- VIPER scores regulon activity in a static two-group contrast, DREM
scores regulator association with a *temporal* split. The two therefore need
different inputs (a time series, not a contrast) and different benchmarks.

Implementation notes recovered from the DREM 2.0 Java source (jernst98/STEM_DREM):
  * batch invocation is `java -jar drem.jar -b settings.txt model.txt
    geneassign.txt TFSCOREDIR`
  * the TF-gene file's three-column form requires a 3-token header
    `TF<TAB>GENE<TAB>SCORE`; a 2-token header is parsed as the grid format
  * with TFSCOREDIR set, DREM writes one `path_*.txt` per out-edge of every
    split, with per-TF enrichment p-values ("Score Split" / "Score Overall")
"""

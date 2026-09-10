"""Score the HLCA prior-only run: gate, specificity, club ranking.

Substitution declared up front: the frozen benchmark requires a control to pass
in BOTH network families, where the families are de-novo (ARACNe) and prior.
This run has no ARACNe arm, so the two families here are the two independent
evidence types WITHIN ChEA3 -- expression-derived (ARCHS4, GTEx, Enrichr
queries) and binding-derived (ENCODE, Literature, ReMap ChIP-seq). The frozen
config is NOT edited; the relaxation is recorded in the output.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, "src")
from basal_to_club.benchmark.predict_club import rank_club_tfs
from basal_to_club.benchmark.score import check_specificity, score_gating
from basal_to_club.utils.io import load_yaml
from basal_to_club.viper.run_metaviper import integrate

RES = Path("results/hlca_prior_only")
EXPR = {"ARCHS4_Coexpression", "GTEx_Coexpression", "Enrichr_Queries"}
FAMILY = lambda n: "prior_expression" if n in EXPR else "prior_chipseq"  # noqa: E731

bench = load_yaml("config/benchmark.yaml")
markers = load_yaml("config/markers.yaml")
long = pd.read_csv(RES / "tf_activity_long.csv")
long["family"] = long["network"].map(FAMILY)

combined = integrate(long, ["tf", "contrast", "method"])
by_family = integrate(long, ["tf", "contrast", "method", "family"])
combined.to_csv(RES / "metaviper_combined.csv", index=False)
by_family.to_csv(RES / "metaviper_by_family.csv", index=False)

# Only contrasts actually computed in this scoped run are gated. The rest are
# reported as not evaluated -- never as passed.
COMPUTED = {"basal_to_goblet", "basal_to_ciliated"}
gating = {k: v for k, v in bench["gating"]["contrasts"].items() if k in COMPUTED}
not_evaluated = sorted(set(bench["gating"]["contrasts"]) - COMPUTED)

acts = {f: by_family[by_family["family"] == f].drop(columns="family")
        for f in sorted(by_family["family"].unique())}

results = {}
for label, req_patch in [("two_family", {}), ("single_family",
                          {"require_both_network_families": False})]:
    b = {"gating": {"contrasts": gating,
                    "requirements": {**bench["gating"]["requirements"], **req_patch}}}
    scores, gate = score_gating(acts, b)
    gate["contrasts_not_evaluated"] = not_evaluated
    gate["relaxation"] = req_patch or None
    results[label] = gate
    scores.to_csv(RES / f"gate_scores_{label}.csv", index=False)

spec = check_specificity(combined, bench["gating"]["requirements"]["specificity"])
spec.to_csv(RES / "specificity_combined.csv", index=False)

ranking = rank_club_tfs(combined, by_family.rename(columns={"family": "family"}),
                        markers.get("club_prior_plausible", []))
ranking.to_csv(RES / "club_tf_ranking.csv")

summary = {"gate": results, "specificity_combined": spec.to_dict("records"),
           "n_tf_ranked": int(len(ranking)),
           "families": {f: int(a["tf"].nunique()) for f, a in acts.items()}}
(RES / "gate_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(results, indent=2)[:2600])
print("\nspecificity:\n", spec.to_string(index=False))

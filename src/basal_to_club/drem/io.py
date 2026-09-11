"""Writers for the three DREM input files."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# Settings keys are those accepted by DREM_IO_Batch.parseDefaults.
SETTINGS_DEFAULTS: dict[str, str] = {
    "Spot_IDs_included_in_the_data_file": "false",
    "Normalize_Data": "Normalize data",          # subtract first timepoint
    "Change_should_be_based_on": "Maximum-Minimum",
    "Maximum_Number_of_Missing_Values": "0",
    "Minimum_Absolute_Log_Ratio_Expression": "1.0",
    "Minimum_Standard_Deviation": "0.0",
    "Filter_Gene_If_It_Has_No_Static_Input_Data": "false",
    "Maximum_number_of_paths_out_of_split": "3",
    "Model_selection_framework": "Penalized Likelihood",
    "Penalized_likelihood_node_penalty": "40",
    "Convergence_Likelihood_%": "0.01",
    "Minimum_score_improvement": "0.0",
    "Random_Seed": "1260",
    "Allow_Path_Merges": "false",
    "Regulator_Types_Used_For_Activity_Scoring": "None",
    "Use_transcription_factor-gene_interaction_data_to_build": "true",
}


def write_expression(series: pd.DataFrame, path: str | Path) -> Path:
    """genes x timepoints (log2 scale) -> DREM expression file.

    Column order is the time order; DREM subtracts the first column when
    Normalize_Data is "Normalize data", so absolute level does not matter.
    """
    path = Path(path)
    if series.isna().any().any():
        raise ValueError("DREM expression input must not contain NaN")
    with open(path, "w") as fh:
        fh.write("Gene\t" + "\t".join(str(c) for c in series.columns) + "\n")
        series.to_csv(fh, sep="\t", header=False, float_format="%.5f")
    return path


def write_tf_gene(edges: pd.DataFrame, path: str | Path,
                  tf_col: str = "tf", target_col: str = "target",
                  score_col: str | None = None) -> Path:
    """TF-target edges -> DREM three-column interaction file.

    The header must carry three tokens or the file is read as the grid format.
    """
    path = Path(path)
    out = edges[[tf_col, target_col]].copy()
    out["score"] = 1.0 if score_col is None else edges[score_col].astype(float).values
    out = out.drop_duplicates(subset=[tf_col, target_col])
    with open(path, "w") as fh:
        fh.write("TF\tGENE\tSCORE\n")
        out.to_csv(fh, sep="\t", header=False, index=False, float_format="%.4f")
    return path


def write_settings(path: str | Path, data_file: str | Path, tf_gene_file: str | Path,
                   **overrides: str) -> Path:
    path = Path(path)
    cfg = dict(SETTINGS_DEFAULTS)
    cfg.update({k: str(v) for k, v in overrides.items()})
    cfg["Data_File"] = str(Path(data_file).resolve())
    cfg["TF_gene_Interactions_File"] = str(Path(tf_gene_file).resolve())
    with open(path, "w") as fh:
        for k, v in cfg.items():
            fh.write(f"{k}\t{v}\n")
    return path

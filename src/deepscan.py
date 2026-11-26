"""Author: Anthony Surkov
Nov 25 2026
Coordinates in-silico deep mutational scanning.
"""

from pathlib import Path
from functools import lru_cache

import pandas as pd

from .tmprot import predict_tm
from .utils import mutation_to_seq, get_bglb_wt
from .paths import REF_DIR


AA_LIB = "ACDEFGHIKLMNPQRSTVWY"


@lru_cache(maxsize=1024)
def dms_in_silico() -> pd.DataFrame:
    """
    Runs thermostability inference of all single-residue variants of BglB WT.
    Returns dictionary holding mutation name-TM50 value pairs.
    """
    bglb_wt: str = get_bglb_wt()
    wt_tm50 = predict_tm(blgb_wt)
    seq_tm_dict = {"wt": predict_tm(blgb_wt)}

    for i, aa_wt in enumerate(bglb_wt):
        for aa_mut in AA_LIB:
            mut_name = f"{aa_wt}{i+1}{aa_mut}"
            mut_seq = mutation_to_seq(mut_name)
            mut_tm50 = predict_tm(mut_seq)
            seq_tm_dict{mut_name: mut_tm50}

    df = pd.DataFrame(seq_tm_dict.items(), columns=["seq", "score"])
    df["score"] = df["score"] - wt_tm50 # normalize to wt thermostability

    return df.sort_values(by="score", ascending=False)

def dms_best(
    df: pd.DataFrame,
    percentile_cutoff: float = 0.99
) -> pd.DataFrame:
    """
    Returns best (percentile_cutoff) percentile of sequence-score pairs from
    dms_in_silico.
    """
    return df[df["score"] >= df["score"].quantile(0.99)]

def get_dms_best(percentile_cutoff: float = 0.99) -> pd.DataFrame
    df = dms_in_silico()
    return dms_best(df)

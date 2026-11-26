from pathlib import Path
import os
import sys

import torch
import pandas as pd
from sklearn.model_selection import train_test_split
import joblib
from transformers import AutoTokenizer, EsmModel

from .tmprot_embeddings import (
    load_model,
    get_device,
    extract_embeddings
)
from .paths import REF_DIR, DATA_DIR


fasta_path = REF_DIR / "bglb.fasta"
csv_path = DATA_DIR / "bglb.csv"
data_dir = DATA_DIR / "tmprot"
model_path = DATA_DIR / "tmprot/model_output/mlp.joblib"


def read_wt(fp: str) -> str:
    with open(fp, "r") as f:
        line = f.readline()
        return f.readline().strip()

def translate_mutations(fp: str, wt: str) -> dict[str, float]:
    df = pd.read_csv(csv_path)
    variants_found = []
    tm_mapping = {}
    for index, row in df.iterrows():
        if pd.isna(row["Tm (C)"]):
            continue
        if row["Variant"] in variants_found:
            continue
        variants_found.append(row["Variant"])
        seq = ""
        if row["Variant"] == "WT":
            seq = wt
        else:
            var = row["Variant"]
            old_aa = var[0]
            new_aa = var[len(var) - 1]
            index = int(var[1: len(var) - 1]) + 3
            if old_aa != wt[index - 1]:
                print(var)
                print(wt[index - 1])
                continue
            seq = wt[:index - 1] + new_aa + wt[index:]
            if len(seq) != len(wt):
                print("wrong sizing")
                continue
        tm_mapping[seq] = float(row["Tm (C)"])
    return tm_mapping

def convert_format(mapping: dict[str, float]) -> pd.DataFrame:
    df = pd.DataFrame(
    [
        {
            "ProteinID": f"BglB_{i}",
            "sequence": sequence,
            "Tm": tm,
            "Length": len(sequence),
        }
        for i, (sequence, tm) in enumerate(mapping.items())
    ])
    return df

def split_data_to_csv(df: pd.DataFrame) -> None:
    train_df, temp_df = train_test_split(df, test_size=0.2, random_state=42,)
    val_df, test_df = train_test_split(temp_df, test_size=0.5, random_state=42,)
    os.makedirs(data_dir, exist_ok=True)
    train_df.to_csv(f"{data_dir}/train.csv", index=False)
    val_df.to_csv(f"{data_dir}/val.csv", index=False)
    test_df.to_csv(f"{data_dir}/test.csv", index=False)

def predict_tm(seq: str) -> float:
    '''
    Predicts Tm given sequence. Import into other files
    '''
    mlp = joblib.load(model_path)
    device = get_device()
    tokenizer, esm_model = load_model(device)
    embedding = extract_embeddings([seq], tokenizer, esm_model, device)[0]
    X = pd.DataFrame(
        [embedding],
        columns=[f"embedding_{i}" for i in range(len(embedding))]
    )
    return float(mlp.predict(X)[0])

if __name__ == "__main__":
    """
    Run: python fine-tuning.py <fasta_path> <csv_path> <data_directory>
    Note: the model path will be defined by this file
    """
    # comment out for testing/if not needed
    fasta_path = sys.argv[1]
    csv_path = sys.argv[2]
    data_dir = sys.argv[3]

    seq = read_wt(fasta_path)
    tms = translate_mutations(csv_path, seq)
    df = convert_format(tms)
    split_data_to_csv(df)

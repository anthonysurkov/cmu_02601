"""
Author: Anthony Surkov
Date: Oct 19 2026
Featurizer functions (part 6)

Inputs: List of sequences of BglB mutations
Outputs: Dictionary of sequences and ProtBERT descriptor tensors
"""

from functools import cache
import re

import torch
import pandas as pd
from transformers import BertTokenizer, BertModel
from tqdm import tqdm


@cache
def init_protbert() -> tuple[BertTokenizer, BertModel]:
    """
    Initializes tokenizer and model from Rostlab's ProtBERT implementation.
    Caches both to avoid reinitialization on sequential uses.
    The tokenizer converts amino acid letter codes into ProtBERT model inputs
    (tokens). The model generates, among other things, context-aware embeddings
    for the tokenized sequence that are then used as downstream features.
    """
    tokenizer = BertTokenizer.from_pretrained(
        "Rostlab/prot_bert",
        do_lower_case=False
    )
    model = BertModel.from_pretrained("Rostlab/prot_bert")
    model.eval()

    return tokenizer, model

@torch.inference_mode()
def get_protbert_embedding(seq: str) -> torch.Tensor:
    """
    Retrieves ProtBERT embeddings for a given sequence. Coerces embeddings to
    [L, 1024] shape where L is the protein length and 1024 is the hard-coded
    Rostlab embedding size.
    """
    tokenizer, model = init_protbert()
    tokenized_seq = tokenizer( # [1, L+2], with (CLS, L tokens, SEP)
        seq,
        return_tensors="pt",
        padding=True
    )
    output = model(**tokenized_seq) # [1, L+2, 1024]; 1024 hardcoded for rostlab
    residue_embeddings = output.last_hidden_state

    # https://huggingface.co/docs/transformers/glossary
    residue_embeddings = residue_embeddings[:, 1:-1, :] # [1, L, 1024]
    return residue_embeddings.squeeze(0) # [L, 1024]

def get_protbert_dict(seqs: list[str]) -> dict[str, torch.Tensor]:
    """
    Retrieves ProtBERT embeddings for a list of sequences.
    """
    embeddings = []
    for seq in tqdm(seqs, desc="Getting protbert dict", total=len(seqs)):
        spaced_seq = " ".join([s for s in seq]) # expected by ProtBERT
        embedding = get_protbert_embedding(spaced_seq)
        embeddings.append(embedding)

    if not embeddings:
        raise ValueError("No embeddings retrieved")

    return dict(zip(seqs, embeddings, strict=True))

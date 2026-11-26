import pandas as pd
from pathlib import Path

from .paths import DATA_DIR


"""
An alignment-based algorithm will find conserved regions across the ancestral sequences of the target protein, which are estimated from extant homologs by ASR.
The homologs will be identified using the BLAST algorithm2 with the Biopython API wrapper.3 The amino-acid substitution model used by ASR will be selected using IQ-TREE’s ModelFinder4.
Our in-house algorithm will select the best candidate from (1) that does not propose changing the conserved regions.
The mutated protein will be evaluated with a machine learning model trained on assay activities of single variant mutants of the protein, with features computed by ProtBERT.5
"""

def main():
    candidates = get_dms_best() # get best theoretical candidates for
                                # thermostability increases


if __name__ == "__main__":
    main()

from .blast import get_homologs
from .deepscan import get_dms_best
from .final_decision import select_viable_mutations
from .iqtree import find_model
from .msa import msa_scored
from .paths import DATA_DIR
from .training_eval import evaluate_kfold, infer_tm50
from .truth import BglBRepository
from .utils import get_bglb_wt


"""
Scientific purpose:
An alignment-based algorithm will find conserved regions across the ancestral
sequences of the target protein, which are estimated from extant homologs by
ASR.
The homologs will be identified using the BLAST algorithm with the Biopython
API wrapper. The amino-acid substitution model used by ASR will be selected
using IQ-TREE’s ModelFinder.
Our in-house algorithm will select the best candidate(s) from conserved regions
that does not propose changing the conserved regions.
The mutated protein will be evaluated with a machine learning model trained on
assay activities of single variant mutants of the protein, with features
computed by ProtBERT.
"""

def main():
    all_mutants = get_dms_best() # get best theoretical candidates for
                                # thermostability increases
    wt_seq = get_bglb_wt()
    bglb_homologs = get_homologs(wt_seq) # identify all homologs of BglB thru
                                         # BLAST API call

    asr_model = find_model(bglb_homologs) # identify best-fitting empirical ASR
                                # model for BglB homologs, via IQ-TREE call

    aligned_homologs, _ = msa_scored(bglb_homologs) # align BglB
                                # homologs, get their per-base agreement scores

    candidates = select_viable_mutations( # select most-viable mutations based
        all_mutants,            # on the degree to which they affect conserved
        conserved_positions     # bases
    )

    bglb_repo = BglBRepository(DATA_DIR / "bglb.csv") # load BglB dataset
    xgboost_models, mean_val_rmse, test_rmse = evaluate_kfold( # train XGBoost
        bglb_repo.yield_train_val_folds() # models over known single-variant
    )                           # mutant data

    candidate_tm50s = infer_tm50(xgboost_models, candidates) # infer tm50's of
    print(candidate_tm50s)      # candidates using models

if __name__ == "__main__":
    main()

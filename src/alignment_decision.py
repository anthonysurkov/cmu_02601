"""
(4a) Alignment decision algorithm

Takes: 
    - list[str] of aligned sequences
    - list[float] of agreement between residues across aligned sequences
Outputs: 
    -dict[int, str] of all positions strongly conserved in ASR (residue agreement
    within some parameter theta)
"""


from pathlib import Path


# Default theta is variable (we can change later)
DEFAULT_THETA = 0.7
GAP = "-"


def _get_amino_residue(protein_strings: list[str], 
                       string_index: int) -> str | None:
    """
    Searches for the most frequent amino acid at a specific alignment location.

    This helper function takes a number of sequences at a particular length, and 
    builds a frequency map (dictionary) and returns the residue with the highest 
    count
    """
    temp_dict: dict[str, int] = {}
    max_count = 0

    # Counting residue frequencies at position string_index across all proteins 
    for protein in protein_strings:
        residue = protein[string_index]
        temp_dict[residue] = temp_dict.get(residue, 0) + 1

    # Finding residue with highest count 
    for key, value in temp_dict.items():
        if value > max_count:
            max_count = value
            most_common_residue = key 

    if most_common_residue == GAP:
        return None
    else:
        return most_common_residue 
    

def find_conserved_positions(protein_strings: list[str], 
                            agreements: list[float], 
                            theta: float = DEFAULT_THETA,) -> dict[int, str]:
    """
    Searching for alignment positions that are strictly above theta.

    This function find_conserved_positions iterates through agreement scores
    (from part 4) and compares with the residues the proteins.
    """
    # Initializing final common residues dictionary 
    common_residues: dict[int, str] = {}

    # Checking agreement scores against theta
    for index, value in enumerate(agreements):
        if value >= theta:
           winner_residue = _get_amino_residue(protein_strings, index)
           if winner_residue is not None:
               common_residues[index] = winner_residue
    return common_residues


# Testing
if __name__ == "__main__":
    # Took first 8 aminos of Ecoli Bglb (first)
    protein_strings = ["MKTKFPKD","MKT-FPKD","MK--FPKD","MKA-FPRD"]
    agreements = [1.0, 0.9, 0.65, 0.4, 0.85, 0.5, 0.95, 0.8]

    conserved_pos = find_conserved_positions(protein_strings, agreements)
    print(conserved_pos)













    
        




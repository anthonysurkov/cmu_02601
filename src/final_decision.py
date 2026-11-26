"""
(5) Final decision algorithm

Takes: 
    - dict[int, str] of thermostable predicted mutations from ML/DMS (Task 1b)
    - dict[int, str] of strongly conserved positions from ASR (Task 4a)
Outputs: 
    - dict[int, str] of filtered, viable mutations for experimental screening
"""


from pathlib import Path


def select_viable_mutations(
    thermo_mutations: dict[int, str], 
    conserved_positions: dict[int, str],
) -> dict[int, str]:
    """
    Filters machine learning thermostable predictions against the ASR consensus.

    This function cross-references predicted mutations with strongly conserved 
    ancestral positions to satisfy viable variant selection.
    """
    # Initializing final viable mutations dictionary 
    viable_mutations: dict[int, str] = {}

    # Checking predicted mutations against conserved positions
    for position, mutation in thermo_mutations.items():
        if position not in conserved_positions:
            viable_mutations[position] = mutation

    return viable_mutations


# Testing
if __name__ == "__main__":
    test_thermo_mutations = {1: "K", 2: "R", 15: "V"}
    test_conserved_positions = {1: "M", 5: "F"}

    viable_pos = select_viable_mutations(test_thermo_mutations, test_conserved_positions)
    print(viable_pos)
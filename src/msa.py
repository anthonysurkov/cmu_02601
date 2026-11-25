"""Author: Anthony Surkov
Date: Nov 25 2026
Smith-waterman aligner functions.

Inputs: list[str] of DNA
Outputs: list[str] of aligned DNA, list[float] of base agreements
"""

import argparse


def find_aligned_substr(dna1: str, dna2: str) -> list[str]:
    """Return aligned pieces using the Smith-Waterman scoring and traceback
    rules.
    Matches contribute score of 3, gaps contribute score of 2, and mismatches
    contribute zero.
    """
    if not dna1 or not dna2:
        return ["", ""]

    n, m = len(dna1), len(dna2)
    matrix = [[0] * (m + 1) for _ in range(n + 1)]
    global_max = 0
    pos_i, pos_j = 0, 0

    for i in range(n):
        for j in range(m):
            align = matrix[i][j] + 3 if dna1[i] == dna2[j] else 0
            top_gap = matrix[i + 1][j] - 2
            left_gap = matrix[i][j + 1] - 2
            score = max(0, align, top_gap, left_gap)
            matrix[i + 1][j + 1] = score
            if score > global_max:
                global_max = score
                pos_i, pos_j = i, j

    if global_max == 0:
        return ["", ""]

    local_dna1 = [dna1[pos_i]]
    local_dna2 = [dna2[pos_j]]

    while matrix[pos_i][pos_j] != 0:
        best = 0
        diag = matrix[pos_i - 1][pos_j - 1]
        up = matrix[pos_i][pos_j - 1]
        left = matrix[pos_i - 1][pos_j]

        if diag >= best:
            best = diag
            pos_i -= 1
            pos_j -= 1
        if up > best:
            best = up
            pos_j -= 1
        if left > best:
            pos_i -= 1

        if pos_i < 0 or pos_j < 0:
            break
        local_dna1.append(dna1[pos_i])
        local_dna2.append(dna2[pos_j])

    return ["".join(reversed(local_dna1)), "".join(reversed(local_dna2))]

def find_substr(a: str, b: str) -> int:
    """Return the first index of b in a, or -1 if it is absent."""
    return a.find(b)

def calculate_base_agreement(dna: list[str]) -> list[float]:
    """Return the most common character's fraction at each aligned position."""
    if not dna:
        return []
    if any(len(sequence) != len(dna[0]) for sequence in dna):
        raise ValueError("Aligned sequences must have equal lengths")

    agreements = []
    for column in zip(*dna):
        counts = {}
        for base in column:
            if base in counts:
                if base != '-':
                    counts[base] += 1
            else:
                counts[base] = 1
        agreements.append(max(counts.values()) / len(dna))
    return agreements

def msa(dna: list[str]) -> list[str]:
    """
    Multiple sequence alignment via incremental Smith-Waterman.
    Each sequence is positioned relative to the preceding input sequence using
    their local alignment.
    """
    if len(dna) < 2:
        raise ValueError("Cannot align less than 2 sequences")

    offsets = [0]
    for i in range(len(dna) - 1):
        previous, current = find_aligned_substr(dna[i], dna[i + 1])
        if not previous or not current:
            raise ValueError(
                f"Cannot align sequences {i} and {i + 1}"
            )
        shift = find_substr(dna[i], previous) - find_substr(dna[i + 1], current)
        offsets.append(offsets[-1] + shift)

    left = min(offsets)
    right = max(
        offset + len(sequence)
        for offset, sequence in zip(offsets, dna)
    )
    return [
        "-" * (offset - left) +
        sequence + "-" * (right - offset - len(sequence))
        for offset, sequence in zip(offsets, dna)
    ]

def msa_scored(dna: list[str]) -> tuple[list[str], list[float]]:
    """
    Intended handle for MSA use.
    Returns aligned sequences and their associated agreement scores.
    """
    aligned = msa(dna)
    scores = calculate_base_agreement(aligned)
    return (aligned, scores)


# debug land
def main(dna1: str, dna2: str) -> None:
    aligned, scores = sw_align_and_score([dna1, dna2])
    print(aligned)
    print(scores)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dna1", type=str, required=True)
    parser.add_argument("--dna2", type=str, required=True)
    args = parser.parse_args()
    main(args.dna1, args.dna2)

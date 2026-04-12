"""
Mutation Detection Module
=========================
Functions for detecting, classifying, and summarizing point mutations
between two DNA sequences. Supports transition/transversion classification
and mutation rate calculation.
"""


def classify_mutation(base1: str, base2: str) -> str:
    """Classify a point mutation as a transition or transversion.

    - Transition: Purine↔Purine (A↔G) or Pyrimidine↔Pyrimidine (C↔T)
    - Transversion: Purine↔Pyrimidine (e.g., A↔C, G↔T)

    Args:
        base1: Original nucleotide base.
        base2: Mutated nucleotide base.

    Returns:
        'transition' or 'transversion'.
    """
    purines = {"A", "G"}
    pyrimidines = {"C", "T"}

    if (base1 in purines and base2 in purines) or \
       (base1 in pyrimidines and base2 in pyrimidines):
        return "transition"
    else:
        return "transversion"


def detect_point_mutations(seq1: str, seq2: str) -> list[dict]:
    """Detect all point mutations between two aligned DNA sequences.

    Compares sequences position by position up to the length of the
    shorter sequence.

    Args:
        seq1: Reference DNA sequence.
        seq2: Query DNA sequence.

    Returns:
        List of mutation dicts with keys: position, original, mutated, type.
    """
    mutations = []
    min_len = min(len(seq1), len(seq2))

    for i in range(min_len):
        if seq1[i] != seq2[i]:
            mutation_type = classify_mutation(seq1[i], seq2[i])
            mutations.append({
                "position": i + 1,
                "original": seq1[i],
                "mutated": seq2[i],
                "type": mutation_type
            })

    return mutations


def calculate_mutation_rate(seq1: str, seq2: str, mutations: list) -> float:
    """Calculate the mutation rate between two sequences.

    Mutation rate = (number of mutations / comparable length) × 100

    Args:
        seq1: Reference DNA sequence.
        seq2: Query DNA sequence.
        mutations: List of detected mutations.

    Returns:
        Mutation rate as a percentage rounded to 2 decimal places.
    """
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    return round(len(mutations) / min_len * 100, 2)


# Keep backward compatibility with the old misspelled name
calculate_muatation_rate = calculate_mutation_rate


def summarize_mutations(mutations: list[dict]) -> dict:
    """Summarize mutations by type (transition vs transversion).

    Args:
        mutations: List of mutation dicts from detect_point_mutations().

    Returns:
        Dict mapping mutation type to count.
    """
    summary = {}
    for m in mutations:
        mtype = m["type"]
        summary[mtype] = summary.get(mtype, 0) + 1
    return summary

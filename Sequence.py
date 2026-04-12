"""
Sequence Analysis Module
========================
Core functions for parsing FASTA files and performing DNA sequence analysis
including GC content calculation, base composition, ORF finding, complement
generation, and protein translation.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def parse_fasta(file_path: str) -> list[dict]:
    """Parse a FASTA file and return a list of sequence records.

    Args:
        file_path: Path to the FASTA format file.

    Returns:
        List of dicts, each containing 'id' and 'sequence' keys.
    """
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append({
            "id": record.id,
            "sequence": str(record.seq).upper()
        })
    return records


def calculate_gc_content(sequence: str) -> float:
    """Calculate the GC content percentage of a DNA sequence.

    GC content is the proportion of Guanine (G) and Cytosine (C) bases
    in the sequence, expressed as a percentage.

    Args:
        sequence: DNA sequence string (uppercase).

    Returns:
        GC content as a percentage rounded to 2 decimal places.
    """
    g = sequence.count("G")
    c = sequence.count("C")
    total = len(sequence)
    if total == 0:
        return 0.0
    return round((g + c) / total * 100, 2)


def calculate_base_composition(sequence: str) -> dict:
    """Calculate the count and percentage of each nucleotide base.

    Args:
        sequence: DNA sequence string (uppercase).

    Returns:
        Dict mapping each base (A, T, G, C) to its count and percentage.
    """
    total = len(sequence)
    bases = ["A", "T", "G", "C"]
    composition = {}
    for base in bases:
        count = sequence.count(base)
        composition[base] = {
            "count": count,
            "percentage": round(count / total * 100, 2)
        }
    return composition


def find_orfs(sequence: str, min_length: int = 30) -> list[dict]:
    """Find all Open Reading Frames (ORFs) in the three forward reading frames.

    An ORF is defined as a stretch of DNA beginning with a start codon (ATG)
    and ending with a stop codon (TAA, TAG, or TGA).

    Args:
        sequence: DNA sequence string (uppercase).
        min_length: Minimum ORF length in base pairs (default: 30).

    Returns:
        List of dicts with keys: frame, start, end, length, sequence.
    """
    stop_codons = {"TAA", "TGA", "TAG"}
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == "ATG":
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop = sequence[j:j+3]
                    if stop in stop_codons:
                        orf_seq = sequence[i:j+3]
                        if len(orf_seq) >= min_length:
                            orfs.append({
                                "frame": frame + 1,
                                "start": i,
                                "end": j + 3,
                                "length": len(orf_seq),
                                "sequence": orf_seq
                            })
                        i = j + 3
                        break
                else:
                    i += 3
            else:
                i += 3
    return orfs


def get_complements(sequence: str) -> str:
    """Return the complementary DNA strand.

    Args:
        sequence: DNA sequence string.

    Returns:
        Complement sequence (A↔T, G↔C).
    """
    seq_obj = Seq(sequence)
    return str(seq_obj.complement())


def get_reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    The reverse complement is the complementary strand read in
    the 3' to 5' direction.

    Args:
        sequence: DNA sequence string.

    Returns:
        Reverse complement sequence.
    """
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())


def translate_sequence(sequence: str) -> str:
    """Translate a DNA sequence to a protein sequence.

    Translation starts from the first ATG (start codon) found in
    the sequence and continues until a stop codon is reached.

    Args:
        sequence: DNA sequence string.

    Returns:
        Protein amino acid sequence string, or an error message.
    """
    try:
        seq_obj = Seq(sequence)
        start = sequence.find("ATG")
        if start == -1:
            return "no start codon found"
        trimmed = seq_obj[start:]
        trimmed = trimmed[:len(trimmed) - (len(trimmed) % 3)]
        protein = trimmed.translate(to_stop=True)
        return str(protein)
    except Exception as e:
        return f"translation error: {e}"


def analyze_sequence(record: dict) -> dict:
    """Perform comprehensive analysis on a single DNA sequence record.

    Runs all analysis functions (GC content, base composition, ORF finding,
    complement, reverse complement, and protein translation) on the sequence.

    Args:
        record: Dict with 'id' and 'sequence' keys.

    Returns:
        Dict containing all analysis results for the sequence.
    """
    seq = record["sequence"]
    return {
        "id": record["id"],
        "sequence": seq,
        "length": len(seq),
        "gc_content": calculate_gc_content(seq),
        "base_composition": calculate_base_composition(seq),
        "orfs": find_orfs(seq),
        "complement": get_complements(seq),
        "reverse_complement": get_reverse_complement(seq),
        "protein": translate_sequence(seq)
    }
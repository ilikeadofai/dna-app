"""Shared DNA constraint constants for educational sequence analysis."""

DNA_BASES = ("A", "C", "G", "T")
DNA_BASE_SET = frozenset(DNA_BASES)

COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}

FORBIDDEN_MOTIFS = (
    "GAATTC",  # EcoRI
    "AAGCTT",  # HindIII
    "GGATCC",  # BamHI
)

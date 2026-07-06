"""Fixed 2-bit DNA codec.

This module provides a simple local, reversible mapping between bytes and DNA
bases. It is intentionally deterministic and does not use external APIs or
lookup tables.
"""

BIT_TO_BASE = {
    "00": "A",
    "01": "C",
    "10": "G",
    "11": "T",
}

BASE_TO_BIT = {base: bits for bits, base in BIT_TO_BASE.items()}
VALID_BASES = frozenset(BASE_TO_BIT)


def bytes_to_fixed_dna(data: bytes) -> str:
    """Encode bytes into DNA using the fixed 2-bit mapping.

    One byte is encoded as four DNA bases using most-significant-bit-first
    ordering. For example, ``0b00011011`` becomes ``ACGT``.
    """
    dna_parts: list[str] = []

    for byte in data:
        bits = f"{byte:08b}"
        dna_parts.extend(
            BIT_TO_BASE[bits[index : index + 2]] for index in range(0, 8, 2)
        )

    return "".join(dna_parts)


def fixed_dna_to_bytes(seq: str) -> bytes:
    """Decode fixed 2-bit DNA back into bytes.

    Raises:
        ValueError: If the sequence contains characters outside ``A/C/G/T`` or
            if its length is not a multiple of four bases.
    """
    if len(seq) % 4 != 0:
        raise ValueError("Fixed 2-bit DNA length must be a multiple of 4 bases")

    invalid = [
        (index, base) for index, base in enumerate(seq) if base not in VALID_BASES
    ]
    if invalid:
        index, base = invalid[0]
        raise ValueError(
            f"Invalid DNA base {base!r} at position {index}; expected only A, C, G, or T"
        )

    decoded = bytearray()
    for index in range(0, len(seq), 4):
        bits = "".join(BASE_TO_BIT[base] for base in seq[index : index + 4])
        decoded.append(int(bits, 2))

    return bytes(decoded)


def text_to_fixed_dna(text: str) -> str:
    """Encode Unicode text as UTF-8 bytes, then fixed 2-bit DNA."""
    return bytes_to_fixed_dna(text.encode("utf-8"))


def fixed_dna_to_text(seq: str) -> str:
    """Decode fixed 2-bit DNA into UTF-8 text."""
    return fixed_dna_to_bytes(seq).decode("utf-8")

"""Rotating ternary DNA codec.

The rotating codec converts bytes to base-3 trits, then maps each trit to one
of the three DNA bases that is not equal to the previous base. This prevents
homopolymers by construction while remaining fully reversible when the original
byte length, start base, and mapping id are known.
"""

from __future__ import annotations

from itertools import permutations

from .constraints import DNA_BASE_SET, DNA_BASES

MAPPING_IDS = tuple(range(6))
START_BASES = DNA_BASES


def bytes_to_trits(data: bytes) -> list[int]:
    """Convert bytes to a fixed-length list of base-3 digits.

    The trit length is the smallest number of trits that can represent every
    possible byte string of the same length. Leading zero bytes are therefore
    preserved when decoding with the original byte length.
    """
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")

    trit_count = _trit_count_for_byte_length(len(data))
    if trit_count == 0:
        return []

    value = int.from_bytes(data, byteorder="big")
    trits = [0] * trit_count

    for index in range(trit_count - 1, -1, -1):
        value, trit = divmod(value, 3)
        trits[index] = trit

    return trits


def trits_to_bytes(trits: list[int], length: int) -> bytes:
    """Convert a fixed-length trit list back to bytes.

    Args:
        trits: Base-3 digits produced by ``bytes_to_trits`` or decoded from a
            rotating DNA sequence.
        length: Original byte length to restore leading zero bytes.

    Raises:
        ValueError: If trits are invalid, the trit count does not match the
            requested byte length, or the represented integer does not fit.
    """
    _validate_byte_length(length)
    _validate_trits(trits)

    expected_trit_count = _trit_count_for_byte_length(length)
    if len(trits) != expected_trit_count:
        raise ValueError(
            f"Expected {expected_trit_count} trits for {length} bytes, got {len(trits)}"
        )

    value = 0
    for trit in trits:
        value = (value * 3) + trit

    max_value = 1 << (8 * length)
    if value >= max_value:
        raise ValueError(f"Trit value does not fit in {length} bytes")

    return value.to_bytes(length, byteorder="big")


def get_candidate_bases(previous_base: str, mapping_id: int = 0) -> list[str]:
    """Return the three allowed next bases for a previous base and mapping id."""
    previous_base = _validate_base(previous_base, name="previous_base")
    _validate_mapping_id(mapping_id)

    candidates = [base for base in DNA_BASES if base != previous_base]
    return list(tuple(permutations(candidates))[mapping_id])


def rotating_encode(data: bytes, start_base: str = "A", mapping_id: int = 0) -> str:
    """Encode bytes into homopolymer-free rotating ternary DNA."""
    start_base = _validate_base(start_base, name="start_base")
    _validate_mapping_id(mapping_id)

    previous_base = start_base
    dna_parts: list[str] = []

    for trit in bytes_to_trits(data):
        candidates = get_candidate_bases(previous_base, mapping_id=mapping_id)
        base = candidates[trit]
        dna_parts.append(base)
        previous_base = base

    return "".join(dna_parts)


def rotating_decode(
    seq: str, length: int, start_base: str = "A", mapping_id: int = 0
) -> bytes:
    """Decode rotating ternary DNA back to bytes."""
    _validate_byte_length(length)
    start_base = _validate_base(start_base, name="start_base")
    _validate_mapping_id(mapping_id)
    _validate_dna_sequence_strict(seq)

    expected_trit_count = _trit_count_for_byte_length(length)
    if len(seq) != expected_trit_count:
        raise ValueError(
            f"Expected rotating DNA length {expected_trit_count} for {length} bytes, got {len(seq)}"
        )

    previous_base = start_base
    trits: list[int] = []

    for index, base in enumerate(seq):
        candidates = get_candidate_bases(previous_base, mapping_id=mapping_id)
        if base not in candidates:
            raise ValueError(
                f"Invalid rotating DNA at position {index}: base {base!r} cannot follow {previous_base!r}"
            )

        trits.append(candidates.index(base))
        previous_base = base

    return trits_to_bytes(trits, length)


def _trit_count_for_byte_length(length: int) -> int:
    _validate_byte_length(length)
    if length == 0:
        return 0

    required_values = 1 << (8 * length)
    capacity = 1
    trit_count = 0

    while capacity < required_values:
        capacity *= 3
        trit_count += 1

    return trit_count


def _validate_byte_length(length: int) -> None:
    if type(length) is not int:
        raise TypeError("length must be an integer")
    if length < 0:
        raise ValueError("length must be non-negative")


def _validate_mapping_id(mapping_id: int) -> None:
    if type(mapping_id) is not int:
        raise TypeError("mapping_id must be an integer")
    if mapping_id not in MAPPING_IDS:
        raise ValueError("mapping_id must be between 0 and 5")


def _validate_base(base: str, name: str) -> str:
    if not isinstance(base, str):
        raise TypeError(f"{name} must be a string")
    if len(base) != 1 or base not in DNA_BASE_SET:
        raise ValueError(f"{name} must be one of A, C, G, or T")
    return base


def _validate_dna_sequence_strict(seq: str) -> None:
    if not isinstance(seq, str):
        raise TypeError("seq must be a string")

    for index, base in enumerate(seq):
        if base not in DNA_BASE_SET:
            raise ValueError(
                f"Invalid DNA base {base!r} at position {index}; expected only A, C, G, or T"
            )


def _validate_trits(trits: list[int]) -> None:
    if not isinstance(trits, list):
        raise TypeError("trits must be a list of integers")

    for index, trit in enumerate(trits):
        if type(trit) is not int:
            raise TypeError(f"Trit at position {index} must be an integer")
        if trit not in (0, 1, 2):
            raise ValueError(f"Trit at position {index} must be 0, 1, or 2")

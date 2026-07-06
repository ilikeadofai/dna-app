"""Reusable DNA sequence metrics.

These metrics are intended for an educational DNA data-storage simulator. In
particular, ``hairpin_like_score`` is a simplified proxy based on short
reverse-complement matches; it is not a thermodynamic secondary-structure
prediction model.
"""

from __future__ import annotations

from .constraints import COMPLEMENT, DNA_BASE_SET, DNA_BASES, FORBIDDEN_MOTIFS


def normalize_dna(seq: str) -> str:
    """Return uppercase DNA with whitespace removed, rejecting invalid bases."""
    normalized = "".join(seq.split()).upper()
    _raise_for_invalid_bases(normalized)
    return normalized


def validate_dna(seq: str) -> None:
    """Validate that a sequence contains only A, C, G, and T after normalization."""
    normalize_dna(seq)


def base_counts(seq: str) -> dict[str, int]:
    """Count A/C/G/T bases, always returning all four keys."""
    normalized = normalize_dna(seq)
    return {base: normalized.count(base) for base in DNA_BASES}


def gc_content(seq: str) -> float:
    """Return GC content as a ratio from 0.0 to 1.0."""
    normalized = normalize_dna(seq)
    if not normalized:
        return 0.0

    counts = base_counts(normalized)
    return (counts["G"] + counts["C"]) / len(normalized)


def max_homopolymer_run(seq: str) -> int:
    """Return the longest run of the same base."""
    normalized = normalize_dna(seq)
    if not normalized:
        return 0

    longest = 1
    current = 1
    previous = normalized[0]

    for base in normalized[1:]:
        if base == previous:
            current += 1
        else:
            longest = max(longest, current)
            current = 1
            previous = base

    return max(longest, current)


def local_gc_stats(seq: str, window: int = 30) -> dict[str, float | int | list[float]]:
    """Return sliding-window GC-content statistics.

    If ``window`` is larger than the sequence, the whole sequence is used as a
    single effective window. Empty input returns zero-like stats.
    """
    if window <= 0:
        raise ValueError("Local GC window must be greater than 0")

    normalized = normalize_dna(seq)
    if not normalized:
        return {
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "window": window,
            "effective_window": 0,
            "values": [],
        }

    effective_window = min(window, len(normalized))
    values = [
        _gc_ratio(normalized[index : index + effective_window])
        for index in range(0, len(normalized) - effective_window + 1)
    ]

    return {
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
        "window": window,
        "effective_window": effective_window,
        "values": values,
    }


def count_forbidden_motifs(
    seq: str, motifs: list[str] | tuple[str, ...] | None = None
) -> int:
    """Count forbidden motif occurrences, including overlapping matches."""
    normalized = normalize_dna(seq)
    motifs_to_count = FORBIDDEN_MOTIFS if motifs is None else motifs

    total = 0
    for motif in motifs_to_count:
        normalized_motif = normalize_dna(motif)
        if not normalized_motif:
            raise ValueError("Forbidden motifs cannot be empty")

        start = 0
        while True:
            match_index = normalized.find(normalized_motif, start)
            if match_index == -1:
                break
            total += 1
            start = match_index + 1

    return total


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    normalized = normalize_dna(seq)
    return "".join(COMPLEMENT[base] for base in reversed(normalized))


def hairpin_like_score(seq: str, min_stem: int = 4, max_stem: int = 8) -> int:
    """Score short reverse-complement stem-like matches.

    This is an educational proxy only. It counts non-overlapping downstream
    reverse-complement matches for stem sizes in the requested range.
    """
    if min_stem <= 0:
        raise ValueError("min_stem must be greater than 0")
    if max_stem < min_stem:
        raise ValueError("max_stem must be greater than or equal to min_stem")

    normalized = normalize_dna(seq)
    score = 0

    for stem_len in range(min_stem, max_stem + 1):
        if len(normalized) < stem_len * 2:
            continue

        for start in range(0, len(normalized) - stem_len + 1):
            stem = normalized[start : start + stem_len]
            stem_reverse_complement = _reverse_complement_validated(stem)
            search_start = start + stem_len

            while search_start <= len(normalized) - stem_len:
                match_index = normalized.find(stem_reverse_complement, search_start)
                if match_index == -1:
                    break
                score += 1
                search_start = match_index + 1

    return score


def repeat_score(
    seq: str, k_values: tuple[int, ...] = (2, 3), min_repeats: int = 3
) -> int:
    """Score tandem repeated k-mer patterns such as ATATAT or CAGCAGCAG."""
    if not k_values:
        raise ValueError("k_values must contain at least one k-mer size")
    if any(k <= 0 for k in k_values):
        raise ValueError("All k-mer sizes must be greater than 0")
    if min_repeats < 2:
        raise ValueError("min_repeats must be at least 2")

    normalized = normalize_dna(seq)
    score = 0

    for k in k_values:
        index = 0
        while index <= len(normalized) - (k * min_repeats):
            unit = normalized[index : index + k]
            repeats = 1
            next_index = index + k

            while normalized[next_index : next_index + k] == unit:
                repeats += 1
                next_index += k

            if repeats >= min_repeats:
                score += repeats - min_repeats + 1
                index = next_index
            else:
                index += 1

    return score


def bits_per_base(original_bits: int, dna_length: int) -> float:
    """Return information density as original bits divided by DNA length."""
    if original_bits < 0:
        raise ValueError("original_bits must be non-negative")
    if dna_length <= 0:
        raise ValueError("dna_length must be greater than 0")

    return original_bits / dna_length


def analyze_sequence(
    seq: str, original_bits: int | None = None, window: int = 30
) -> dict:
    """Return a combined metrics dictionary for a DNA sequence."""
    normalized = normalize_dna(seq)
    counts = base_counts(normalized)
    local_gc = local_gc_stats(normalized, window=window)
    dna_length = len(normalized)

    density = None
    if original_bits is not None:
        if dna_length == 0:
            if original_bits != 0:
                raise ValueError(
                    "Cannot calculate bits per base for an empty DNA sequence"
                )
            density = 0.0
        else:
            density = bits_per_base(original_bits, dna_length)

    return {
        "length": dna_length,
        "base_counts": counts,
        "gc_content": gc_content(normalized),
        "gc_percent": gc_content(normalized) * 100,
        "local_gc": local_gc,
        "max_homopolymer_run": max_homopolymer_run(normalized),
        "forbidden_motif_count": count_forbidden_motifs(normalized),
        "hairpin_like_score": hairpin_like_score(normalized),
        "repeat_score": repeat_score(normalized),
        "bits_per_base": density,
        "original_bits": original_bits,
    }


def _raise_for_invalid_bases(seq: str) -> None:
    for index, base in enumerate(seq):
        if base not in DNA_BASE_SET:
            raise ValueError(
                f"Invalid DNA base {base!r} at position {index}; expected only A, C, G, or T"
            )


def _gc_ratio(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _reverse_complement_validated(seq: str) -> str:
    return "".join(COMPLEMENT[base] for base in reversed(seq))

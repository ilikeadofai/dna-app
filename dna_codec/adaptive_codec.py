"""Adaptive optimized DNA codec.

The adaptive codec generates multiple reversible rotating-ternary candidates for
a byte payload, scores each candidate using DNA-friendly educational metrics,
and stores the selected parameters in a deterministic strand header.
"""

from __future__ import annotations

import hashlib
import zlib
from dataclasses import dataclass
from typing import cast

from .checksum import crc32_bytes, verify_crc32
from .chunking import (
    COMPACT_HEADER_DNA_BASES,
    DEFAULT_CHUNK_SIZE,
    StrandHeader,
    header_from_dna,
    header_to_dna,
    reassemble_chunks,
    split_chunks,
)
from .constraints import DNA_BASE_SET, DNA_BASES
from .metrics import (
    analyze_sequence,
    count_forbidden_motifs,
    gc_content,
    hairpin_like_score,
    local_gc_stats,
    max_homopolymer_run,
    normalize_dna,
    repeat_score,
)
from .ternary_codec import MAPPING_IDS, rotating_decode, rotating_encode

ADAPTIVE_FORMAT_VERSION = 1
ADAPTIVE_CODEC_ID = 3
ADAPTIVE_PREAMBLE = "ACGTACGTACGTACGT"
ADAPTIVE_HEADER_DNA_BASES = len(ADAPTIVE_PREAMBLE) + COMPACT_HEADER_DNA_BASES
MAX_SEED_TRIALS = 256
DEFAULT_SEED_TRIALS = 16
DEFAULT_LOCAL_GC_WINDOW = 30

DEFAULT_SCORE_WEIGHTS: dict[str, float] = {
    "global_gc": 100.0,
    "local_gc": 80.0,
    "homopolymer": 50.0,
    "forbidden_motif": 30.0,
    "hairpin": 10.0,
    "repeat": 10.0,
}


class AdaptiveCodecError(ValueError):
    """Base error for adaptive codec failures."""


class ChecksumMismatchError(AdaptiveCodecError):
    """Raised when decoded payload bytes fail CRC32 verification."""


@dataclass(frozen=True)
class AdaptiveCandidate:
    payload_dna: str
    seed: int
    mapping_id: int
    start_base: str
    score: dict


@dataclass(frozen=True)
class AdaptiveEncoded:
    dna: str
    payload_dna: str
    header: StrandHeader
    score: dict


@dataclass(frozen=True)
class AdaptiveDecoded:
    data: bytes
    header: StrandHeader


def prng_bytes(seed: int, length: int) -> bytes:
    """Return a deterministic hash-based byte stream for candidate whitening."""
    _validate_non_negative_int(seed, "seed")
    _validate_non_negative_int(length, "length")

    if length == 0:
        return b""

    output = bytearray()
    counter = 0
    while len(output) < length:
        block = hashlib.sha256(f"{seed}:{counter}".encode("ascii")).digest()
        output.extend(block)
        counter += 1

    return bytes(output[:length])


def xor_whiten(data: bytes, seed: int) -> bytes:
    """XOR bytes with a deterministic PRNG stream.

    This is not encryption. It only creates reversible candidate diversity for
    adaptive sequence selection.
    """
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")

    stream = prng_bytes(seed, len(data))
    return bytes(byte ^ mask for byte, mask in zip(data, stream))


def score_sequence(
    seq: str,
    window: int = DEFAULT_LOCAL_GC_WINDOW,
    weights: dict[str, float] | None = None,
) -> dict:
    """Score a DNA sequence using educational DNA-friendliness metrics."""
    if window <= 0:
        raise ValueError("window must be greater than 0")

    normalized = normalize_dna(seq)
    resolved_weights = _resolve_weights(weights)
    sequence_metrics = analyze_sequence(normalized, window=window)

    global_gc_penalty = abs(gc_content(normalized) - 0.5) if normalized else 0.0
    local_values = cast(
        list[float], local_gc_stats(normalized, window=window)["values"]
    )
    local_gc_penalty = _mean_distance_outside_range(local_values, lower=0.4, upper=0.6)
    homopolymer_penalty = float(max(0, max_homopolymer_run(normalized) - 1) ** 2)
    forbidden_motif_penalty = float(count_forbidden_motifs(normalized))
    hairpin_penalty = float(hairpin_like_score(normalized))
    repeat_penalty = float(repeat_score(normalized))

    total = (
        resolved_weights["global_gc"] * global_gc_penalty
        + resolved_weights["local_gc"] * local_gc_penalty
        + resolved_weights["homopolymer"] * homopolymer_penalty
        + resolved_weights["forbidden_motif"] * forbidden_motif_penalty
        + resolved_weights["hairpin"] * hairpin_penalty
        + resolved_weights["repeat"] * repeat_penalty
    )

    return {
        "total": total,
        "global_gc_penalty": global_gc_penalty,
        "local_gc_penalty": local_gc_penalty,
        "homopolymer_penalty": homopolymer_penalty,
        "forbidden_motif_penalty": forbidden_motif_penalty,
        "hairpin_penalty": hairpin_penalty,
        "repeat_penalty": repeat_penalty,
        "metrics": sequence_metrics,
        "weights": resolved_weights,
    }


def choose_best_candidate(
    data: bytes,
    seed_trials: int = DEFAULT_SEED_TRIALS,
    window: int = DEFAULT_LOCAL_GC_WINDOW,
    weights: dict[str, float] | None = None,
) -> AdaptiveCandidate:
    """Generate reversible candidates and return the lowest-scoring sequence."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")
    _validate_positive_int(seed_trials, "seed_trials")
    if seed_trials > MAX_SEED_TRIALS:
        raise ValueError(f"seed_trials must be less than or equal to {MAX_SEED_TRIALS}")
    if window <= 0:
        raise ValueError("window must be greater than 0")

    best: AdaptiveCandidate | None = None

    for seed in range(seed_trials):
        whitened = xor_whiten(data, seed)
        for mapping_id in MAPPING_IDS:
            for start_base in DNA_BASES:
                payload_dna = rotating_encode(
                    whitened,
                    start_base=start_base,
                    mapping_id=mapping_id,
                )
                score = score_sequence(payload_dna, window=window, weights=weights)
                candidate = AdaptiveCandidate(
                    payload_dna=payload_dna,
                    seed=seed,
                    mapping_id=mapping_id,
                    start_base=start_base,
                    score=score,
                )

                if best is None or candidate.score["total"] < best.score["total"]:
                    best = candidate

    if best is None:
        raise AdaptiveCodecError("No adaptive candidates were generated")

    return best


def adaptive_encode(
    data: bytes,
    *,
    seed_trials: int = DEFAULT_SEED_TRIALS,
    chunk_index: int = 0,
    total_chunks: int = 1,
    original_length: int | None = None,
    compressed: bool = False,
    window: int = DEFAULT_LOCAL_GC_WINDOW,
    weights: dict[str, float] | None = None,
) -> AdaptiveEncoded:
    """Encode one independently decodable adaptive DNA strand."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")
    if type(compressed) is not bool:
        raise TypeError("compressed must be a boolean")

    _validate_non_negative_int(chunk_index, "chunk_index")
    _validate_positive_int(total_chunks, "total_chunks")
    if chunk_index >= total_chunks:
        raise ValueError("chunk_index must be less than total_chunks")

    if original_length is None:
        original_length = len(data)
    _validate_non_negative_int(original_length, "original_length")

    payload_bytes = zlib.compress(data) if compressed else data
    best = choose_best_candidate(
        payload_bytes,
        seed_trials=seed_trials,
        window=window,
        weights=weights,
    )

    header = StrandHeader(
        version=ADAPTIVE_FORMAT_VERSION,
        codec_id=ADAPTIVE_CODEC_ID,
        chunk_index=chunk_index,
        total_chunks=total_chunks,
        original_length=original_length,
        payload_length=len(payload_bytes),
        seed=best.seed,
        mapping_id=best.mapping_id,
        start_base=best.start_base,
        compressed=compressed,
        crc32=crc32_bytes(payload_bytes),
    )
    dna = _frame_strand(header, best.payload_dna)

    return AdaptiveEncoded(
        dna=dna,
        payload_dna=best.payload_dna,
        header=header,
        score=best.score,
    )


def adaptive_decode(seq: str) -> AdaptiveDecoded:
    """Decode one adaptive DNA strand and verify its CRC32 checksum."""
    header, payload_dna = parse_adaptive_strand(seq)

    whitened_payload = rotating_decode(
        payload_dna,
        length=header.payload_length,
        start_base=header.start_base,
        mapping_id=header.mapping_id,
    )
    payload_bytes = xor_whiten(whitened_payload, header.seed)

    if not verify_crc32(payload_bytes, header.crc32):
        raise ChecksumMismatchError("CRC32 checksum mismatch for adaptive payload")

    if header.compressed:
        try:
            decoded_data = zlib.decompress(payload_bytes)
        except zlib.error as exc:
            raise AdaptiveCodecError(
                "Compressed adaptive payload could not be decompressed"
            ) from exc
    else:
        decoded_data = payload_bytes

    if header.total_chunks == 1 and len(decoded_data) != header.original_length:
        raise AdaptiveCodecError(
            "Decoded data length does not match header original_length"
        )

    return AdaptiveDecoded(data=decoded_data, header=header)


def parse_adaptive_strand(seq: str) -> tuple[StrandHeader, str]:
    """Parse a framed adaptive strand into header metadata and payload DNA."""
    _validate_strict_dna(seq)

    if not seq.startswith(ADAPTIVE_PREAMBLE):
        raise AdaptiveCodecError("Adaptive DNA strand has an invalid preamble")

    header_dna_start = len(ADAPTIVE_PREAMBLE)
    header_dna_end = header_dna_start + COMPACT_HEADER_DNA_BASES
    if len(seq) < header_dna_end:
        raise AdaptiveCodecError("Adaptive DNA strand has truncated compact header DNA")

    header_dna = seq[header_dna_start:header_dna_end]
    header = header_from_dna(header_dna)
    _validate_adaptive_header(header)

    return header, seq[header_dna_end:]


def bytes_to_adaptive_dna(data: bytes, **kwargs) -> str:
    """Convenience wrapper returning only the framed adaptive DNA string."""
    return adaptive_encode(data, **kwargs).dna


def adaptive_dna_to_bytes(seq: str) -> bytes:
    """Convenience wrapper returning only decoded bytes."""
    return adaptive_decode(seq).data


def text_to_adaptive_dna(text: str, **kwargs) -> str:
    """Encode UTF-8 text into an adaptive DNA strand."""
    if not isinstance(text, str):
        raise TypeError("text must be a string")
    return bytes_to_adaptive_dna(text.encode("utf-8"), **kwargs)


def adaptive_dna_to_text(seq: str) -> str:
    """Decode an adaptive DNA strand into UTF-8 text."""
    return adaptive_dna_to_bytes(seq).decode("utf-8")


def adaptive_encode_chunks(
    data: bytes,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    seed_trials: int = DEFAULT_SEED_TRIALS,
    compressed: bool = False,
    window: int = DEFAULT_LOCAL_GC_WINDOW,
    weights: dict[str, float] | None = None,
) -> list[AdaptiveEncoded]:
    """Encode bytes into multiple independently decodable adaptive strands."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")

    chunks = split_chunks(data, chunk_size=chunk_size)
    if not chunks:
        chunks = [b""]

    total_chunks = len(chunks)
    return [
        adaptive_encode(
            chunk,
            seed_trials=seed_trials,
            chunk_index=index,
            total_chunks=total_chunks,
            original_length=len(data),
            compressed=compressed,
            window=window,
            weights=weights,
        )
        for index, chunk in enumerate(chunks)
    ]


def adaptive_decode_chunks(strands: list[str]) -> bytes:
    """Decode and reassemble adaptive strands in any order."""
    if not isinstance(strands, list):
        raise TypeError("strands must be a list of DNA strings")
    if not strands:
        raise ValueError("At least one adaptive strand is required")

    decoded = [adaptive_decode(strand) for strand in strands]
    total_chunks = decoded[0].header.total_chunks
    original_length = decoded[0].header.original_length

    for item in decoded:
        if item.header.total_chunks != total_chunks:
            raise AdaptiveCodecError("Decoded strands disagree on total_chunks")
        if item.header.original_length != original_length:
            raise AdaptiveCodecError("Decoded strands disagree on original_length")

    reassembled = reassemble_chunks(
        [(item.header.chunk_index, item.data) for item in decoded],
        total_chunks=total_chunks,
    )
    if len(reassembled) != original_length:
        raise AdaptiveCodecError(
            "Reassembled data length does not match original_length"
        )

    return reassembled


def _frame_strand(header: StrandHeader, payload_dna: str) -> str:
    return ADAPTIVE_PREAMBLE + header_to_dna(header) + payload_dna


def _validate_adaptive_header(header: StrandHeader) -> None:
    if header.version != ADAPTIVE_FORMAT_VERSION:
        raise AdaptiveCodecError(
            f"Unsupported adaptive format version: {header.version}"
        )
    if header.codec_id != ADAPTIVE_CODEC_ID:
        raise AdaptiveCodecError(f"Unsupported adaptive codec id: {header.codec_id}")


def _resolve_weights(weights: dict[str, float] | None) -> dict[str, float]:
    resolved = dict(DEFAULT_SCORE_WEIGHTS)
    if weights is None:
        return resolved
    if not isinstance(weights, dict):
        raise TypeError("weights must be a dictionary")

    unknown = set(weights) - set(DEFAULT_SCORE_WEIGHTS)
    if unknown:
        raise ValueError(f"Unknown score weight(s): {sorted(unknown)}")

    for name, value in weights.items():
        if type(value) not in (int, float):
            raise TypeError(f"Weight {name!r} must be numeric")
        if value < 0:
            raise ValueError(f"Weight {name!r} must be non-negative")
        resolved[name] = float(value)

    return resolved


def _mean_distance_outside_range(
    values: list[float], lower: float, upper: float
) -> float:
    if not values:
        return 0.0

    distances = []
    for value in values:
        if value < lower:
            distances.append(lower - value)
        elif value > upper:
            distances.append(value - upper)
        else:
            distances.append(0.0)

    return sum(distances) / len(distances)


def _validate_strict_dna(seq: str) -> None:
    if not isinstance(seq, str):
        raise TypeError("seq must be a string")

    for index, base in enumerate(seq):
        if base not in DNA_BASE_SET:
            raise AdaptiveCodecError(
                f"Invalid DNA base {base!r} at position {index}; expected only A, C, G, or T"
            )


def _validate_non_negative_int(value: int, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if value < 0:
        raise ValueError(f"{name} must be non-negative")


def _validate_positive_int(value: int, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if value <= 0:
        raise ValueError(f"{name} must be greater than 0")

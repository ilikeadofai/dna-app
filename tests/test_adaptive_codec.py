import random
from dataclasses import replace

import pytest

from dna_codec.adaptive_codec import (
    ADAPTIVE_CODEC_ID,
    ADAPTIVE_FORMAT_VERSION,
    ADAPTIVE_HEADER_DNA_BASES,
    ADAPTIVE_PREAMBLE,
    AdaptiveCodecError,
    ChecksumMismatchError,
    adaptive_decode,
    adaptive_decode_chunks,
    adaptive_dna_to_bytes,
    adaptive_dna_to_text,
    adaptive_encode,
    adaptive_encode_chunks,
    bytes_to_adaptive_dna,
    choose_best_candidate,
    parse_adaptive_strand,
    prng_bytes,
    score_sequence,
    text_to_adaptive_dna,
    xor_whiten,
)
from dna_codec.checksum import crc32_bytes
from dna_codec.chunking import COMPACT_HEADER_DNA_BASES, header_to_dna
from dna_codec.ternary_codec import (
    get_candidate_bases,
    rotating_decode,
    rotating_encode,
)


def reframe(header, payload_dna):
    return ADAPTIVE_PREAMBLE + header_to_dna(header) + payload_dna


def mutate_payload_to_checksum_failure(encoded):
    payload_start = len(encoded.dna) - len(encoded.payload_dna)
    payload = list(encoded.payload_dna)

    for index, current_base in enumerate(payload):
        previous_base = encoded.header.start_base if index == 0 else payload[index - 1]
        candidates = get_candidate_bases(previous_base, encoded.header.mapping_id)

        for replacement in candidates:
            if replacement == current_base:
                continue
            if index + 1 < len(payload):
                next_base = payload[index + 1]
                if next_base not in get_candidate_bases(
                    replacement, encoded.header.mapping_id
                ):
                    continue

            mutated_payload = payload.copy()
            mutated_payload[index] = replacement
            mutated_dna = encoded.dna[:payload_start] + "".join(mutated_payload)

            try:
                adaptive_decode(mutated_dna)
            except ChecksumMismatchError:
                return mutated_dna
            except ValueError:
                continue

    raise AssertionError("Could not produce a valid-transition checksum mutation")


def test_prng_bytes_is_deterministic_and_seeded():
    assert prng_bytes(0, 0) == b""
    assert prng_bytes(7, 64) == prng_bytes(7, 64)
    assert prng_bytes(7, 64) != prng_bytes(8, 64)
    assert len(prng_bytes(7, 65)) == 65


def test_prng_bytes_rejects_invalid_params():
    with pytest.raises(ValueError):
        prng_bytes(-1, 8)

    with pytest.raises(ValueError):
        prng_bytes(1, -1)

    for bad_seed in [True, 1.5, "1"]:
        with pytest.raises(TypeError):
            prng_bytes(bad_seed, 8)

    for bad_length in [True, 1.5, "8"]:
        with pytest.raises(TypeError):
            prng_bytes(1, bad_length)


def test_xor_whiten_is_reversible():
    data = b"\x00\x00hello adaptive dna\xff"
    whitened = xor_whiten(data, seed=3)

    assert whitened != data
    assert xor_whiten(whitened, seed=3) == data
    assert xor_whiten(b"", seed=3) == b""

    with pytest.raises(TypeError):
        xor_whiten("not bytes", seed=3)


def test_score_sequence_returns_expected_keys_and_prefers_balanced_gc():
    balanced = score_sequence("ACGTACGT", window=4)
    all_a = score_sequence("AAAAAAAA", window=4)

    expected_keys = {
        "total",
        "global_gc_penalty",
        "local_gc_penalty",
        "homopolymer_penalty",
        "forbidden_motif_penalty",
        "hairpin_penalty",
        "repeat_penalty",
        "metrics",
        "weights",
    }
    assert set(balanced) == expected_keys
    assert balanced["total"] < all_a["total"]


def test_score_sequence_penalty_components_and_custom_weights():
    motif_score = score_sequence("GAATTC", window=6)
    repeat = score_sequence("ATATAT", window=6)
    zero_weighted = score_sequence(
        "AAAAAA",
        window=3,
        weights={
            "global_gc": 0,
            "local_gc": 0,
            "homopolymer": 0,
            "forbidden_motif": 0,
            "hairpin": 0,
            "repeat": 0,
        },
    )

    assert motif_score["forbidden_motif_penalty"] == 1
    assert repeat["repeat_penalty"] > 0
    assert zero_weighted["total"] == pytest.approx(0.0)

    with pytest.raises(ValueError):
        score_sequence("ACGT", window=0)

    with pytest.raises(ValueError):
        score_sequence("ACGN")

    with pytest.raises(ValueError):
        score_sequence("ACGT", weights={"unknown": 1.0})

    with pytest.raises(ValueError):
        score_sequence("ACGT", weights={"global_gc": -1.0})

    with pytest.raises(TypeError):
        score_sequence("ACGT", weights={"global_gc": True})


def test_choose_best_candidate_is_deterministic_and_decodable():
    data = b"hello dna"

    first = choose_best_candidate(data, seed_trials=2, window=8)
    second = choose_best_candidate(data, seed_trials=2, window=8)

    assert first == second
    assert first.seed in {0, 1}
    assert first.mapping_id in range(6)
    assert first.start_base in "ACGT"
    assert set(first.payload_dna).issubset(set("ACGT"))
    assert isinstance(first.score["total"], float)

    whitened = rotating_decode(
        first.payload_dna,
        length=len(data),
        start_base=first.start_base,
        mapping_id=first.mapping_id,
    )
    assert xor_whiten(whitened, first.seed) == data


def test_choose_best_candidate_rejects_invalid_params():
    with pytest.raises(TypeError):
        choose_best_candidate("not bytes")

    with pytest.raises(ValueError):
        choose_best_candidate(b"abc", seed_trials=0)

    with pytest.raises(ValueError):
        choose_best_candidate(b"abc", seed_trials=257)

    with pytest.raises(TypeError):
        choose_best_candidate(b"abc", seed_trials=True)

    with pytest.raises(ValueError):
        choose_best_candidate(b"abc", window=0)


def test_adaptive_roundtrip_korean_emoji_text_helpers():
    text = "안녕하세요 DNA 저장 테스트입니다 🧬"
    dna = text_to_adaptive_dna(text, seed_trials=2, window=10)

    assert dna.startswith(ADAPTIVE_PREAMBLE)
    assert adaptive_dna_to_text(dna) == text


def test_adaptive_roundtrip_random_bytes_leading_zeroes_and_all_values():
    payloads = [
        b"",
        b"\x00",
        b"\x00\x00abc",
        bytes(range(256)),
        random.Random(99).randbytes(64),
    ]

    for data in payloads:
        dna = bytes_to_adaptive_dna(data, seed_trials=1, window=12)
        assert adaptive_dna_to_bytes(dna) == data


def test_adaptive_roundtrip_compressed_payload():
    data = (b"ACGT" * 64) + "압축 테스트 🧬".encode("utf-8")
    encoded = adaptive_encode(data, seed_trials=1, compressed=True, window=12)
    decoded = adaptive_decode(encoded.dna)

    assert decoded.data == data
    assert decoded.header.compressed is True
    assert decoded.header.original_length == len(data)
    assert decoded.header.payload_length <= len(data)


def test_adaptive_payload_length_matches_rotating_for_uncompressed_single_chunk():
    data = "안녕하세요 DNA 저장 테스트입니다 🧬".encode("utf-8")
    encoded = adaptive_encode(data, seed_trials=2, window=10)
    reference_rotating = rotating_encode(data, start_base="A", mapping_id=0)

    assert len(encoded.payload_dna) == len(reference_rotating)
    assert len(encoded.dna) == len(encoded.payload_dna) + ADAPTIVE_HEADER_DNA_BASES


def test_adaptive_metadata_and_parse_strand():
    data = b"metadata payload"
    encoded = adaptive_encode(data, seed_trials=2, window=8)
    header, payload_dna = parse_adaptive_strand(encoded.dna)

    assert encoded.dna.startswith(ADAPTIVE_PREAMBLE)
    assert len(encoded.dna[:ADAPTIVE_HEADER_DNA_BASES]) == ADAPTIVE_HEADER_DNA_BASES
    assert (
        len(
            encoded.dna[
                len(ADAPTIVE_PREAMBLE) : len(ADAPTIVE_PREAMBLE)
                + COMPACT_HEADER_DNA_BASES
            ]
        )
        == COMPACT_HEADER_DNA_BASES
    )
    assert header == encoded.header
    assert payload_dna == encoded.payload_dna
    assert header.version == ADAPTIVE_FORMAT_VERSION
    assert header.codec_id == ADAPTIVE_CODEC_ID
    assert header.payload_length == len(data)
    assert header.original_length == len(data)
    assert header.crc32 == crc32_bytes(data)
    assert header.seed == encoded.header.seed
    assert header.mapping_id == encoded.header.mapping_id
    assert header.start_base == encoded.header.start_base
    assert encoded.score["metrics"]["length"] == len(encoded.payload_dna)


def test_adaptive_decode_detects_checksum_mutation():
    encoded = adaptive_encode(b"checksum mutation payload", seed_trials=2, window=8)
    mutated_dna = mutate_payload_to_checksum_failure(encoded)

    with pytest.raises(ChecksumMismatchError):
        adaptive_decode(mutated_dna)


def test_adaptive_decode_rejects_bad_preamble_and_truncated_fields():
    encoded = adaptive_encode(b"abc", seed_trials=1, window=4)

    with pytest.raises(AdaptiveCodecError, match="preamble"):
        adaptive_decode("T" + encoded.dna[1:])

    with pytest.raises(AdaptiveCodecError, match="truncated compact header"):
        adaptive_decode(ADAPTIVE_PREAMBLE + "ACGT")

    truncated_header = encoded.dna[: len(ADAPTIVE_PREAMBLE) + 4]
    with pytest.raises(AdaptiveCodecError, match="truncated compact header"):
        adaptive_decode(truncated_header)

    with pytest.raises(AdaptiveCodecError):
        adaptive_decode(encoded.dna.lower())


def test_adaptive_decode_rejects_unsupported_version_and_codec_id():
    encoded = adaptive_encode(b"abc", seed_trials=1, window=4)

    bad_version = replace(encoded.header, version=2)
    with pytest.raises(AdaptiveCodecError, match="version"):
        adaptive_decode(reframe(bad_version, encoded.payload_dna))

    bad_codec = replace(encoded.header, codec_id=4)
    with pytest.raises(AdaptiveCodecError, match="codec id"):
        adaptive_decode(reframe(bad_codec, encoded.payload_dna))


def test_adaptive_encode_rejects_invalid_params():
    with pytest.raises(TypeError):
        adaptive_encode("not bytes")

    with pytest.raises(TypeError):
        adaptive_encode(b"abc", compressed=1)

    with pytest.raises(ValueError):
        adaptive_encode(b"abc", chunk_index=1, total_chunks=1)

    with pytest.raises(ValueError):
        adaptive_encode(b"abc", original_length=-1)


def test_adaptive_chunks_roundtrip_out_of_order_and_empty():
    data = b"abcdefghijklmnopqrstuvwxyz0123456789"
    encoded_chunks = adaptive_encode_chunks(data, chunk_size=7, seed_trials=1, window=6)
    shuffled = [item.dna for item in reversed(encoded_chunks)]

    assert len(encoded_chunks) > 1
    assert adaptive_decode_chunks(shuffled) == data

    empty_chunks = adaptive_encode_chunks(b"", chunk_size=7, seed_trials=1, window=6)
    assert len(empty_chunks) == 1
    assert adaptive_decode_chunks([empty_chunks[0].dna]) == b""


def test_adaptive_decode_chunks_detects_missing_duplicate_and_bad_input():
    data = b"abcdefghijklmnopqrstuvwxyz"
    encoded_chunks = adaptive_encode_chunks(data, chunk_size=5, seed_trials=1, window=5)

    with pytest.raises(ValueError, match="missing"):
        adaptive_decode_chunks([item.dna for item in encoded_chunks[:-1]])

    duplicate = [encoded_chunks[0].dna, encoded_chunks[0].dna]
    with pytest.raises(ValueError, match="duplicate"):
        adaptive_decode_chunks(duplicate)

    with pytest.raises(TypeError):
        adaptive_decode_chunks("not a list")

    with pytest.raises(ValueError):
        adaptive_decode_chunks([])

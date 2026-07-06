import random

import pytest

from dna_codec.metrics import max_homopolymer_run
from dna_codec.ternary_codec import (
    bytes_to_trits,
    get_candidate_bases,
    rotating_decode,
    rotating_encode,
    trits_to_bytes,
)


def test_bytes_to_trits_empty_roundtrip():
    assert bytes_to_trits(b"") == []
    assert trits_to_bytes([], length=0) == b""


def test_trits_roundtrip_all_single_byte_values():
    for value in range(256):
        data = bytes([value])
        trits = bytes_to_trits(data)

        assert all(trit in {0, 1, 2} for trit in trits)
        assert trits_to_bytes(trits, length=1) == data


def test_trits_roundtrip_random_and_leading_zero_bytes():
    rng = random.Random(12345)
    payloads = [
        b"\x00",
        b"\x00\x00",
        b"\x00\x01",
        b"\x00\x00abc",
        bytes([0]) * 32,
    ]
    payloads.extend(rng.randbytes(size) for size in [1, 2, 3, 8, 31, 128])

    for data in payloads:
        assert trits_to_bytes(bytes_to_trits(data), length=len(data)) == data


def test_trits_to_bytes_rejects_invalid_trits_length_and_overflow():
    with pytest.raises(ValueError):
        trits_to_bytes([-1, 0, 0, 0, 0, 0], length=1)

    with pytest.raises(ValueError):
        trits_to_bytes([3, 0, 0, 0, 0, 0], length=1)

    with pytest.raises(TypeError):
        trits_to_bytes([0, 1, "2"], length=1)

    with pytest.raises(TypeError):
        trits_to_bytes([0, 1, 2.0], length=1)

    with pytest.raises(ValueError, match="Expected"):
        trits_to_bytes([0, 1, 2], length=1)

    # Six trits are required for one byte. This value is 256, just outside one-byte range.
    with pytest.raises(ValueError, match="does not fit"):
        trits_to_bytes([1, 0, 0, 1, 1, 1], length=1)


def test_trits_to_bytes_rejects_invalid_length():
    with pytest.raises(ValueError):
        trits_to_bytes([], length=-1)

    with pytest.raises(TypeError):
        trits_to_bytes([], length=1.5)

    with pytest.raises(TypeError):
        trits_to_bytes([], length=True)


def test_get_candidate_bases_mapping_zero_expected_order():
    assert get_candidate_bases("A", mapping_id=0) == ["C", "G", "T"]
    assert get_candidate_bases("C", mapping_id=0) == ["A", "G", "T"]
    assert get_candidate_bases("G", mapping_id=0) == ["A", "C", "T"]
    assert get_candidate_bases("T", mapping_id=0) == ["A", "C", "G"]


def test_get_candidate_bases_all_mapping_permutations():
    for previous_base in "ACGT":
        seen = set()
        for mapping_id in range(6):
            candidates = get_candidate_bases(previous_base, mapping_id=mapping_id)
            seen.add(tuple(candidates))

            assert len(candidates) == 3
            assert previous_base not in candidates
            assert len(set(candidates)) == 3
            assert set(candidates).issubset(set("ACGT"))

        assert len(seen) == 6


def test_get_candidate_bases_rejects_invalid_params():
    for invalid_base in ["N", "X", "", "AA", "a"]:
        with pytest.raises(ValueError):
            get_candidate_bases(invalid_base)

    for invalid_mapping in [-1, 6]:
        with pytest.raises(ValueError):
            get_candidate_bases("A", mapping_id=invalid_mapping)

    for invalid_mapping in [1.5, True, "0"]:
        with pytest.raises(TypeError):
            get_candidate_bases("A", mapping_id=invalid_mapping)


def test_rotating_empty_roundtrip():
    assert rotating_encode(b"") == ""
    assert rotating_decode("", length=0) == b""


def test_rotating_roundtrip_korean_emoji():
    data = "안녕하세요 DNA 저장 테스트입니다 🧬".encode("utf-8")
    dna = rotating_encode(data, start_base="A", mapping_id=0)

    assert rotating_decode(dna, length=len(data), start_base="A", mapping_id=0) == data
    assert max_homopolymer_run(dna) <= 1


@pytest.mark.parametrize("start_base", list("ACGT"))
@pytest.mark.parametrize("mapping_id", range(6))
def test_rotating_roundtrip_all_mapping_and_start_base_combinations(
    start_base, mapping_id
):
    rng = random.Random(1000 + mapping_id + ord(start_base))
    payloads = [
        b"",
        b"\x00",
        b"\x00\x00",
        b"\x00\x01",
        b"\x00\x00abc",
        bytes([0]) * 32,
        bytes([255]) * 64,
        bytes(range(32)),
    ]
    payloads.extend(rng.randbytes(size) for size in [1, 2, 3, 31, 128])

    for data in payloads:
        dna = rotating_encode(data, start_base=start_base, mapping_id=mapping_id)

        assert set(dna).issubset(set("ACGT"))
        if dna:
            assert dna[0] != start_base
            assert max_homopolymer_run(dna) <= 1
            assert all(left != right for left, right in zip(dna, dna[1:]))

        assert (
            rotating_decode(
                dna, length=len(data), start_base=start_base, mapping_id=mapping_id
            )
            == data
        )


def test_rotating_roundtrip_large_random_payload():
    data = random.Random(2024).randbytes(1024)
    dna = rotating_encode(data, start_base="T", mapping_id=5)

    assert rotating_decode(dna, length=len(data), start_base="T", mapping_id=5) == data
    assert max_homopolymer_run(dna) <= 1


def test_rotating_encode_rejects_invalid_params():
    with pytest.raises(TypeError):
        rotating_encode(bytearray(b"abc"))

    for invalid_start in ["N", "", "AA", "a"]:
        with pytest.raises(ValueError):
            rotating_encode(b"abc", start_base=invalid_start)

    with pytest.raises(ValueError):
        rotating_encode(b"abc", mapping_id=6)

    with pytest.raises(TypeError):
        rotating_encode(b"abc", mapping_id=True)


def test_rotating_decode_rejects_invalid_params_and_dna():
    dna = rotating_encode(b"abc")

    with pytest.raises(ValueError):
        rotating_decode(dna, length=-1)

    with pytest.raises(TypeError):
        rotating_decode(dna, length=1.5)

    for invalid_start in ["N", "", "AA", "a"]:
        with pytest.raises(ValueError):
            rotating_decode(dna, length=3, start_base=invalid_start)

    with pytest.raises(ValueError):
        rotating_decode(dna, length=3, mapping_id=6)

    with pytest.raises(TypeError):
        rotating_decode(dna, length=3, mapping_id=True)

    for invalid_seq in ["ACGN", "acgt", "ACGT\n"]:
        with pytest.raises(ValueError):
            rotating_decode(invalid_seq, length=1)


def test_rotating_decode_rejects_wrong_length_impossible_sequence_and_overflow():
    dna = rotating_encode(b"abc")

    with pytest.raises(ValueError, match="Expected rotating DNA length"):
        rotating_decode(dna, length=2)

    with pytest.raises(ValueError, match="cannot follow"):
        rotating_decode("AAAAAA", length=1, start_base="A", mapping_id=0)

    # Valid rotating transitions for mapping 0/start A, but trits decode to 256.
    overflow_dna = "GACGCG"
    with pytest.raises(ValueError, match="does not fit"):
        rotating_decode(overflow_dna, length=1, start_base="A", mapping_id=0)

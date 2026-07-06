import os

import pytest

from dna_codec.binary_codec import (
    bytes_to_fixed_dna,
    fixed_dna_to_bytes,
    fixed_dna_to_text,
    text_to_fixed_dna,
)


def test_fixed_known_mapping_single_byte():
    assert bytes_to_fixed_dna(b"\x1b") == "ACGT"
    assert fixed_dna_to_bytes("ACGT") == b"\x1b"


def test_fixed_known_mapping_extreme_bytes():
    assert bytes_to_fixed_dna(b"\x00") == "AAAA"
    assert bytes_to_fixed_dna(b"\xff") == "TTTT"
    assert fixed_dna_to_bytes("AAAATTTT") == b"\x00\xff"


def test_fixed_empty_roundtrip():
    assert bytes_to_fixed_dna(b"") == ""
    assert fixed_dna_to_bytes("") == b""
    assert text_to_fixed_dna("") == ""
    assert fixed_dna_to_text("") == ""


def test_fixed_dna_length_is_four_bases_per_byte():
    data = b"abc123"
    dna = bytes_to_fixed_dna(data)
    assert len(dna) == 4 * len(data)


def test_fixed_roundtrip_korean_emoji():
    text = "안녕하세요 DNA 저장 테스트입니다 🧬"
    assert fixed_dna_to_text(text_to_fixed_dna(text)) == text


def test_fixed_roundtrip_mixed_unicode_symbols():
    text = "Hello 한글 😀 🧬 © € 𝄞"
    assert fixed_dna_to_text(text_to_fixed_dna(text)) == text


def test_fixed_roundtrip_all_byte_values():
    data = bytes(range(256))
    assert fixed_dna_to_bytes(bytes_to_fixed_dna(data)) == data


def test_fixed_roundtrip_random_bytes():
    for size in [1, 2, 3, 16, 255, 1024]:
        data = os.urandom(size)
        assert fixed_dna_to_bytes(bytes_to_fixed_dna(data)) == data


def test_invalid_dna_character_rejected():
    with pytest.raises(ValueError, match="X"):
        fixed_dna_to_bytes("ACXT")


def test_lowercase_dna_rejected():
    with pytest.raises(ValueError):
        fixed_dna_to_bytes("acgt")


def test_whitespace_dna_rejected():
    with pytest.raises(ValueError):
        fixed_dna_to_bytes("ACGT\n")


def test_dna_length_must_be_multiple_of_four():
    with pytest.raises(ValueError, match="multiple of 4"):
        fixed_dna_to_bytes("ACG")


def test_fixed_dna_to_text_invalid_utf8_raises_unicode_decode_error():
    with pytest.raises(UnicodeDecodeError):
        fixed_dna_to_text("TTTT")

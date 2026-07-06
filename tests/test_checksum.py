import pytest

from dna_codec.checksum import MAX_CRC32, crc32_bytes, verify_crc32


def test_crc32_empty_bytes_is_zero():
    assert crc32_bytes(b"") == 0


def test_crc32_known_standard_vector():
    assert crc32_bytes(b"123456789") == 0xCBF43926


def test_crc32_returns_unsigned_32_bit_int():
    value = crc32_bytes(bytes(range(256)))

    assert isinstance(value, int)
    assert 0 <= value <= MAX_CRC32


def test_crc32_rejects_non_bytes():
    with pytest.raises(TypeError):
        crc32_bytes("not bytes")


def test_verify_crc32_accepts_matching_checksum():
    data = b"hello dna"
    assert verify_crc32(data, crc32_bytes(data)) is True


def test_verify_crc32_rejects_mutated_payload():
    data = b"hello dna"
    expected = crc32_bytes(data)

    assert verify_crc32(b"hello dnb", expected) is False


def test_verify_crc32_rejects_wrong_expected_value():
    assert verify_crc32(b"abc", 0) is False


def test_verify_crc32_rejects_invalid_expected_values():
    assert verify_crc32(b"abc", -1) is False
    assert verify_crc32(b"abc", MAX_CRC32 + 1) is False
    assert verify_crc32(b"abc", 1.5) is False
    assert verify_crc32(b"abc", True) is False

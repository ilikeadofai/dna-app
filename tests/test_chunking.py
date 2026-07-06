import json

import pytest

from dna_codec.checksum import crc32_bytes
from dna_codec.chunking import (
    COMPACT_HEADER_BYTES,
    COMPACT_HEADER_DNA_BASES,
    DEFAULT_CHUNK_SIZE,
    StrandHeader,
    deserialize_header,
    header_from_dna,
    header_to_dna,
    reassemble_chunks,
    serialize_header,
    split_chunks,
)


def make_header(**overrides):
    values = {
        "version": 1,
        "codec_id": 2,
        "chunk_index": 0,
        "total_chunks": 3,
        "original_length": 300,
        "payload_length": 128,
        "seed": 0,
        "mapping_id": 1,
        "start_base": "A",
        "compressed": False,
        "crc32": 0xCBF43926,
    }
    values.update(overrides)
    return StrandHeader(**values)


def test_default_chunk_size_is_512_bytes():
    assert DEFAULT_CHUNK_SIZE == 512


def test_split_chunks_empty_data():
    assert split_chunks(b"", chunk_size=128) == []


def test_split_chunks_smaller_than_chunk_size():
    assert split_chunks(b"abc", chunk_size=128) == [b"abc"]


def test_split_chunks_exact_boundary_no_empty_trailing_chunk():
    assert split_chunks(b"abcdef", chunk_size=3) == [b"abc", b"def"]


def test_split_chunks_with_remainder():
    assert split_chunks(b"abcdefg", chunk_size=3) == [b"abc", b"def", b"g"]


def test_split_chunks_chunk_size_one():
    assert split_chunks(b"abc", chunk_size=1) == [b"a", b"b", b"c"]


def test_split_chunks_rejects_invalid_input_and_chunk_size():
    with pytest.raises(TypeError):
        split_chunks("abc", chunk_size=1)

    for chunk_size in [0, -1]:
        with pytest.raises(ValueError):
            split_chunks(b"abc", chunk_size=chunk_size)

    for chunk_size in [1.5, True]:
        with pytest.raises(TypeError):
            split_chunks(b"abc", chunk_size=chunk_size)


def test_reassemble_chunks_in_order():
    chunks = [(0, b"abc"), (1, b"def"), (2, b"g")]
    assert reassemble_chunks(chunks, total_chunks=3) == b"abcdefg"


def test_reassemble_chunks_out_of_order():
    chunks = [(2, b"g"), (0, b"abc"), (1, b"def")]
    assert reassemble_chunks(chunks, total_chunks=3) == b"abcdefg"


def test_reassemble_chunks_empty_total_zero():
    assert reassemble_chunks([], total_chunks=0) == b""


def test_reassemble_chunks_missing_chunk_detection():
    chunks = [(0, b"abc"), (2, b"g")]

    with pytest.raises(ValueError, match="missing"):
        reassemble_chunks(chunks, total_chunks=3)


def test_reassemble_chunks_duplicate_chunk_detection():
    chunks = [(0, b"abc"), (1, b"def"), (1, b"DEF")]

    with pytest.raises(ValueError, match="duplicate"):
        reassemble_chunks(chunks, total_chunks=2)


def test_reassemble_chunks_rejects_invalid_index_and_total():
    with pytest.raises(ValueError):
        reassemble_chunks([(-1, b"bad")], total_chunks=1)

    with pytest.raises(ValueError):
        reassemble_chunks([(0, b"abc"), (1, b"extra")], total_chunks=1)

    with pytest.raises(ValueError):
        reassemble_chunks([], total_chunks=-1)

    with pytest.raises(TypeError):
        reassemble_chunks([], total_chunks=True)


def test_reassemble_chunks_rejects_invalid_chunk_items():
    with pytest.raises(TypeError):
        reassemble_chunks([[0, b"abc"]], total_chunks=1)

    with pytest.raises(TypeError):
        reassemble_chunks([(0, "abc")], total_chunks=1)

    with pytest.raises(TypeError):
        reassemble_chunks([(False, b"abc")], total_chunks=1)


def test_strand_header_serialization_roundtrip():
    header = make_header()

    encoded = serialize_header(header)
    decoded = deserialize_header(encoded)

    assert decoded == header


def test_strand_header_serialization_is_deterministic_and_compact_binary():
    header = make_header(
        chunk_index=1,
        total_chunks=2,
        seed=4,
        mapping_id=3,
        start_base="G",
        compressed=True,
        crc32=12345,
    )

    first = serialize_header(header)
    second = serialize_header(header)

    assert first == second
    assert len(first) == COMPACT_HEADER_BYTES
    assert deserialize_header(first) == header


def test_strand_header_fixed_dna_roundtrip():
    payload = b"chunk payload"
    header = make_header(
        total_chunks=1,
        payload_length=len(payload),
        original_length=len(payload),
        crc32=crc32_bytes(payload),
    )

    dna = header_to_dna(header)

    assert set(dna).issubset(set("ACGT"))
    assert len(dna) == COMPACT_HEADER_DNA_BASES
    assert header_from_dna(dna) == header


def test_strand_header_rejects_invalid_start_base():
    with pytest.raises(ValueError):
        make_header(start_base="X")

    with pytest.raises(ValueError):
        make_header(start_base="a")


def test_strand_header_rejects_invalid_chunk_index():
    with pytest.raises(ValueError):
        make_header(chunk_index=2, total_chunks=2)

    with pytest.raises(ValueError):
        make_header(chunk_index=-1)


def test_strand_header_rejects_invalid_crc32_range():
    with pytest.raises(ValueError):
        make_header(crc32=0x1_0000_0000)

    with pytest.raises(ValueError):
        make_header(crc32=-1)


def test_strand_header_rejects_invalid_integer_and_bool_fields():
    with pytest.raises(ValueError):
        make_header(version=0)

    with pytest.raises(TypeError):
        make_header(total_chunks=True)

    with pytest.raises(TypeError):
        make_header(compressed=0)

    with pytest.raises(TypeError):
        make_header(mapping_id=True)

    with pytest.raises(ValueError):
        make_header(mapping_id=6)


def test_debug_json_header_serialization_roundtrip_and_field_validation():
    header = make_header()
    debug_json = serialize_header(header, debug_json=True)

    assert deserialize_header(debug_json, debug_json=True) == header
    assert list(json.loads(debug_json).keys()) == sorted(json.loads(debug_json).keys())

    header_dict = json.loads(debug_json)
    header_dict.pop("seed")
    with pytest.raises(ValueError, match="missing"):
        deserialize_header(json.dumps(header_dict).encode("utf-8"), debug_json=True)

    header_dict = json.loads(debug_json)
    header_dict["extra"] = 1
    with pytest.raises(ValueError, match="unexpected"):
        deserialize_header(json.dumps(header_dict).encode("utf-8"), debug_json=True)


def test_deserialize_header_rejects_invalid_compact_data():
    with pytest.raises(ValueError, match="exactly"):
        deserialize_header(b"not-json")

    with pytest.raises(ValueError, match="exactly"):
        deserialize_header(b"[]")

    with pytest.raises(TypeError):
        deserialize_header("not bytes")


def test_serialize_header_rejects_non_header():
    with pytest.raises(TypeError):
        serialize_header({"version": 1})

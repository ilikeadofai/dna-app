"""Chunking and compact strand metadata helpers.

DNA storage workflows split payloads into separately addressable strands. Headers
are serialized as a compact fixed-width binary structure and encoded to DNA with
the fixed 2-bit codec.
"""

from __future__ import annotations

import json
import struct
from dataclasses import asdict, dataclass, fields
from typing import Any

from .binary_codec import bytes_to_fixed_dna, fixed_dna_to_bytes
from .checksum import MAX_CRC32
from .constraints import DNA_BASE_SET, DNA_BASES
from .ternary_codec import MAPPING_IDS

DEFAULT_CHUNK_SIZE = 512
HEADER_ENCODING = "utf-8"
COMPACT_HEADER_BYTES = 20
COMPACT_HEADER_DNA_BASES = COMPACT_HEADER_BYTES * 4
HEADER_FLAGS_COMPRESSED = 0b0000_0001
_HEADER_STRUCT = struct.Struct(">BBHHIHBBBBI")
_START_BASE_TO_CODE = {base: index for index, base in enumerate(DNA_BASES)}
_CODE_TO_START_BASE = {index: base for base, index in _START_BASE_TO_CODE.items()}
MAX_UINT8 = 0xFF
MAX_UINT16 = 0xFFFF
MAX_UINT32 = 0xFFFFFFFF


@dataclass(frozen=True)
class StrandHeader:
    version: int
    codec_id: int
    chunk_index: int
    total_chunks: int
    original_length: int
    payload_length: int
    seed: int
    mapping_id: int
    start_base: str
    compressed: bool
    crc32: int

    def __post_init__(self) -> None:
        _validate_uint_range(self.version, "version", minimum=1, maximum=MAX_UINT8)
        _validate_uint_range(self.codec_id, "codec_id", maximum=MAX_UINT8)
        _validate_uint_range(self.chunk_index, "chunk_index", maximum=MAX_UINT16)
        _validate_uint_range(
            self.total_chunks, "total_chunks", minimum=1, maximum=MAX_UINT16
        )
        if self.chunk_index >= self.total_chunks:
            raise ValueError("chunk_index must be less than total_chunks")

        _validate_uint_range(
            self.original_length, "original_length", maximum=MAX_UINT32
        )
        _validate_uint_range(self.payload_length, "payload_length", maximum=MAX_UINT16)
        _validate_uint_range(self.seed, "seed", maximum=MAX_UINT8)

        if type(self.mapping_id) is not int:
            raise TypeError("mapping_id must be an integer")
        if self.mapping_id not in MAPPING_IDS:
            raise ValueError("mapping_id must be between 0 and 5")

        if not isinstance(self.start_base, str):
            raise TypeError("start_base must be a string")
        if len(self.start_base) != 1 or self.start_base not in DNA_BASE_SET:
            raise ValueError("start_base must be one of A, C, G, or T")

        if type(self.compressed) is not bool:
            raise TypeError("compressed must be a boolean")

        _validate_uint_range(self.crc32, "crc32", maximum=MAX_CRC32)


def split_chunks(data: bytes, chunk_size: int = DEFAULT_CHUNK_SIZE) -> list[bytes]:
    """Split bytes into chunks without adding an empty trailing chunk."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")
    _validate_positive_int(chunk_size, "chunk_size")

    return [
        data[index : index + chunk_size] for index in range(0, len(data), chunk_size)
    ]


def reassemble_chunks(chunks: list[tuple[int, bytes]], total_chunks: int) -> bytes:
    """Reassemble chunk bytes by index, detecting missing or duplicate chunks."""
    _validate_non_negative_int(total_chunks, "total_chunks")

    ordered: dict[int, bytes] = {}
    for item in chunks:
        if not isinstance(item, tuple) or len(item) != 2:
            raise TypeError("chunks must contain (chunk_index, payload) tuples")

        chunk_index, payload = item
        _validate_non_negative_int(chunk_index, "chunk_index")
        if chunk_index >= total_chunks:
            raise ValueError("chunk_index must be less than total_chunks")
        if not isinstance(payload, bytes):
            raise TypeError("chunk payloads must be bytes")
        if chunk_index in ordered:
            raise ValueError(f"duplicate chunk index: {chunk_index}")

        ordered[chunk_index] = payload

    missing = [index for index in range(total_chunks) if index not in ordered]
    if missing:
        raise ValueError(f"missing chunk index(es): {missing}")

    return b"".join(ordered[index] for index in range(total_chunks))


def serialize_header(header: StrandHeader, *, debug_json: bool = False) -> bytes:
    """Serialize a strand header.

    By default this returns a compact 20-byte binary header. ``debug_json=True``
    is available only for human inspection and is not used by DNA packaging.
    """
    if not isinstance(header, StrandHeader):
        raise TypeError("header must be a StrandHeader")

    if debug_json:
        return json.dumps(asdict(header), sort_keys=True, separators=(",", ":")).encode(
            HEADER_ENCODING
        )

    flags = HEADER_FLAGS_COMPRESSED if header.compressed else 0
    return _HEADER_STRUCT.pack(
        header.version,
        header.codec_id,
        header.chunk_index,
        header.total_chunks,
        header.original_length,
        header.payload_length,
        header.seed,
        header.mapping_id,
        _START_BASE_TO_CODE[header.start_base],
        flags,
        header.crc32,
    )


def deserialize_header(data: bytes, *, debug_json: bool = False) -> StrandHeader:
    """Deserialize compact binary header bytes into a validated StrandHeader."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")

    if debug_json:
        return _deserialize_json_header(data)

    if len(data) != COMPACT_HEADER_BYTES:
        raise ValueError(
            f"Compact strand header must be exactly {COMPACT_HEADER_BYTES} bytes"
        )

    try:
        (
            version,
            codec_id,
            chunk_index,
            total_chunks,
            original_length,
            payload_length,
            seed,
            mapping_id,
            start_base_code,
            flags,
            crc32,
        ) = _HEADER_STRUCT.unpack(data)
    except struct.error as exc:
        raise ValueError("Invalid compact strand header") from exc

    if start_base_code not in _CODE_TO_START_BASE:
        raise ValueError("Invalid start_base code in compact strand header")
    unknown_flags = flags & ~HEADER_FLAGS_COMPRESSED
    if unknown_flags:
        raise ValueError("Invalid flags in compact strand header")

    return StrandHeader(
        version=version,
        codec_id=codec_id,
        chunk_index=chunk_index,
        total_chunks=total_chunks,
        original_length=original_length,
        payload_length=payload_length,
        seed=seed,
        mapping_id=mapping_id,
        start_base=_CODE_TO_START_BASE[start_base_code],
        compressed=bool(flags & HEADER_FLAGS_COMPRESSED),
        crc32=crc32,
    )


def header_to_dna(header: StrandHeader) -> str:
    """Encode compact header bytes as fixed 2-bit DNA."""
    return bytes_to_fixed_dna(serialize_header(header))


def header_from_dna(seq: str) -> StrandHeader:
    """Decode fixed 2-bit compact header DNA into a validated StrandHeader."""
    return deserialize_header(fixed_dna_to_bytes(seq))


def _deserialize_json_header(data: bytes) -> StrandHeader:
    try:
        decoded = json.loads(data.decode(HEADER_ENCODING))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("Invalid serialized strand header") from exc

    if not isinstance(decoded, dict):
        raise ValueError("Serialized strand header must be a JSON object")

    expected_fields = {field.name for field in fields(StrandHeader)}
    actual_fields = set(decoded)
    if actual_fields != expected_fields:
        missing = sorted(expected_fields - actual_fields)
        unexpected = sorted(actual_fields - expected_fields)
        raise ValueError(
            f"Invalid strand header fields; missing={missing}, unexpected={unexpected}"
        )

    return StrandHeader(**decoded)


def _validate_non_negative_int(value: Any, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if value < 0:
        raise ValueError(f"{name} must be non-negative")


def _validate_positive_int(value: Any, name: str) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if value <= 0:
        raise ValueError(f"{name} must be greater than 0")


def _validate_uint_range(
    value: Any, name: str, *, minimum: int = 0, maximum: int
) -> None:
    if type(value) is not int:
        raise TypeError(f"{name} must be an integer")
    if value < minimum or value > maximum:
        if minimum == 0:
            raise ValueError(f"{name} must be between 0 and {maximum}")
        raise ValueError(f"{name} must be between {minimum} and {maximum}")

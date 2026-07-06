"""CRC32 checksum helpers for DNA storage chunks."""

from __future__ import annotations

import zlib

MAX_CRC32 = 0xFFFFFFFF


def crc32_bytes(data: bytes) -> int:
    """Return an unsigned CRC32 checksum for bytes."""
    if not isinstance(data, bytes):
        raise TypeError("data must be bytes")

    return zlib.crc32(data) & MAX_CRC32


def verify_crc32(data: bytes, expected: int) -> bool:
    """Return True when bytes match an expected unsigned CRC32 value."""
    if type(expected) is not int:
        return False
    if expected < 0 or expected > MAX_CRC32:
        return False

    return crc32_bytes(data) == expected

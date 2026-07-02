# Implementation Plan

## 1. Current repository state

The repository currently has a Streamlit app in `web.py`.

The existing app:

- supports text input and file upload
- extracts text from txt, md, docx, and pdf files
- encodes text into DNA using a lookup table
- requests missing characters from an external API
- decodes DNA using a reverse lookup table

The refactor should remove external API and lookup-table dependency and replace them with local reversible codecs.

## 2. Target structure

Recommended structure:

```text
dna-app/
├─ web.py
├─ dna_codec/
│  ├─ __init__.py
│  ├─ binary_codec.py
│  ├─ ternary_codec.py
│  ├─ adaptive_codec.py
│  ├─ metrics.py
│  ├─ constraints.py
│  ├─ chunking.py
│  └─ checksum.py
├─ tests/
│  ├─ test_binary_codec.py
│  ├─ test_ternary_codec.py
│  ├─ test_adaptive_codec.py
│  ├─ test_chunking.py
│  └─ test_metrics.py
├─ requirements.txt
└─ README.md
```

If the project needs to stay simpler, `dna_codec.py` may be used as a single module first, then split later.

## 3. Phase 1: fixed 2-bit local codec

### Goal

Replace external API and lookup table encoding with local UTF-8 based reversible encoding.

### Required functions

```python
BIT_TO_BASE = {"00": "A", "01": "C", "10": "G", "11": "T"}
BASE_TO_BIT = {"A": "00", "C": "01", "G": "10", "T": "11"}

def bytes_to_fixed_dna(data: bytes) -> str: ...
def fixed_dna_to_bytes(seq: str) -> bytes: ...
def text_to_fixed_dna(text: str) -> str: ...
def fixed_dna_to_text(seq: str) -> str: ...
```

### Rules

- Use UTF-8.
- Support Korean, emoji, symbols, and arbitrary Unicode.
- Reject invalid DNA characters.
- Fixed DNA length should be `4 * len(data)` because one byte becomes four DNA bases.

### Acceptance tests

- `"안녕하세요 DNA 저장 테스트입니다 🧬"` round-trips successfully.
- Random bytes round-trip successfully.
- Invalid DNA characters raise a clear error.

## 4. Phase 2: sequence metrics

### Goal

Implement reusable sequence analysis functions.

### Required functions

```python
def normalize_dna(seq: str) -> str: ...
def validate_dna(seq: str) -> None: ...
def base_counts(seq: str) -> dict[str, int]: ...
def gc_content(seq: str) -> float: ...
def max_homopolymer_run(seq: str) -> int: ...
def local_gc_stats(seq: str, window: int = 30) -> dict[str, float]: ...
def count_forbidden_motifs(seq: str, motifs: list[str] | None = None) -> int: ...
def reverse_complement(seq: str) -> str: ...
def hairpin_like_score(seq: str, min_stem: int = 4, max_stem: int = 8) -> int: ...
def repeat_score(seq: str, k_values: tuple[int, ...] = (2, 3), min_repeats: int = 3) -> int: ...
def bits_per_base(original_bits: int, dna_length: int) -> float: ...
def analyze_sequence(seq: str, original_bits: int | None = None, window: int = 30) -> dict: ...
```

### Notes

- `gc_content` should return a ratio from 0.0 to 1.0.
- UI may display it as a percentage.
- Local GC should return at least min, max, mean, and window size.
- Hairpin scoring is simplified and educational. Do not claim it is a full secondary-structure predictor.

## 5. Phase 3: rotating ternary codec

### Goal

Create a reversible DNA codec that prevents homopolymer repetition by converting bytes into trits and mapping each trit to a base that is not equal to the previous base.

### Core idea

At each step, remove the previous base from A/C/G/T, leaving three possible bases. A trit value 0, 1, or 2 selects one of those bases.

Example:

```text
previous A -> candidates C, G, T
previous C -> candidates A, G, T
previous G -> candidates A, C, T
previous T -> candidates A, C, G
```

### Required functions

```python
def bytes_to_trits(data: bytes) -> list[int]: ...
def trits_to_bytes(trits: list[int], length: int) -> bytes: ...
def get_candidate_bases(previous_base: str, mapping_id: int = 0) -> list[str]: ...
def rotating_encode(data: bytes, start_base: str = "A", mapping_id: int = 0) -> str: ...
def rotating_decode(seq: str, length: int, start_base: str = "A", mapping_id: int = 0) -> bytes: ...
```

### Important detail

Converting bytes to a big integer can remove leading zero bytes. Therefore, the original byte length must be stored in metadata and supplied to decoding.

### Acceptance tests

- Random bytes round-trip successfully.
- Korean text round-trips successfully.
- `max_homopolymer_run(rotating_encode(data)) <= 1`.

## 6. Phase 4: checksum and chunking

### Goal

Make encoded DNA strands independently decodable and able to detect errors.

### Chunk design

Default chunk size:

```python
DEFAULT_CHUNK_SIZE = 128
```

Each chunk should have metadata:

```python
@dataclass
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
```

### Required functions

```python
def crc32_bytes(data: bytes) -> int: ...
def verify_crc32(data: bytes, expected: int) -> bool: ...
def split_chunks(data: bytes, chunk_size: int = 128) -> list[bytes]: ...
def reassemble_chunks(chunks: list[tuple[int, bytes]], total_chunks: int) -> bytes: ...
```

### Header encoding

For simplicity and robustness, header may be encoded using the fixed 2-bit codec with a fixed byte length or a compact JSON-like metadata prefix.

Recommended educational approach:

- Use a clear textual header format first.
- Encode header bytes with fixed 2-bit codec.
- Add a fixed preamble to identify strand boundaries.
- Keep payload separately encoded by the selected codec.

Example logical structure:

```text
[Preamble][Header length][Header DNA][Payload DNA]
```

Avoid overengineering. The key requirement is reversible decoding and clear metadata.

## 7. Phase 5: adaptive optimized codec

### Goal

Generate multiple reversible DNA candidates and choose the most DNA-friendly one based on metrics.

### Candidate parameters

Try combinations of:

- `seed` from 0 to `seed_trials - 1`
- `mapping_id` from 0 to 5
- `start_base` from A/C/G/T

### Deterministic whitening

Whitening should create different candidate payloads without losing reversibility.

Recommended approach:

```python
def prng_bytes(seed: int, length: int) -> bytes: ...
def xor_whiten(data: bytes, seed: int) -> bytes: ...
```

Use a deterministic PRNG from Python standard library or hash-based stream generation. This is not encryption; it is only used to create candidate diversity.

Decoding must apply the same XOR again to recover the original payload.

### Candidate scoring

Required function:

```python
def score_sequence(seq: str, window: int = 30, weights: dict | None = None) -> dict:
    """Return score and sub-scores for a DNA sequence."""
```

Recommended sub-scores:

- global GC distance from 0.5
- local GC deviation beyond target range
- homopolymer penalty
- forbidden motif penalty
- hairpin-like score
- repeat score

### Adaptive encode pseudocode

```python
best = None

for seed in range(seed_trials):
    candidate_bytes = xor_whiten(payload, seed)

    for mapping_id in range(6):
        for start_base in "ACGT":
            dna = rotating_encode(candidate_bytes, start_base=start_base, mapping_id=mapping_id)
            score = score_sequence(dna)

            if best is None or score["total"] < best.score["total"]:
                best = Candidate(dna=dna, seed=seed, mapping_id=mapping_id, start_base=start_base, score=score)
```

### Acceptance tests

- Adaptive codec round-trips Korean and emoji text.
- Adaptive codec round-trips random bytes.
- If a base in payload is mutated, checksum should fail.
- Adaptive mode should produce metrics and selected parameters.

## 8. Phase 6: Streamlit UI

The existing `web.py` should be refactored, not discarded entirely. Keep useful features such as file upload and text extraction.

### Tabs

Use three main tabs:

1. Encode
2. Decode
3. Analyze DNA

### Encode tab

Controls:

- Direct text input
- File upload for txt/md/docx/pdf
- Codec selector:
  - Fixed 2-bit
  - Rotating ternary
  - Adaptive optimized
- Compression toggle
- Chunk size
- Seed trials
- Local GC window

Output:

- DNA result text area
- Download button
- Metrics table
- Base composition chart
- Comparison panel that runs all three codecs on the same input

### Decode tab

Inputs:

- DNA text area
- DNA txt file upload

Output:

- Decoded text
- Chunk status
- Checksum status
- Decode errors, if any

### Analyze tab

Input:

- Arbitrary DNA sequence

Output:

- Base counts
- GC content
- Local GC min/max
- Max homopolymer run
- Forbidden motifs
- Hairpin-like score
- Repeat score

## 9. Phase 7: tests

Use `pytest`.

Required tests:

```python
def test_fixed_roundtrip_korean_emoji(): ...
def test_fixed_roundtrip_random_bytes(): ...
def test_invalid_dna_rejected(): ...
def test_rotating_roundtrip_random_bytes(): ...
def test_rotating_no_homopolymer(): ...
def test_adaptive_roundtrip_korean_emoji(): ...
def test_crc_detects_mutation(): ...
def test_chunk_reassembly_out_of_order(): ...
def test_missing_chunk_detection(): ...
def test_metrics_gc_content(): ...
def test_forbidden_motif_counter(): ...
def test_repeat_score(): ...
```

## 10. README requirements

Update README with:

- project purpose
- how this differs from the old lookup/API version
- explanation of three encoding modes
- biological constraints considered
- how to run Streamlit
- how to run tests
- limitations

## 11. Development order

Recommended order:

1. Implement fixed 2-bit codec.
2. Implement metrics.
3. Implement rotating ternary codec.
4. Add unit tests for the above.
5. Add chunking and CRC32.
6. Implement adaptive candidate generation and scoring.
7. Refactor Streamlit UI.
8. Add comparison panel and CSV export.
9. Update README.
10. Run all tests.

Do not start with the adaptive optimizer before the simpler codecs pass tests.

## 12. Definition of done

The implementation is done when:

- No external API is required.
- No lookup JSON file is required for core encoding/decoding.
- `pytest` passes.
- `streamlit run web.py` works.
- Korean, English, symbols, and emoji round-trip successfully.
- The app compares all three encoding modes.
- The adaptive mode stores all required reversible parameters.
- The app clearly labels itself as an educational simulator.

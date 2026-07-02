# Codex Prompt

Use this prompt after the project documentation files have been added to the repository.

---

## Prompt to send to Codex

You are working in the repository `ilikeadofai/dna-app`.

Read these documentation files first:

- `docs/PROJECT_OVERVIEW.md`
- `docs/IMPLEMENTATION_PLAN.md`
- `docs/BIOLOGY_REPORT_NOTES.md`

Then implement the project according to the documentation.

## Goal

Refactor the current Streamlit DNA app from an external API / lookup-table based character encoder into a local, reversible, educational DNA data storage simulator.

The app must compare three encoding modes:

1. Fixed 2-bit UTF-8 DNA encoding
2. Rotating ternary DNA encoding that avoids homopolymers
3. Adaptive optimized DNA encoding that generates multiple reversible candidates and selects the best one using DNA-friendly constraints

## Existing code

The current core app is `web.py`. It currently uses:

- Streamlit UI
- text input
- file upload
- PDF/docx/text extraction
- external API lookup for unsupported characters
- local JSON lookup tables

Remove the external API and lookup-table dependency from the core encoding/decoding logic.

Keep useful UI features such as:

- direct text input
- txt/md/docx/pdf file upload if dependencies are available
- DNA output text area
- download buttons
- decoding from DNA input or DNA txt upload

## Required file structure

Create this structure if reasonable:

```text
dna_codec/
  __init__.py
  binary_codec.py
  ternary_codec.py
  adaptive_codec.py
  metrics.py
  constraints.py
  chunking.py
  checksum.py

tests/
  test_binary_codec.py
  test_ternary_codec.py
  test_adaptive_codec.py
  test_chunking.py
  test_metrics.py
```

Update:

- `web.py`
- `requirements.txt`, if needed
- `README.md`

## Functional requirements

### A. Fixed 2-bit codec

Implement a local UTF-8 codec:

```text
00 -> A
01 -> C
10 -> G
11 -> T
```

Requirements:

- Encode arbitrary UTF-8 text into DNA.
- Decode DNA back exactly.
- Support Korean, English, numbers, symbols, and emoji.
- Support arbitrary bytes internally.
- Reject invalid DNA characters clearly.
- No external API.
- No lookup JSON file required.

### B. Rotating ternary codec

Implement a codec that avoids homopolymers:

- Convert bytes into base-3 trits.
- At each DNA position, remove the previous base from A/C/G/T.
- Select one of the remaining three bases using the trit value.
- Support `mapping_id` permutations of trit-to-candidate mapping.
- Support `start_base` in A/C/G/T.
- Decode back exactly using original byte length, mapping id, and start base.

Important:

- Preserve leading zero bytes by storing and using original byte length.
- `max_homopolymer_run(rotating_encode(data))` should be <= 1.

### C. Sequence metrics

Implement reusable metrics:

```python
normalize_dna(seq)
validate_dna(seq)
base_counts(seq)
gc_content(seq)
max_homopolymer_run(seq)
local_gc_stats(seq, window=30)
count_forbidden_motifs(seq, motifs=None)
reverse_complement(seq)
hairpin_like_score(seq, min_stem=4, max_stem=8)
repeat_score(seq, k_values=(2, 3), min_repeats=3)
bits_per_base(original_bits, dna_length)
analyze_sequence(seq, original_bits=None, window=30)
```

The hairpin score is a simplified educational proxy, not a full secondary-structure predictor.

### D. Chunking and checksum

Implement chunking and CRC32 error detection.

Requirements:

- Split input bytes into chunks, default chunk size 128 bytes.
- Each chunk should be independently decodable.
- Reassemble chunks by chunk index.
- Detect missing chunks.
- Use CRC32 to detect mutated payloads.

Recommended metadata:

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

The exact header serialization can be simple, but it must be deterministic and decodable.

### E. Adaptive optimized codec

Implement adaptive encoding by generating multiple reversible DNA candidates for each chunk and selecting the lowest-score candidate.

Try candidate combinations:

- `seed` from 0 to `seed_trials - 1`
- `mapping_id` from 0 to 5
- `start_base` from A/C/G/T

Use deterministic XOR whitening to create different candidates:

- The whitening must be reversible.
- Store selected seed in metadata.
- This is not encryption; it only creates candidate diversity.

Score candidates using:

- global GC-content distance from 50%
- local GC-content deviation using sliding window
- max homopolymer run
- forbidden motif count
- simplified reverse-complement / hairpin-like score
- repeated k-mer score

Return selected metadata and sub-scores so the UI can display why the candidate was selected.

### F. Streamlit UI

Refactor `web.py` into three tabs:

1. Encode
2. Decode
3. Analyze DNA

#### Encode tab

Controls:

- direct text input
- file upload for txt/md/docx/pdf if available
- codec selector:
  - Fixed 2-bit
  - Rotating ternary
  - Adaptive optimized
- compression toggle
- chunk size
- seed trials
- local GC window

Outputs:

- DNA output
- download DNA file
- sequence metrics
- base composition chart
- comparison panel that runs the same input through all three codecs
- CSV download for the comparison table

#### Decode tab

Inputs:

- DNA text area
- DNA txt file upload

Outputs:

- decoded text
- checksum status
- chunk status
- clear error messages if decoding fails

#### Analyze DNA tab

Inputs:

- arbitrary DNA sequence

Outputs:

- base counts
- GC content
- local GC min/max
- max homopolymer run
- forbidden motif count
- hairpin-like score
- repeat score

### G. Tests

Add pytest tests:

- fixed 2-bit roundtrip for Korean + emoji
- fixed 2-bit roundtrip for random bytes
- invalid DNA characters rejected
- rotating ternary roundtrip for random bytes
- rotating ternary max homopolymer run <= 1
- adaptive roundtrip for Korean + emoji
- adaptive roundtrip for random bytes
- CRC detects mutation
- chunk reassembly works out of order
- missing chunk detection works
- metrics return expected values
- forbidden motif counter works
- repeat score works

### H. README

Update README with:

- project purpose
- how this differs from the old API/lookup version
- explanation of three encoding modes
- DNA constraints considered
- how to run Streamlit
- how to run tests
- limitations

## Development order

Do not start with the adaptive optimizer first. Use this order:

1. Implement fixed 2-bit codec.
2. Implement metrics.
3. Implement rotating ternary codec.
4. Add tests for the above.
5. Add chunking and CRC32.
6. Implement adaptive optimizer.
7. Refactor Streamlit UI.
8. Add comparison table and CSV export.
9. Update README.
10. Run tests and fix failures.

## Acceptance criteria

The task is complete when:

- `pytest` passes.
- `streamlit run web.py` works.
- The app no longer needs an external API.
- The app no longer needs lookup JSON files for core encoding/decoding.
- Korean and emoji text can be encoded and decoded.
- The three modes can be compared in the UI.
- Adaptive mode displays selected seed, mapping id, start base, metrics, and score.
- Mutated DNA can be detected by checksum in modes that use metadata/checksum.
- README clearly states that this is an educational simulator, not a production DNA storage or diagnostic tool.

## Style requirements

- Keep code readable and typed where reasonable.
- Prefer small pure functions for codecs and metrics.
- Keep Streamlit UI code separate from core encoding logic.
- Use clear Korean labels in the UI where appropriate.
- Do not make unsupported biological claims.
- Any simplified biological metric must be labeled as simplified or educational.

# DNA 데이터 저장 시뮬레이터

A Streamlit-based educational simulator for encoding digital data as DNA sequences and comparing different DNA storage strategies.

This project started as a lookup-table/API-based character-to-DNA converter. It has been refactored into a fully local, reversible DNA data storage simulator that works with UTF-8 text, Korean, English, symbols, emoji, and arbitrary bytes.

> This is an educational simulator. It is not a production DNA storage system, not a DNA synthesis tool, and not a medical or diagnostic tool.

## What changed from the old version

The old app depended on:

- a character lookup table
- `prebuilt_lookup.json`
- an external API for missing characters

The current app does **not** need any lookup JSON files or external encoding API. Text is converted locally using this pipeline:

```text
Unicode text -> UTF-8 bytes -> DNA sequence
```

## Encoding modes

### 1. Fixed 2-bit encoding

Baseline codec using a direct 2-bit mapping:

```text
00 -> A
01 -> C
10 -> G
11 -> T
```

Properties:

- shortest simple representation
- exactly 4 DNA bases per byte
- best raw density among the implemented modes
- does not avoid DNA-unfriendly patterns like long homopolymers

### 2. Rotating ternary encoding

Converts bytes into base-3 trits. At each step, it chooses from the three bases that are **not** equal to the previous base.

Properties:

- avoids adjacent repeated bases by construction
- `max_homopolymer_run` should be `1`
- longer than Fixed 2-bit
- requires metadata for decoding:
  - original byte length
  - start base
  - mapping id

### 3. Adaptive optimized encoding

Generates many reversible rotating-ternary candidates and selects the best-scoring one using DNA-friendliness metrics.

Candidate parameters:

- deterministic whitening seed
- ternary mapping id
- start base

Adaptive mode stores compact metadata and checksum information in each strand:

```text
[preamble][compact binary header encoded as DNA][payload DNA]
```

The compact binary header is 20 bytes, encoded as 80 DNA bases with the Fixed 2-bit codec. With the 16-base preamble, the current per-strand non-payload overhead is 96 DNA bases.

Adaptive mode is designed to improve DNA-friendliness, not to beat Fixed 2-bit raw density. The default adaptive candidate search uses 32 seed trials, each combined with 6 mapping IDs and 4 start bases.

## DNA constraints analyzed

The app calculates educational DNA metrics including:

- total DNA length
- A/C/G/T base counts
- GC content
- local GC min/mean/max
- max homopolymer run
- forbidden motif count
- simplified hairpin-like score
- repeat score
- bits per base

The forbidden motif list currently includes examples such as:

```text
GAATTC  EcoRI
AAGCTT  HindIII
GGATCC  BamHI
```

The hairpin-like score is a simplified educational proxy based on short reverse-complement patterns. It is **not** a real thermodynamic secondary-structure prediction.

## Fair comparison table

The Streamlit comparison table separates payload efficiency from packaging overhead:

- `payload_dna_length`
- `header_dna_length`
- `total_dna_length`
- `bits_per_base_payload`
- `bits_per_base_total`
- `gc_percent`
- `local_gc_min`
- `local_gc_max`
- `max_homopolymer`
- `forbidden_motif_count`
- `hairpin_like_score`
- `repeat_score`

This avoids unfairly comparing bare Fixed/Rotating payloads against fully packaged Adaptive strands.

Expected density behavior:

- Fixed 2-bit usually has the best raw density.
- Rotating ternary is longer but avoids homopolymers.
- Adaptive payload length is close to Rotating ternary payload length for uncompressed data.
- Adaptive total length includes compact metadata and CRC32 checksum overhead.

## Project structure

```text
dna_codec/
  __init__.py
  adaptive_codec.py
  binary_codec.py
  checksum.py
  chunking.py
  constraints.py
  metrics.py
  ternary_codec.py

tests/
  test_adaptive_codec.py
  test_binary_codec.py
  test_checksum.py
  test_chunking.py
  test_metrics.py
  test_ternary_codec.py

web.py
requirements.txt
```

## Installation

Create or activate a Python virtual environment, then install dependencies:

```sh
python -m pip install -r requirements.txt
```

If you are using the existing local virtual environment in this repository:

```sh
./.venv/bin/python -m pip install -r requirements.txt
```

## Run the app

```sh
streamlit run web.py
```

Or with the local virtual environment:

```sh
./.venv/bin/python -m streamlit run web.py
```

## Run tests

```sh
pytest
```

Or with the local virtual environment:

```sh
./.venv/bin/python -m pytest tests
```

## Continuous integration

GitHub Actions runs the test suite automatically on pushes and pull requests to `main` using `.github/workflows/tests.yml`.

## Usage notes

### Encode tab

Use this tab to:

- type text directly
- upload `.txt`, `.md`, `.docx`, or `.pdf` files
- choose an encoding mode
- download DNA output
- compare all three codecs

For Adaptive mode, the output may contain multiple newline-separated DNA strands if chunking is needed.

### Decode tab

Use this tab to decode DNA back into text.

- Fixed 2-bit requires no external metadata.
- Rotating ternary requires the original byte length, start base, and mapping id.
- Adaptive optimized is self-describing and verifies CRC32 checksums.

### Analyze tab

Use this tab to analyze any DNA sequence using the metrics listed above.

## Limitations

- No real DNA synthesis is performed.
- No real sequencing simulation is performed.
- CRC32 detects errors but does not correct them.
- There is no Reed-Solomon or other forward error correction yet.
- Primer design is not implemented.
- Hairpin scoring is simplified and educational.
- Forbidden motifs are examples, not a complete synthetic biology constraint set.
- Compression and whitening are not encryption.
- Adaptive optimization improves selected metrics but does not guarantee biological suitability.

## Suggested future work

- Add error-correction codes such as Reed-Solomon.
- Add sequencing error simulation.
- Add primer design constraints.
- Add downloadable experiment reports for biology class use.
- Add visual explanations of why a candidate sequence was selected.

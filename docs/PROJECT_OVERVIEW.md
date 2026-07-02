# Adaptive DNA-Aware Encoding Project Overview

## 1. Project purpose

This project refactors the existing Streamlit DNA app into an educational DNA data storage simulator.

The current app encodes text into DNA using a character lookup table and an external API. The new project should replace that approach with local, reversible, deterministic encoding methods that work for arbitrary Unicode text. The main goal is not only to convert text to DNA, but to compare several DNA encoding strategies by considering the biological and engineering constraints of DNA storage.

The central research question is:

> DNA uses four bases, A, T, G, and C, so it can represent digital information. However, real DNA sequences are affected by constraints such as homopolymers, GC content, local GC imbalance, repeated motifs, and secondary-structure-like patterns. Can we design an adaptive encoding algorithm that stores information while producing more DNA-friendly sequences than a naive fixed 2-bit mapping?

This is an educational simulator, not a production DNA synthesis tool, a medical diagnostic tool, or a real-world archival DNA storage system.

## 2. Core idea

A naive digital-to-DNA mapping can use four bases as four symbols:

```text
00 -> A
01 -> C
10 -> G
11 -> T
```

This is compact, because one DNA base can theoretically store 2 bits of information. However, this method ignores DNA-specific constraints. For example, it can produce long runs such as `AAAAAA`, skewed GC ratios, or repetitive motifs.

The improved algorithm should compare three encoding modes:

1. **Fixed 2-bit encoding**
   - Baseline method.
   - Highest simple information density.
   - Does not avoid biologically unfavorable patterns.

2. **Rotating ternary encoding**
   - Converts bytes into base-3 trits.
   - At each DNA position, avoids selecting the same base as the previous base.
   - Prevents homopolymers by construction.
   - Lower density than fixed 2-bit encoding.

3. **Adaptive optimized encoding**
   - Generates multiple reversible candidate sequences for each chunk.
   - Scores each candidate using DNA-friendly metrics.
   - Selects the candidate with the best score.
   - Stores all parameters required for decoding in metadata.

## 3. Biological and engineering constraints to consider

The adaptive algorithm should consider the following constraints.

### 3.1 Homopolymer runs

Long repeated bases such as `AAAAAA` or `GGGGGG` can be difficult for some synthesis or sequencing workflows.

Metric:

- `max_homopolymer_run(seq)`

Goal:

- Fixed 2-bit mode may produce long runs.
- Rotating ternary mode should usually have max run = 1.
- Adaptive mode should strongly penalize long runs.

### 3.2 Global GC content

GC content is the proportion of G and C bases in a sequence.

Metric:

- `gc_content(seq)`

Goal:

- Prefer roughly balanced GC content.
- A practical educational target is around 40-60%, ideally near 50%.

### 3.3 Local GC content

A sequence can have good global GC content while still containing local AT-rich or GC-rich regions. Therefore, local sliding-window GC should also be measured.

Metric:

- `local_gc_stats(seq, window=30)`

Goal:

- Keep local GC min/max closer to the target range.

### 3.4 Forbidden motifs

Some short motifs may be undesirable in synthetic DNA workflows. This project should include a small example list for educational purposes, such as common restriction enzyme recognition sequences.

Example motifs:

```python
FORBIDDEN_MOTIFS = [
    "GAATTC",  # EcoRI
    "AAGCTT",  # HindIII
    "GGATCC",  # BamHI
]
```

Metric:

- `count_forbidden_motifs(seq)`

Goal:

- Penalize candidates containing forbidden motifs.

### 3.5 Reverse-complement / hairpin-like patterns

Full DNA secondary structure prediction is out of scope. However, the app can implement a simplified educational approximation by detecting short reverse-complement palindromes.

Metric:

- `hairpin_like_score(seq)`

Goal:

- Penalize sequences with many short self-complementary patterns.

### 3.6 Repeated low-complexity motifs

Patterns such as `ATATATAT` or `CAGCAGCAG` can reduce sequence complexity and may complicate alignment or reconstruction.

Metric:

- `repeat_score(seq)`

Goal:

- Penalize repeated 2-mer or 3-mer patterns.

### 3.7 Chunking and addressing

Realistic DNA storage does not place an entire file into one very long strand. It divides data into shorter DNA oligo-like chunks. Since DNA strands are unordered, each chunk needs an index.

Design:

```text
[Preamble][Header][Payload][Checksum]
```

The header should include at least:

- version
- codec id
- chunk index
- total chunks
- original length or payload length
- seed
- mapping id
- start base
- compression flag
- CRC32 checksum

### 3.8 Error detection

This project does not need full error correction such as Reed-Solomon. However, it should implement simple error detection using CRC32.

Goal:

- If a DNA character is changed, the decoder should detect a checksum mismatch.
- The app should report checksum success/failure clearly.

## 4. Adaptive encoding concept

The adaptive encoder should be reversible. Therefore, it must not make random choices that cannot be reproduced during decoding.

Recommended approach:

1. Split input bytes into chunks.
2. For each chunk, generate candidate DNA sequences by trying:
   - multiple deterministic whitening seeds
   - multiple mapping permutations
   - multiple start bases
3. Analyze each candidate sequence.
4. Score each candidate.
5. Select the lowest-score candidate.
6. Store the selected seed, mapping id, start base, and other metadata in the header.
7. Decode using the stored metadata.

The adaptive part is candidate selection, not non-deterministic decoding.

## 5. Candidate generation parameters

Suggested parameters:

- `seed_trials`: 16 or 32 by default
- `mapping_id`: 0-5, representing all permutations of trit values [0, 1, 2]
- `start_base`: A, C, G, T
- `chunk_size`: 128 bytes by default

Example candidate count:

```text
32 seeds x 6 mappings x 4 start bases = 768 candidates per chunk
```

This is still small enough for a high-school-level Streamlit demo.

## 6. Scoring function

A simple weighted score is enough:

```text
score =
  global_gc_penalty
+ local_gc_penalty
+ homopolymer_penalty
+ forbidden_motif_penalty
+ hairpin_penalty
+ repeat_penalty
```

The exact weights can be adjusted. The important part is to show the trade-off between information density and DNA-friendliness.

Recommended first weights:

```python
weights = {
    "global_gc": 100.0,
    "local_gc": 80.0,
    "homopolymer": 50.0,
    "forbidden_motif": 30.0,
    "hairpin": 10.0,
    "repeat": 10.0,
}
```

## 7. Comparison metrics

The app should show comparison tables for all three modes.

Metrics:

- DNA length
- original bytes
- bits/base
- A/C/G/T counts
- GC content
- local GC minimum and maximum
- maximum homopolymer run
- forbidden motif count
- hairpin-like score
- repeat score
- decode success
- checksum status, if applicable

The report should compare compactness and DNA-friendliness.

## 8. Expected conclusion

Expected qualitative result:

- Fixed 2-bit encoding is usually shortest.
- Rotating ternary encoding avoids homopolymers but increases length.
- Adaptive optimized encoding can improve GC balance, reduce motifs, and reduce repeated patterns by testing multiple reversible candidates.

The project should emphasize that more biological constraints usually reduce raw density, but can improve reliability and reduce the need for redundancy or repair.

## 9. Limitations

The project should clearly state limitations:

- It does not synthesize real DNA.
- It does not perform real sequencing simulation.
- It does not model all biochemical constraints.
- Hairpin scoring is simplified.
- Forbidden motifs are examples only.
- CRC32 only detects errors; it does not correct them.
- This is an educational simulator for understanding DNA data storage principles.

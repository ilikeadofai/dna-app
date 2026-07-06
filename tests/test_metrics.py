import pytest

from dna_codec.constraints import COMPLEMENT, FORBIDDEN_MOTIFS
from dna_codec.metrics import (
    analyze_sequence,
    base_counts,
    bits_per_base,
    count_forbidden_motifs,
    gc_content,
    hairpin_like_score,
    local_gc_stats,
    max_homopolymer_run,
    normalize_dna,
    repeat_score,
    reverse_complement,
    validate_dna,
)


def test_constraints_defaults():
    assert FORBIDDEN_MOTIFS == ("GAATTC", "AAGCTT", "GGATCC")
    assert COMPLEMENT == {"A": "T", "C": "G", "G": "C", "T": "A"}


def test_normalize_dna_uppercases_and_strips_whitespace():
    assert normalize_dna(" acg t\n") == "ACGT"
    assert normalize_dna("") == ""


def test_validate_dna_rejects_invalid_characters():
    validate_dna("ACGTacgt\n")

    for invalid in ["ACNT", "ACGU", "ACG-1"]:
        with pytest.raises(ValueError):
            validate_dna(invalid)


def test_base_counts_all_keys_and_empty():
    assert base_counts("AACGT") == {"A": 2, "C": 1, "G": 1, "T": 1}
    assert base_counts("") == {"A": 0, "C": 0, "G": 0, "T": 0}


def test_gc_content_ratio():
    assert gc_content("AAGC") == pytest.approx(0.5)
    assert gc_content("GGCC") == pytest.approx(1.0)
    assert gc_content("AATT") == pytest.approx(0.0)
    assert gc_content("") == pytest.approx(0.0)


def test_max_homopolymer_run():
    assert max_homopolymer_run("") == 0
    assert max_homopolymer_run("A") == 1
    assert max_homopolymer_run("AAACCCCGTT") == 4


def test_local_gc_stats_sliding_window():
    stats = local_gc_stats("AAAACCCC", window=4)

    assert stats["values"] == pytest.approx([0.0, 0.25, 0.5, 0.75, 1.0])
    assert stats["min"] == pytest.approx(0.0)
    assert stats["max"] == pytest.approx(1.0)
    assert stats["mean"] == pytest.approx(0.5)
    assert stats["window"] == 4
    assert stats["effective_window"] == 4


def test_local_gc_stats_short_sequence_empty_and_bad_window():
    short_stats = local_gc_stats("ACG", window=30)
    assert short_stats["values"] == pytest.approx([2 / 3])
    assert short_stats["min"] == pytest.approx(2 / 3)
    assert short_stats["max"] == pytest.approx(2 / 3)
    assert short_stats["mean"] == pytest.approx(2 / 3)
    assert short_stats["window"] == 30
    assert short_stats["effective_window"] == 3

    empty_stats = local_gc_stats("", window=30)
    assert empty_stats["values"] == []
    assert empty_stats["min"] == pytest.approx(0.0)
    assert empty_stats["max"] == pytest.approx(0.0)
    assert empty_stats["mean"] == pytest.approx(0.0)

    for bad_window in [0, -1]:
        with pytest.raises(ValueError):
            local_gc_stats("ACGT", window=bad_window)


def test_count_forbidden_motifs_defaults():
    assert count_forbidden_motifs("GAATTCAAGCTTGGATCC") == 3


def test_count_forbidden_motifs_custom_overlap():
    assert count_forbidden_motifs("AAAAA", motifs=["AAA"]) == 3


def test_count_forbidden_motifs_empty_and_invalid_motif():
    assert count_forbidden_motifs("GAATTC", motifs=[]) == 0

    with pytest.raises(ValueError):
        count_forbidden_motifs("AAAA", motifs=["ANN"])

    with pytest.raises(ValueError):
        count_forbidden_motifs("AAAA", motifs=[""])


def test_reverse_complement():
    assert reverse_complement("") == ""
    assert reverse_complement("ACGA") == "TCGT"
    assert reverse_complement("aagctt") == "AAGCTT"

    with pytest.raises(ValueError):
        reverse_complement("ACGN")


def test_hairpin_like_score_basic_cases():
    assert hairpin_like_score("AAAATTTT", min_stem=4, max_stem=4) > 0
    assert hairpin_like_score("AAAACCCC", min_stem=4, max_stem=4) == 0
    assert hairpin_like_score("ACG", min_stem=4, max_stem=4) == 0

    with pytest.raises(ValueError):
        hairpin_like_score("AAAATTTT", min_stem=0, max_stem=4)

    with pytest.raises(ValueError):
        hairpin_like_score("AAAATTTT", min_stem=5, max_stem=4)


def test_repeat_score_basic_cases():
    assert repeat_score("ATATAT", k_values=(2,), min_repeats=3) > 0
    assert repeat_score("CAGCAGCAG", k_values=(3,), min_repeats=3) > 0
    assert repeat_score("ACGTACGA", k_values=(2, 3), min_repeats=3) == 0
    assert repeat_score("AC", k_values=(2,), min_repeats=3) == 0

    with pytest.raises(ValueError):
        repeat_score("ATATAT", k_values=(), min_repeats=3)

    with pytest.raises(ValueError):
        repeat_score("ATATAT", k_values=(0,), min_repeats=3)

    with pytest.raises(ValueError):
        repeat_score("ATATAT", k_values=(2,), min_repeats=1)


def test_bits_per_base():
    assert bits_per_base(16, 8) == pytest.approx(2.0)
    assert bits_per_base(0, 8) == pytest.approx(0.0)

    with pytest.raises(ValueError):
        bits_per_base(1, 0)

    with pytest.raises(ValueError):
        bits_per_base(-1, 8)


def test_analyze_sequence_integrates_metrics():
    result = analyze_sequence("ACGT", original_bits=8, window=2)

    assert result["length"] == 4
    assert result["base_counts"] == {"A": 1, "C": 1, "G": 1, "T": 1}
    assert result["gc_content"] == pytest.approx(0.5)
    assert result["gc_percent"] == pytest.approx(50.0)
    assert result["local_gc"]["window"] == 2
    assert result["max_homopolymer_run"] == 1
    assert result["forbidden_motif_count"] == 0
    assert result["bits_per_base"] == pytest.approx(2.0)
    assert result["original_bits"] == 8


def test_analyze_sequence_empty_is_stable():
    result = analyze_sequence("", original_bits=0, window=30)

    assert result["length"] == 0
    assert result["base_counts"] == {"A": 0, "C": 0, "G": 0, "T": 0}
    assert result["gc_content"] == pytest.approx(0.0)
    assert result["local_gc"]["values"] == []
    assert result["bits_per_base"] == pytest.approx(0.0)


def test_public_metrics_reject_invalid_sequence():
    invalid = "ACGN"

    with pytest.raises(ValueError):
        base_counts(invalid)
    with pytest.raises(ValueError):
        gc_content(invalid)
    with pytest.raises(ValueError):
        max_homopolymer_run(invalid)
    with pytest.raises(ValueError):
        local_gc_stats(invalid)
    with pytest.raises(ValueError):
        count_forbidden_motifs(invalid)
    with pytest.raises(ValueError):
        hairpin_like_score(invalid)
    with pytest.raises(ValueError):
        repeat_score(invalid)
    with pytest.raises(ValueError):
        analyze_sequence(invalid)

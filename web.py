from __future__ import annotations

import csv
import io
from typing import Any

import streamlit as st

try:
    import docx
except ImportError:  # pragma: no cover - depends on optional runtime package
    docx = None

try:
    import fitz  # PyMuPDF
except ImportError:  # pragma: no cover - depends on optional runtime package
    fitz = None

try:
    import pandas as pd
except ImportError:  # pragma: no cover - Streamlit can still show tables without charts
    pd = None

from dna_codec import (
    DEFAULT_CHUNK_SIZE,
    DEFAULT_SEED_TRIALS,
    AdaptiveCodecError,
    ChecksumMismatchError,
    adaptive_decode_chunks,
    adaptive_encode_chunks,
    analyze_sequence,
    bytes_to_fixed_dna,
    fixed_dna_to_bytes,
    normalize_dna,
    parse_adaptive_strand,
    rotating_decode,
    rotating_encode,
)

CODEC_FIXED = "Fixed 2-bit"
CODEC_ROTATING = "Rotating ternary"
CODEC_ADAPTIVE = "Adaptive optimized"
CODEC_OPTIONS = [CODEC_FIXED, CODEC_ROTATING, CODEC_ADAPTIVE]


st.set_page_config(
    page_title="DNA 데이터 저장 시뮬레이터", page_icon="🧬", layout="centered"
)


def extract_text(uploaded_file) -> str:
    """Extract text from txt/md/docx/pdf uploads."""
    extension = uploaded_file.name.rsplit(".", 1)[-1].lower()

    if extension in {"txt", "md"}:
        return uploaded_file.read().decode("utf-8")

    if extension == "docx":
        if docx is None:
            raise RuntimeError(
                "python-docx가 설치되어 있지 않아 DOCX 파일을 읽을 수 없습니다."
            )
        document = docx.Document(uploaded_file)
        return "\n".join(paragraph.text for paragraph in document.paragraphs)

    if extension == "pdf":
        if fitz is None:
            raise RuntimeError(
                "PyMuPDF가 설치되어 있지 않아 PDF 파일을 읽을 수 없습니다."
            )
        pdf = fitz.open(stream=uploaded_file.read(), filetype="pdf")
        page_texts: list[str] = [str(page.get_text()) for page in pdf]
        return "\n".join(page_texts)

    raise ValueError(f"지원하지 않는 파일 형식입니다: {extension}")


def read_text_input(uploaded_file, manual_text: str) -> tuple[str, str]:
    if uploaded_file is not None:
        text = extract_text(uploaded_file)
        return text, uploaded_file.name
    return manual_text, "직접 입력"


def read_dna_input(uploaded_file, manual_dna: str) -> tuple[str, str]:
    if uploaded_file is not None:
        dna = uploaded_file.read().decode("utf-8")
        return dna, uploaded_file.name
    return manual_dna, "직접 입력"


def split_adaptive_strands(raw_dna: str) -> list[str]:
    """Return newline-separated adaptive strands, normalized per strand."""
    lines = [line.strip() for line in raw_dna.splitlines() if line.strip()]
    if not lines and raw_dna.strip():
        lines = [raw_dna.strip()]

    if len(lines) <= 1:
        normalized = normalize_dna(raw_dna)
        return [normalized] if normalized else []

    return [normalize_dna(line) for line in lines]


def format_percent(value: float) -> str:
    return f"{value * 100:.1f}%"


def metric_rows(
    seq: str, original_bits: int | None = None, window: int = 30
) -> dict[str, Any]:
    normalized = normalize_dna(seq)
    bits = 0 if original_bits == 0 else original_bits
    return analyze_sequence(normalized, original_bits=bits, window=window)


def render_metrics(
    seq: str, original_bits: int | None = None, window: int = 30
) -> None:
    metrics = metric_rows(seq, original_bits=original_bits, window=window)
    counts = metrics["base_counts"]
    local_gc = metrics["local_gc"]

    cols = st.columns(4)
    cols[0].metric("DNA bases", metrics["length"])
    cols[1].metric("GC", f"{metrics['gc_percent']:.1f}%")
    cols[2].metric("Max run", metrics["max_homopolymer_run"])
    bits_per_base = metrics["bits_per_base"]
    cols[3].metric(
        "Bits/base", "-" if bits_per_base is None else f"{bits_per_base:.2f}"
    )

    st.write(
        {
            "A": counts["A"],
            "C": counts["C"],
            "G": counts["G"],
            "T": counts["T"],
            "Local GC min": f"{local_gc['min'] * 100:.1f}%",
            "Local GC mean": f"{local_gc['mean'] * 100:.1f}%",
            "Local GC max": f"{local_gc['max'] * 100:.1f}%",
            "Forbidden motifs": metrics["forbidden_motif_count"],
            "Hairpin-like score": metrics["hairpin_like_score"],
            "Repeat score": metrics["repeat_score"],
        }
    )

    if pd is not None and metrics["length"] > 0:
        chart_data = pd.DataFrame(
            {
                "Base": ["A", "C", "G", "T"],
                "Count": [counts["A"], counts["C"], counts["G"], counts["T"]],
            }
        )
        st.bar_chart(chart_data, x="Base", y="Count", use_container_width=True)


def encode_selected_codec(
    data: bytes,
    codec: str,
    *,
    rotating_start_base: str,
    rotating_mapping_id: int,
    adaptive_chunk_size: int,
    adaptive_seed_trials: int,
    adaptive_compressed: bool,
    window: int,
) -> tuple[str, dict[str, Any], str]:
    if codec == CODEC_FIXED:
        dna = bytes_to_fixed_dna(data)
        return (
            dna,
            {
                "codec": codec,
                "original_bytes": len(data),
                "metadata": "No external metadata required.",
            },
            dna,
        )

    if codec == CODEC_ROTATING:
        dna = rotating_encode(
            data, start_base=rotating_start_base, mapping_id=rotating_mapping_id
        )
        return (
            dna,
            {
                "codec": codec,
                "original_bytes": len(data),
                "start_base": rotating_start_base,
                "mapping_id": rotating_mapping_id,
                "metadata": "Decoding requires original byte length, start base, and mapping id.",
            },
            dna,
        )

    encoded_chunks = adaptive_encode_chunks(
        data,
        chunk_size=adaptive_chunk_size,
        seed_trials=adaptive_seed_trials,
        compressed=adaptive_compressed,
        window=window,
    )
    dna = "\n".join(item.dna for item in encoded_chunks)
    payload_dna_length = sum(len(item.payload_dna) for item in encoded_chunks)
    total_dna_length = sum(len(item.dna) for item in encoded_chunks)
    payload_dna = "".join(item.payload_dna for item in encoded_chunks)
    return (
        dna,
        {
            "codec": codec,
            "original_bytes": len(data),
            "chunks": len(encoded_chunks),
            "compressed": adaptive_compressed,
            "payload_dna_length": payload_dna_length,
            "header_dna_length": total_dna_length - payload_dna_length,
            "total_dna_length": total_dna_length,
            "first_chunk_seed": encoded_chunks[0].header.seed,
            "first_chunk_mapping_id": encoded_chunks[0].header.mapping_id,
            "first_chunk_start_base": encoded_chunks[0].header.start_base,
            "metadata": "Self-describing strands with compact binary header and CRC32 checksum.",
        },
        payload_dna,
    )


def decode_selected_codec(
    raw_dna: str,
    codec: str,
    *,
    rotating_length: int,
    rotating_start_base: str,
    rotating_mapping_id: int,
) -> tuple[bytes, dict[str, Any]]:
    if codec == CODEC_FIXED:
        seq = normalize_dna(raw_dna)
        return fixed_dna_to_bytes(seq), {"codec": codec, "checksum": "Not used"}

    if codec == CODEC_ROTATING:
        seq = normalize_dna(raw_dna)
        return rotating_decode(
            seq,
            length=rotating_length,
            start_base=rotating_start_base,
            mapping_id=rotating_mapping_id,
        ), {
            "codec": codec,
            "byte_length": rotating_length,
            "start_base": rotating_start_base,
            "mapping_id": rotating_mapping_id,
            "checksum": "Not used",
        }

    strands = split_adaptive_strands(raw_dna)
    if not strands:
        raise ValueError("DNA 입력이 비어 있습니다.")

    headers = [parse_adaptive_strand(strand)[0] for strand in strands]
    data = adaptive_decode_chunks(strands)

    return data, {
        "codec": codec,
        "checksum": "CRC32 verified",
        "chunks": len(strands),
        "total_chunks": headers[0].total_chunks,
        "compressed": any(header.compressed for header in headers),
    }


def comparison_rows(
    data: bytes,
    *,
    adaptive_chunk_size: int,
    adaptive_seed_trials: int,
    adaptive_compressed: bool,
    window: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    original_bits = len(data) * 8

    fixed_payload_dna = bytes_to_fixed_dna(data)
    rows.append(
        build_comparison_row(
            CODEC_FIXED,
            fixed_payload_dna,
            header_dna_length=0,
            original_bits=original_bits,
            window=window,
        )
    )

    rotating_payload_dna = rotating_encode(data, start_base="A", mapping_id=0)
    rows.append(
        build_comparison_row(
            CODEC_ROTATING,
            rotating_payload_dna,
            header_dna_length=0,
            original_bits=original_bits,
            window=window,
        )
    )

    adaptive_chunks = adaptive_encode_chunks(
        data,
        chunk_size=adaptive_chunk_size,
        seed_trials=adaptive_seed_trials,
        compressed=adaptive_compressed,
        window=window,
    )
    adaptive_payload_dna = "".join(item.payload_dna for item in adaptive_chunks)
    adaptive_total_dna_length = sum(len(item.dna) for item in adaptive_chunks)
    adaptive_header_dna_length = adaptive_total_dna_length - len(adaptive_payload_dna)
    rows.append(
        build_comparison_row(
            CODEC_ADAPTIVE,
            adaptive_payload_dna,
            header_dna_length=adaptive_header_dna_length,
            original_bits=original_bits,
            window=window,
        )
    )

    return rows


def build_comparison_row(
    codec: str,
    payload_dna: str,
    *,
    header_dna_length: int,
    original_bits: int,
    window: int,
) -> dict[str, Any]:
    normalized_payload = normalize_dna(payload_dna)
    metrics = analyze_sequence(normalized_payload, window=window)
    payload_dna_length = metrics["length"]
    total_dna_length = payload_dna_length + header_dna_length
    return {
        "codec": codec,
        "payload_dna_length": payload_dna_length,
        "header_dna_length": header_dna_length,
        "total_dna_length": total_dna_length,
        "bits_per_base_payload": round(original_bits / payload_dna_length, 3)
        if payload_dna_length
        else 0.0,
        "bits_per_base_total": round(original_bits / total_dna_length, 3)
        if total_dna_length
        else 0.0,
        "gc_percent": round(metrics["gc_percent"], 1),
        "local_gc_min": round(metrics["local_gc"]["min"] * 100, 1),
        "local_gc_max": round(metrics["local_gc"]["max"] * 100, 1),
        "max_homopolymer": metrics["max_homopolymer_run"],
        "forbidden_motif_count": metrics["forbidden_motif_count"],
        "hairpin_like_score": metrics["hairpin_like_score"],
        "repeat_score": metrics["repeat_score"],
    }


def rows_to_csv(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return ""
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def render_app_header() -> None:
    st.title("DNA 데이터 저장 시뮬레이터")
    st.caption(
        "UTF-8 데이터를 DNA 염기서열로 변환하고, GC 함량·homopolymer·반복 패턴을 비교하는 교육용 도구입니다."
    )
    st.info(
        "이 앱은 교육용 시뮬레이터이며 실제 DNA 합성, 시퀀싱, 진단 도구가 아닙니다."
    )


def render_encode_tab() -> None:
    st.subheader("Encode")
    uploaded_file = st.file_uploader(
        "파일 업로드", type=["txt", "md", "docx", "pdf"], key="encode_file"
    )
    manual_text = st.text_area("텍스트 입력", height=160, key="encode_text")

    codec = st.selectbox("Encoding mode", CODEC_OPTIONS, key="encode_codec")
    window = st.number_input(
        "Local GC window", min_value=1, max_value=300, value=30, step=1
    )

    rotating_start_base = "A"
    rotating_mapping_id = 0
    adaptive_chunk_size = DEFAULT_CHUNK_SIZE
    adaptive_seed_trials = DEFAULT_SEED_TRIALS
    adaptive_compressed = False

    if codec == CODEC_ROTATING:
        cols = st.columns(2)
        rotating_start_base = cols[0].selectbox(
            "Start base", ["A", "C", "G", "T"], key="encode_rot_start"
        )
        rotating_mapping_id = cols[1].number_input(
            "Mapping ID", min_value=0, max_value=5, value=0, step=1
        )

    if codec == CODEC_ADAPTIVE:
        cols = st.columns(3)
        adaptive_chunk_size = cols[0].number_input(
            "Chunk size", min_value=1, max_value=4096, value=DEFAULT_CHUNK_SIZE, step=16
        )
        adaptive_seed_trials = cols[1].number_input(
            "Seed trials", min_value=1, max_value=64, value=DEFAULT_SEED_TRIALS, step=1
        )
        adaptive_compressed = cols[2].checkbox("Compression", value=False)

    try:
        text, source = read_text_input(uploaded_file, manual_text)
    except Exception as exc:
        st.error(f"파일을 읽을 수 없습니다: {exc}")
        return

    if text:
        st.caption(f"Input source: {source} · {len(text.encode('utf-8'))} bytes")
        with st.expander("Input preview", expanded=False):
            st.text_area("Preview", text, height=120, disabled=True)
    else:
        st.info("텍스트를 입력하거나 파일을 업로드하세요.")

    encode_clicked = st.button("Encode", type="primary", disabled=not bool(text))
    if encode_clicked:
        try:
            data = text.encode("utf-8")
            with st.spinner("Encoding..."):
                dna, metadata, payload_dna = encode_selected_codec(
                    data,
                    codec,
                    rotating_start_base=rotating_start_base,
                    rotating_mapping_id=int(rotating_mapping_id),
                    adaptive_chunk_size=int(adaptive_chunk_size),
                    adaptive_seed_trials=int(adaptive_seed_trials),
                    adaptive_compressed=adaptive_compressed,
                    window=int(window),
                )

            st.success("Encoding complete")
            st.text_area("DNA output", dna, height=240)
            st.download_button("Download DNA", dna, file_name="encoded_dna.txt")

            st.write("Metadata")
            st.json(metadata)
            st.write("Payload metrics")
            render_metrics(payload_dna, original_bits=len(data) * 8, window=int(window))
            if metadata.get("header_dna_length", 0):
                with st.expander("Packaged strand metrics", expanded=False):
                    render_metrics(dna, original_bits=len(data) * 8, window=int(window))
        except Exception as exc:
            st.error(f"Encoding failed: {exc}")

    st.divider()
    st.write("Comparison")
    st.caption(
        "Payload metrics are compared fairly across codecs. Header/checksum overhead is shown separately."
    )
    st.info(
        "Fixed 2-bit is expected to have the best raw density because one DNA base can encode at most 2 bits in a simple 4-symbol alphabet. Adaptive optimized mode aims to improve DNA-friendliness, not beat the theoretical raw density of Fixed 2-bit."
    )
    if st.button("Compare all codecs", disabled=not bool(text)):
        try:
            data = text.encode("utf-8")
            with st.spinner("Building comparison..."):
                rows = comparison_rows(
                    data,
                    adaptive_chunk_size=int(adaptive_chunk_size),
                    adaptive_seed_trials=int(adaptive_seed_trials),
                    adaptive_compressed=adaptive_compressed,
                    window=int(window),
                )
            st.dataframe(rows, use_container_width=True)
            st.download_button(
                "Download comparison CSV",
                rows_to_csv(rows),
                file_name="dna_comparison.csv",
            )
        except Exception as exc:
            st.error(f"Comparison failed: {exc}")


def render_decode_tab() -> None:
    st.subheader("Decode")
    uploaded_file = st.file_uploader("DNA 파일 업로드", type=["txt"], key="decode_file")
    manual_dna = st.text_area("DNA 입력", height=180, key="decode_dna")
    codec = st.selectbox("Decoding mode", CODEC_OPTIONS, key="decode_codec")

    rotating_length = 0
    rotating_start_base = "A"
    rotating_mapping_id = 0
    if codec == CODEC_ROTATING:
        cols = st.columns(3)
        rotating_length = cols[0].number_input(
            "Original byte length", min_value=0, value=0, step=1
        )
        rotating_start_base = cols[1].selectbox(
            "Start base", ["A", "C", "G", "T"], key="decode_rot_start"
        )
        rotating_mapping_id = cols[2].number_input(
            "Mapping ID",
            min_value=0,
            max_value=5,
            value=0,
            step=1,
            key="decode_rot_mapping",
        )

    try:
        raw_dna, source = read_dna_input(uploaded_file, manual_dna)
    except Exception as exc:
        st.error(f"DNA 파일을 읽을 수 없습니다: {exc}")
        return

    if raw_dna.strip():
        st.caption(f"Input source: {source}")
    else:
        st.info("DNA 시퀀스를 입력하거나 txt 파일을 업로드하세요.")

    if st.button("Decode", type="primary", disabled=not bool(raw_dna.strip())):
        try:
            data, metadata = decode_selected_codec(
                raw_dna,
                codec,
                rotating_length=int(rotating_length),
                rotating_start_base=rotating_start_base,
                rotating_mapping_id=int(rotating_mapping_id),
            )
            text = data.decode("utf-8")
            st.success("Decode complete")
            st.write("Status")
            st.json(metadata)
            st.text_area("Decoded text", text, height=220)
            st.download_button(
                "Download decoded text", text, file_name="decoded_text.txt"
            )
        except ChecksumMismatchError as exc:
            st.error(f"Checksum failed: {exc}")
        except UnicodeDecodeError:
            st.error("DNA는 bytes로 복원되었지만 UTF-8 텍스트로 해석할 수 없습니다.")
        except AdaptiveCodecError as exc:
            st.error(f"Adaptive decode failed: {exc}")
        except Exception as exc:
            st.error(f"Decode failed: {exc}")


def render_analyze_tab() -> None:
    st.subheader("Analyze DNA")
    uploaded_file = st.file_uploader(
        "DNA 파일 업로드", type=["txt"], key="analyze_file"
    )
    manual_dna = st.text_area("DNA 입력", height=180, key="analyze_dna")
    window = st.number_input(
        "Local GC window",
        min_value=1,
        max_value=300,
        value=30,
        step=1,
        key="analyze_window",
    )

    try:
        raw_dna, source = read_dna_input(uploaded_file, manual_dna)
    except Exception as exc:
        st.error(f"DNA 파일을 읽을 수 없습니다: {exc}")
        return

    if raw_dna.strip():
        st.caption(f"Input source: {source}")
    else:
        st.info("분석할 DNA 시퀀스를 입력하거나 txt 파일을 업로드하세요.")

    if st.button("Analyze", type="primary", disabled=not bool(raw_dna.strip())):
        try:
            normalized = normalize_dna(raw_dna)
            render_metrics(normalized, window=int(window))
            st.caption(
                "Hairpin-like score는 교육용 단순 지표이며 실제 2차 구조 예측이 아닙니다."
            )
        except Exception as exc:
            st.error(f"Analysis failed: {exc}")


def main() -> None:
    render_app_header()
    encode_tab, decode_tab, analyze_tab = st.tabs(["Encode", "Decode", "Analyze"])

    with encode_tab:
        render_encode_tab()
    with decode_tab:
        render_decode_tab()
    with analyze_tab:
        render_analyze_tab()


if __name__ == "__main__":
    main()

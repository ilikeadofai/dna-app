import streamlit as st
import requests
import json
import os
import fitz  # PyMuPDF
import docx

API_URL = "https://epvdcbyvjhcta77pqfye2shufu0fqgwu.lambda-url.eu-west-1.on.aws/"
LOOKUP_FILE = "dna_lookup.json"

# ▒▒▒ 로컬 룩업 테이블 ▒▒▒
def load_lookup_table():
    lookup = {}
    # 사용자 정의 lookup
    if os.path.exists(LOOKUP_FILE):
        with open(LOOKUP_FILE, "r", encoding="utf-8") as f:
            lookup.update(json.load(f))
    # prebuilt lookup 자동 병합
    if os.path.exists("prebuilt_lookup.json"):
        with open("prebuilt_lookup.json", "r", encoding="utf-8") as f:
            lookup.update(json.load(f))
    return lookup

def save_lookup_table(lookup):
    with open(LOOKUP_FILE, "w", encoding="utf-8") as f:
        json.dump(lookup, f, ensure_ascii=False, indent=2)

# ▒▒▒ 누락 문자 API 요청 및 테이블 업데이트 ▒▒▒
def update_lookup_table(text, lookup):
    updated = False
    missing = [ch for ch in set(text) if ch not in lookup]

    for ch in missing:
        try:
            res = requests.post(API_URL, json={"text": ch}, timeout=5)
            res.raise_for_status()
            dna = res.text.strip()
            lookup[ch] = dna
            updated = True
        except Exception as e:
            st.warning(f"⚠️ 문자 '{ch}' 인코딩 실패: {e}")

    if updated:
        save_lookup_table(lookup)

    return lookup

# ▒▒▒ 인코딩 / 복호화 함수 ▒▒▒
def encode_to_dna(text, lookup):
    return ''.join(lookup.get(ch, '') for ch in text)

def build_reverse_lookup(lookup):
    return {v: k for k, v in lookup.items()}

def decode_from_dna(dna_seq, reverse_lookup):
    decoded = []
    max_len = max(len(k) for k in reverse_lookup)  # 가장 긴 시퀀스 기준
    i = 0
    while i < len(dna_seq):
        match = False
        for l in range(max_len, 7, -1):  # 24 → 16 → 8 → ...
            segment = dna_seq[i:i+l]
            if segment in reverse_lookup:
                decoded.append(reverse_lookup[segment])
                i += l
                match = True
                break
        if not match:
            decoded.append('?')
            i += 8  # fallback skip
    return ''.join(decoded)


# ▒▒▒ 파일 유형별 텍스트 추출 ▒▒▒
def extract_text(file, file_type):
    if file_type == "txt" or file_type == "md":
        return file.read().decode("utf-8")
    elif file_type == "docx":
        doc = docx.Document(file)
        return "\n".join([para.text for para in doc.paragraphs])
    elif file_type == "pdf":
        pdf = fitz.open(stream=file.read(), filetype="pdf")
        text = ""
        for page in pdf:
            text += page.get_text()
        return text
    else:
        return ""

# ▒▒▒ UI 시작 ▒▒▒
st.set_page_config(page_title="DNA 암호화기", page_icon="🧬")
st.title("🧬 문서 기반 DNA 암호화 / 복호화기")

lookup_table = load_lookup_table()
mode = st.radio("모드를 선택하세요", ["암호화", "복호화"])

# ▒▒▒ 암호화 모드 ▒▒▒
if mode == "암호화":
    st.subheader("🔼 파일 업로드 또는 직접 입력")
    file = st.file_uploader("지원 파일: .txt, .md, .docx, .pdf", type=["txt", "md", "docx", "pdf"])
    manual_input = st.text_area("또는 텍스트 직접 입력", "")

    text = ""
    if file:
        ext = file.name.split('.')[-1].lower()
        text = extract_text(file, ext)
        st.success(f"✅ 텍스트 추출 완료 ({ext})")
    elif manual_input.strip():
        text = manual_input.strip()

    if text:
        st.text_area("📄 원문 미리보기", text, height=200)
        if st.button("DNA로 암호화"):
            with st.spinner("DNA로 변환 중..."):
                lookup_table = update_lookup_table(text, lookup_table)
                dna_result = encode_to_dna(text, lookup_table)

            st.success("✅ 암호화 완료!")
            st.text_area("🧬 DNA 시퀀스", dna_result, height=200)
            st.download_button("DNA 파일 다운로드", dna_result, file_name="encrypted_dna.txt")
    else:
        st.info("파일 업로드 또는 직접 입력을 해주세요.")

# ▒▒▒ 복호화 모드 ▒▒▒
else:
    st.subheader("🔼 DNA 파일 업로드 또는 DNA 시퀀스 입력")
    dna_file = st.file_uploader("DNA 파일 (.txt)", type=["txt"])
    dna_input = st.text_area("또는 DNA 시퀀스를 직접 입력", "")

    dna_text = ""
    if dna_file:
        dna_text = dna_file.read().decode("utf-8").replace("\n", "")
    elif dna_input.strip():
        dna_text = dna_input.strip()

    if dna_text:
        st.text_area("🧬 DNA 시퀀스 미리보기", dna_text, height=200)
        if st.button("텍스트로 복호화"):
            reverse_lookup = build_reverse_lookup(lookup_table)  # <- 여기서 생성
            decoded_text = decode_from_dna(dna_text, reverse_lookup)

            st.success("✅ 복호화 완료!")
            st.text_area("📄 복원된 텍스트", decoded_text, height=200)
            st.download_button("텍스트 저장", decoded_text, file_name="decrypted_text.txt")
    else:
        st.info("DNA 시퀀스를 입력하거나 파일을 업로드하세요.")

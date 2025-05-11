import requests
import json
import time

API_URL = "https://epvdcbyvjhcta77pqfye2shufu0fqgwu.lambda-url.eu-west-1.on.aws/"

# 문자셋: 영문, 숫자, 기호, 이모지, 한글 전체
ascii_chars = list("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
symbols = list(" !@#$%^&*()-_=+[]{}|;:',.<>?/`~\\\"")
emoji = ["😀", "🧬", "✨", "🔥", "❤️"]
hangul = [chr(code) for code in range(44032, 55204)]  # '가' ~ '힣'

charset = ascii_chars + symbols + emoji + hangul

lookup = {}
failed = []

for i, ch in enumerate(charset, 1):
    try:
        res = requests.post(API_URL, json={"text": ch}, timeout=10)
        res.raise_for_status()
        dna = res.text.strip()
        lookup[ch] = dna
        print(f"[{i}/{len(charset)}] ✅ '{ch}' → {dna}")
    except Exception as e:
        print(f"[{i}/{len(charset)}] ❌ '{ch}' 실패: {e}")
        failed.append(ch)

# 저장
with open("prebuilt_lookup.json", "w", encoding="utf-8") as f:
    json.dump(lookup, f, ensure_ascii=False, indent=2)

with open("failed_chars.txt", "w", encoding="utf-8") as f:
    f.write("".join(failed))

print(f"\n✅ 완료: 총 {len(lookup)}개 성공, {len(failed)}개 실패")

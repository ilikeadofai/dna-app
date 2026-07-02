# Biology Report Notes

## 1. Suggested report title

Recommended title:

> DNA 염기서열의 특성을 고려한 적응형 정보 저장 알고리즘 설계

Alternative titles:

- Homopolymer 회피와 GC 함량 균형을 고려한 DNA 데이터 인코딩 알고리즘
- DNA의 분자적 특성을 활용한 정보 저장 방식 탐구
- DNA를 쓰고 읽는 교육용 웹앱 설계: 정보 저장과 변이 해석의 관점에서

## 2. Research motivation

생명과학 수업에서 DNA, 유전자, 염색체의 관계를 배우며 DNA가 생명 정보를 저장하는 물질이라는 점에 관심을 가지게 되었다. 컴퓨터공학에서도 정보의 저장과 해석이 핵심이므로, DNA 염기서열을 디지털 정보 저장 매체로 볼 수 있다는 점이 흥미로웠다.

1학년 때 제작했던 DNA 변환 웹앱은 문자열을 DNA 염기서열 형태로 바꾸는 기능을 구현했지만, 문자별 lookup table과 외부 API에 의존하는 한계가 있었다. 이번 탐구에서는 이를 개선하여 유니코드 문자를 UTF-8 바이트로 변환하고, DNA의 A, T, G, C 염기 배열로 직접 인코딩하는 방식을 설계하였다.

나아가 단순히 0과 1을 A, T, G, C에 대응시키는 방식은 DNA의 실제 분자적 특성을 충분히 고려하지 못한다고 보았다. 따라서 homopolymer, GC 함량, local GC 불균형, 반복 motif, 자기상보적 구조 가능성 등을 고려한 적응형 DNA 정보 저장 알고리즘을 설계하고자 하였다.

## 3. Biology concepts connected to the project

### 3.1 DNA and nucleotides

DNA는 뉴클레오타이드가 연결된 고분자이며, 각 뉴클레오타이드는 염기, 당, 인산으로 구성된다. DNA의 염기는 A, T, G, C 네 종류가 있으며, 이 염기 배열이 유전 정보를 저장한다.

이 탐구에서는 이 네 종류의 염기를 디지털 정보를 표현하는 기호로 보고, 컴퓨터 데이터와 DNA 염기서열 사이의 변환 규칙을 설계하였다.

### 3.2 Genes and chromosomes

염색체는 DNA와 단백질로 이루어진 구조이며, 유전자는 DNA 중 특정 기능을 담당하는 염기서열 구간이다. 즉 DNA 염기서열의 순서가 생명 정보의 핵심이다.

이 프로젝트는 DNA 염기서열이 정보를 담는다는 생명과학 개념을 컴퓨터공학의 정보 저장 개념과 연결한 사례이다.

### 3.3 GC content

GC 함량은 DNA 서열 중 G와 C가 차지하는 비율이다. GC 함량은 DNA의 안정성, 녹는 온도, 증폭 및 합성 조건과 관련될 수 있다. 따라서 디지털 데이터를 DNA로 저장할 때도 GC 함량이 지나치게 높거나 낮지 않도록 고려할 필요가 있다.

### 3.4 Homopolymer

Homopolymer는 같은 염기가 연속해서 반복되는 구간이다.

Examples:

```text
AAAAAA
CCCCCC
GGGGGG
TTTTTT
```

이러한 구간은 합성 또는 염기서열 분석 과정에서 오류 가능성을 높일 수 있으므로, 본 탐구에서는 이를 줄이는 인코딩 방식을 설계하였다.

### 3.5 Mutation and sequence comparison

DNA 서열의 변화는 SNP, 삽입, 결실 등으로 나타날 수 있다. 이 프로젝트의 핵심은 데이터 저장이지만, 이후 확장으로 기준 서열과 변이 서열을 비교하여 변이 유형을 교육용으로 분석하는 기능을 추가할 수 있다.

## 4. Computer science concepts connected to the project

### 4.1 Encoding and decoding

문자열은 컴퓨터 내부에서 UTF-8 바이트로 표현된다. 이 바이트를 비트 또는 3진수 형태로 변환한 뒤 DNA 염기서열로 대응시킬 수 있다.

Encoding:

```text
text -> UTF-8 bytes -> DNA sequence
```

Decoding:

```text
DNA sequence -> bytes -> text
```

### 4.2 Fixed mapping

가장 단순한 방식은 2비트를 하나의 염기에 대응시키는 것이다.

```text
00 -> A
01 -> C
10 -> G
11 -> T
```

이 방식은 짧고 효율적이지만 DNA의 생물학적 제약을 고려하지 않는다.

### 4.3 Rotating ternary encoding

데이터를 3진수로 바꾸고, 직전 염기와 같은 염기를 제외한 세 염기 중 하나를 선택하면 같은 염기의 반복을 피할 수 있다. 이 방식은 저장 길이는 증가하지만 homopolymer를 원천적으로 줄일 수 있다.

### 4.4 Adaptive optimization

적응형 방식은 여러 후보 서열을 만든 뒤, DNA 친화적인 지표를 기준으로 가장 좋은 후보를 선택한다.

고려 지표:

- GC 함량
- local GC 함량
- 최장 homopolymer 길이
- 금지 motif 개수
- hairpin-like score
- 반복 pattern score

## 5. Suggested report structure

### 1. 탐구 동기

- DNA가 유전 정보를 저장한다는 생명과학 개념을 배움.
- 컴퓨터공학의 데이터 저장과 DNA 정보 저장 사이의 공통점을 느낌.
- 기존 DNA 변환 웹앱을 개선하여 실제 DNA의 특성을 고려한 정보 저장 알고리즘을 설계하고자 함.

### 2. 이론적 배경

- DNA, 뉴클레오타이드, 염기서열
- 염색체와 유전자
- GC 함량
- Homopolymer
- DNA 저장 기술의 기본 아이디어

### 3. 알고리즘 설계

비교 대상:

1. Fixed 2-bit encoding
2. Rotating ternary encoding
3. Adaptive optimized encoding

Adaptive algorithm:

- 입력 데이터를 UTF-8 바이트로 변환
- 선택적으로 압축
- chunk로 분할
- 여러 seed, mapping, start base 후보 생성
- 각 후보의 DNA metric 계산
- 점수가 가장 좋은 후보 선택
- header와 checksum을 포함하여 복호화 가능성 보장

### 4. 구현 내용

- Streamlit 웹앱
- 로컬 인코딩/복호화
- DNA sequence metrics
- 비교 실험 표
- CSV 다운로드

### 5. 실험 결과

비교 표 예시:

| 입력 데이터 | 방식 | DNA 길이 | bits/base | GC% | local GC 범위 | max run | motif count | 복호화 |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 한글 문장 | Fixed 2-bit | | | | | | | |
| 한글 문장 | Rotating ternary | | | | | | | |
| 한글 문장 | Adaptive | | | | | | | |
| 코드 일부 | Fixed 2-bit | | | | | | | |
| 코드 일부 | Adaptive | | | | | | | |

### 6. 분석

예상 분석 방향:

- Fixed 2-bit 방식은 가장 짧지만 DNA-friendly metric이 나쁠 수 있다.
- Rotating ternary 방식은 길이가 늘어나지만 homopolymer를 줄일 수 있다.
- Adaptive 방식은 여러 후보 중 가장 좋은 서열을 선택하여 GC 함량과 반복 패턴 측면에서 개선될 수 있다.
- 저장 효율과 안정성 사이에는 trade-off가 있다.

### 7. 한계

- 실제 DNA 합성이나 시퀀싱을 수행한 것은 아니다.
- 2차 구조 예측은 단순화된 점수이다.
- CRC32는 오류 검출만 가능하고 오류 정정은 하지 못한다.
- 실제 DNA 저장에는 더 정교한 오류 정정, primer 설계, 합성 비용, sequencing error model이 필요하다.

### 8. 결론 및 추후 탐구

- DNA를 단순한 4진 저장 매체로만 보는 것은 부족하다.
- DNA의 분자적 특성을 고려해야 실제 저장에 적합한 서열을 만들 수 있다.
- 추후에는 Reed-Solomon 같은 오류 정정, 실제 sequencing error simulation, 변이 분석 교육용 workflow를 추가할 수 있다.

## 6. Report-ready paragraph

아래 문단은 보고서 도입부 또는 결론에 사용할 수 있다.

> DNA는 A, T, G, C 네 종류의 염기 배열을 통해 생명 정보를 저장한다. 이 점에서 DNA는 컴퓨터의 디지털 데이터처럼 정보를 표현할 수 있는 매체로 볼 수 있다. 그러나 단순히 0과 1을 염기서열로 치환하는 방식은 DNA의 실제 분자적 특성을 충분히 고려하지 못한다. 예를 들어 같은 염기가 길게 반복되는 homopolymer, 지나치게 높거나 낮은 GC 함량, 특정 반복 motif, 자기상보적 구조 가능성은 DNA 합성 및 염기서열 분석 과정에서 오류 가능성을 높일 수 있다. 따라서 본 탐구에서는 기존의 고정 2비트 매핑 방식과 비교하여, 여러 후보 DNA 서열을 생성하고 GC 함량, local GC 균형, homopolymer 길이, 금지 motif 개수 등을 기준으로 가장 적합한 서열을 선택하는 적응형 DNA 정보 저장 알고리즘을 설계하였다. 이를 통해 생명과학의 DNA 개념이 컴퓨터공학의 데이터 저장 및 알고리즘 설계와 연결될 수 있음을 확인하였다.

## 7. Important wording

Use these expressions in the report:

- DNA를 생명 정보 저장 매체이자 디지털 정보 표현 매체로 바라보았다.
- 단순한 2비트 치환이 아니라 DNA의 분자적 특성을 고려한 인코딩 방식을 설계하였다.
- 저장 효율과 합성/해독 안정성 사이의 trade-off를 비교하였다.
- 본 웹앱은 실제 진단 또는 생산용 DNA 저장 도구가 아니라 교육용 시뮬레이터이다.
- AI나 알고리즘의 결과는 실제 생물학적 검증을 대체할 수 없다.

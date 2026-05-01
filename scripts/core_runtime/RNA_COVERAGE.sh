#!/bin/bash
# Usage:
#   ./RNA_COVERAGE.sh SAMPLE_NAME FILTERED_BAM STAR_INDEX_DIR CHROM_SIZES OUTDIR THREADS

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 SAMPLE_NAME FILTERED_BAM STAR_INDEX_DIR CHROM_SIZES OUTDIR THREADS" >&2
    exit 1
fi

sample_name="${1}"
INBAM="${2}"
path_refDB="${3}"
path_refCHROMSIZES="${4}"
outdir="${5}"
threads="${6}"
STAR_BIN="${STAR_BIN:-STAR}"
BEDGRAPH_TO_BIGWIG_BIN="${BEDGRAPH_TO_BIGWIG_BIN:-bedGraphToBigWig}"

if [[ ! -d "${path_refDB}" ]]; then
    echo "ERROR: STAR index directory missing: ${path_refDB}" >&2
    exit 1
fi

if [[ ! -s "${path_refCHROMSIZES}" ]]; then
    echo "ERROR: chromosome sizes file missing or empty: ${path_refCHROMSIZES}" >&2
    exit 1
fi

if [[ ! -s "${INBAM}" ]]; then
    echo "ERROR: Input filtered BAM missing or empty: ${INBAM}" >&2
    exit 1
fi

echo "Using STAR_BIN=${STAR_BIN}"
echo "Using BEDGRAPH_TO_BIGWIG_BIN=${BEDGRAPH_TO_BIGWIG_BIN}"
echo "Using STAR index directory=${path_refDB}"
echo "Using chromosome sizes=${path_refCHROMSIZES}"

"${STAR_BIN}" \
  --runMode inputAlignmentsFromBAM \
  --runThreadN "${threads}" \
  --genomeDir "${path_refDB}" \
  --inputBAMfile "${INBAM}" \
  --outWigType bedGraph \
  --outWigStrand Stranded \
  --outWigNorm RPM \
  --outWigReferencesPrefix chr \
  --outFileNamePrefix "${outdir}/${sample_name}.stranded_"

"${STAR_BIN}" \
  --runMode inputAlignmentsFromBAM \
  --runThreadN "${threads}" \
  --genomeDir "${path_refDB}" \
  --inputBAMfile "${INBAM}" \
  --outWigType bedGraph \
  --outWigStrand Unstranded \
  --outWigNorm RPM \
  --outWigReferencesPrefix chr \
  --outFileNamePrefix "${outdir}/${sample_name}.unstranded_"

for f in "${outdir}/${sample_name}.stranded_Signal.Unique.str"*.bg; do
  sort -k1,1 -k2,2n "$f" > "${f%.bedGraph}.sorted.bg"
  "${BEDGRAPH_TO_BIGWIG_BIN}" "${f%.bedGraph}.sorted.bg" "${path_refCHROMSIZES}" "${f%.bg}.bw"
done

for f in "${outdir}/${sample_name}.unstranded_Signal.Unique.str"*.bg; do
  sort -k1,1 -k2,2n "$f" > "${f%.bedGraph}.sorted.bg"
  "${BEDGRAPH_TO_BIGWIG_BIN}" "${f%.bedGraph}.sorted.bg" "${path_refCHROMSIZES}" "${f%.bg}.bw"
done

rm -f "${outdir}/${sample_name}"*.bg \
      "${outdir}/${sample_name}"*.bg

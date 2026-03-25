#!/bin/bash
# Usage:
#   ./RNA_COVERAGE.sh SAMPLE_NAME FILTERED_BAM REF_BASE_DIR OUTDIR THREADS SPECIES
# SPECIES: human | mouse

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 SAMPLE_NAME FILTERED_BAM REF_BASE_DIR OUTDIR THREADS SPECIES(human|mouse)" >&2
    exit 1
fi

sample_name="${1}"
INBAM="${2}"
path_ref="${3}"
outdir="${4}"
threads="${5}"
species="${6}"
STAR_BIN="${STAR_BIN:-STAR}"
BEDGRAPH_TO_BIGWIG_BIN="${BEDGRAPH_TO_BIGWIG_BIN:-bedGraphToBigWig}"

path_refDB="${path_ref}/GRCh38_TrES/star"
path_refCHROMSIZES="${path_ref}/hg38.chrom.sizes"

if [[ "${species}" == "mouse" ]]; then
    path_refDB="${path_ref}/GRCm39_TrES/star"
    path_refCHROMSIZES="${path_ref}/mm39.chrom.sizes"
elif [[ "${species}" != "human" ]]; then
    echo "ERROR: SPECIES must be 'human' or 'mouse' (got: ${species})" >&2
    exit 1
fi

if [[ ! -s "${INBAM}" ]]; then
    echo "ERROR: Input filtered BAM missing or empty: ${INBAM}" >&2
    exit 1
fi

echo "Using STAR_BIN=${STAR_BIN}"
echo "Using BEDGRAPH_TO_BIGWIG_BIN=${BEDGRAPH_TO_BIGWIG_BIN}"

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

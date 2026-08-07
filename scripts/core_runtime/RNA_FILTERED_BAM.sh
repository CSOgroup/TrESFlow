#!/bin/bash
# Usage:
#   ./RNA_FILTERED_BAM.sh SAMPLE_NAME SOLO_DIR ALIGNED_BAM CANONICAL_CONTIGS OUTDIR THREADS

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 SAMPLE_NAME SOLO_DIR ALIGNED_BAM CANONICAL_CONTIGS OUTDIR THREADS" >&2
    exit 1
fi

sample_name="${1}"
solo_dir="${2}"
INBAM="${3}"
canonical_contigs="${4}"
outdir="${5}"
threads="${6}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

BARCODES="${solo_dir}/filtered/barcodes.tsv"
OUTBAM="${outdir}/${sample_name}.filtered_cells.internal.bam"

if [[ ! -s "${BARCODES}" ]]; then
    echo "ERROR: Filtered STARsolo barcodes missing or empty: ${BARCODES}" >&2
    exit 1
fi

if [[ ! -s "${INBAM}" ]]; then
    echo "ERROR: Input aligned BAM missing or empty: ${INBAM}" >&2
    exit 1
fi

echo "Using SAMTOOLS_BIN=${SAMTOOLS_BIN}"

bash "${script_dir}/FilterCanonicalBam.sh" \
    "${INBAM}" \
    "${OUTBAM}" \
    "${canonical_contigs}" \
    "${threads}" \
    uncompressed \
    --exclude-flags 0x100 \
    --require-flags 0x1,0x2 \
    --tag-file RG:"${BARCODES}"

rm "${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"*

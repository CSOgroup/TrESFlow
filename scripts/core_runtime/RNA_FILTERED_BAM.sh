#!/bin/bash
# Usage:
#   ./RNA_FILTERED_BAM.sh SAMPLE_NAME SOLO_DIR ALIGNED_BAM OUTDIR THREADS

set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 SAMPLE_NAME SOLO_DIR ALIGNED_BAM OUTDIR THREADS" >&2
    exit 1
fi

sample_name="${1}"
solo_dir="${2}"
INBAM="${3}"
outdir="${4}"
threads="${5}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

BARCODES="${solo_dir}/filtered/barcodes.tsv"
OUTBAM="${outdir}/${sample_name}.filtered_cells.bam"

if [[ ! -s "${BARCODES}" ]]; then
    echo "ERROR: Filtered STARsolo barcodes missing or empty: ${BARCODES}" >&2
    exit 1
fi

if [[ ! -s "${INBAM}" ]]; then
    echo "ERROR: Input aligned BAM missing or empty: ${INBAM}" >&2
    exit 1
fi

echo "Using SAMTOOLS_BIN=${SAMTOOLS_BIN}"

"${SAMTOOLS_BIN}" view --threads "${threads}" --with-header --exclude-flags 0x100 --require-flags 0x1,0x2 --tag-file RG:"${BARCODES}" --output "${OUTBAM}" "${INBAM}"
"${SAMTOOLS_BIN}" index --threads "${threads}" "${INBAM}"

rm "${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"*

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
PYTHON3_BIN="${PYTHON3_BIN:-python3}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

BARCODES="${solo_dir}/filtered/barcodes.tsv"
OUTBAM="${outdir}/${sample_name}.filtered_cells.internal.bam"
RETENTION_METRICS="${outdir}/${sample_name}.rna_filter_retention.tsv"

if [[ ! -s "${BARCODES}" ]]; then
    echo "ERROR: Filtered STARsolo barcodes missing or empty: ${BARCODES}" >&2
    exit 1
fi

if [[ ! -s "${INBAM}" ]]; then
    echo "ERROR: Input aligned BAM missing or empty: ${INBAM}" >&2
    exit 1
fi

echo "Using SAMTOOLS_BIN=${SAMTOOLS_BIN}"

# Audit the same nested predicates used by the existing final filter in one
# extra sequential read of the transient STAR BAM. This produces counts only;
# it does not feed or alter the BAM data path.
"${SAMTOOLS_BIN}" view "${INBAM}" \
  | "${PYTHON3_BIN}" "${script_dir}/SummarizeRnaRetention.py" \
      --split-id "${sample_name}" \
      --canonical-contigs "${canonical_contigs}" \
      --called-barcodes "${BARCODES}" \
      --output "${RETENTION_METRICS}"

bash "${script_dir}/FilterCanonicalBam.sh" \
    "${INBAM}" \
    "${OUTBAM}" \
    "${canonical_contigs}" \
    "${threads}" \
    uncompressed \
    --exclude-flags 0x100 \
    --require-flags 0x1,0x2 \
    --tag-file RG:"${BARCODES}"

expected_pairs="$(awk -F '\t' '$2 == "called_cell_pairs" { print $3 }' "${RETENTION_METRICS}")"
observed_pairs="$("${SAMTOOLS_BIN}" view --count --require-flags 0x40 --exclude-flags 0x900 "${OUTBAM}")"
if [[ -z "${expected_pairs}" || "${observed_pairs}" != "${expected_pairs}" ]]; then
    echo "ERROR: RNA final pair count does not match retention audit for ${sample_name}: expected=${expected_pairs:-missing}, observed=${observed_pairs}" >&2
    exit 1
fi

rm "${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"*

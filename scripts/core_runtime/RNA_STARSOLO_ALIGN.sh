#!/bin/bash
# Usage:
#   ./RNA_STARSOLO_ALIGN.sh SAMPLE_NAME TAGGED.usam STAR_INDEX_DIR OUTDIR THREADS

set -euo pipefail

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 SAMPLE_NAME TAGGED.usam STAR_INDEX_DIR OUTDIR THREADS" >&2
    exit 1
fi

sample_name="${1}"
USAM_IN="${2}"
path_refDB="${3}"
outdir="${4}"
threads="${5}"
STAR_BIN="${STAR_BIN:-STAR}"

if [[ ! -d "${path_refDB}" ]]; then
    echo "ERROR: STAR index directory missing: ${path_refDB}" >&2
    exit 1
fi

if [[ ! -s "${USAM_IN}" ]]; then
    echo "ERROR: Input SAM missing or empty: ${USAM_IN}" >&2
    exit 1
fi

UMIlen=10

# String mode reads XI and UM independently. soloCBlen remains a valid
# nucleotide-mode placeholder (1..31), never the namespaced string length.
echo "Using STAR index directory=${path_refDB}; full-cell-v1 String barcode XI; separate UM (${UMIlen} nt)"
echo "Using STAR_BIN=${STAR_BIN}"

ulimit -n 32000 2>/dev/null || echo "WARNING: using ulimit -n = $(ulimit -n)" >&2

"${STAR_BIN}" \
  --genomeDir "${path_refDB}" \
  --runThreadN "${threads}" \
  --readFilesIn "${USAM_IN}" \
  --readFilesType SAM PE \
  --twopassMode Basic \
  --outFilterType BySJout \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000 \
  --sjdbScore 2 \
  --limitSjdbInsertNsj 4000000 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --outFilterScoreMinOverLread 0.33 \
  --outFilterMatchNminOverLread 0.33 \
  --outFilterMultimapNmax 20 \
  --outSAMstrandField intronMotif \
  --outFilterIntronMotifs RemoveNoncanonical \
  --soloType CB_UMI_Simple \
  --soloCBtype String \
  --soloInputSAMattrBarcodeSeq XI UM \
  --soloInputSAMattrBarcodeQual - \
  --soloCBwhitelist None \
  --soloCBstart 1 --soloCBlen 24 \
  --soloUMIstart 1 --soloUMIlen "${UMIlen}" \
  --soloBarcodeReadLength 0 \
  --soloFeatures GeneFull \
  --soloStrand Forward \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR 15000 0.99 30 20000 90000 50 0.01 20000 0.05 10000 \
  --soloCellReadStats Standard \
  --outSAMtype BAM SortedByCoordinate \
  --outBAMcompression 0 \
  --outSAMattributes NH HI nM AS UB GX GN NM MD jM jI MC  \
  --outSAMunmapped None \
  --soloOutFileNames "Solo.out" "features.tsv" "barcodes.tsv" "matrix.mtx" \
  --outFileNamePrefix "${outdir}/${sample_name}."

# STAR's SAM transport preserves XI/CB on all alignment records, while UB output
# also emits CB. Validate identity and retain one copy of each cell attribute.
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
PYTHON3_BIN="${PYTHON3_BIN:-python3}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
aligned_bam="${outdir}/${sample_name}.Aligned.sortedByCoord.out.bam"
"${SAMTOOLS_BIN}" view -h "${aligned_bam}" \
  | "${PYTHON3_BIN}" "${script_dir}/NormalizeRnaBamTags.py" \
  | "${SAMTOOLS_BIN}" view -u -o "${aligned_bam}.identity.tmp.bam" -
mv "${aligned_bam}.identity.tmp.bam" "${aligned_bam}"

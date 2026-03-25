#!/bin/bash
# Usage:
#   ./RNA_STARSOLO_ALIGN.sh SAMPLE_NAME TAGGED.usam REF_BASE_DIR OUTDIR THREADS SPECIES
# SPECIES: human | mouse

set -euo pipefail

if [[ $# -lt 6 ]]; then
    echo "Usage: $0 SAMPLE_NAME TAGGED.usam REF_BASE_DIR OUTDIR THREADS SPECIES(human|mouse)" >&2
    exit 1
fi

sample_name="${1}"
USAM_IN="${2}"
path_ref="${3}"
outdir="${4}"
threads="${5}"
species="${6}"
STAR_BIN="${STAR_BIN:-STAR}"

path_refDB="${path_ref}/GRCh38_TrES/star"

if [[ "${species}" == "mouse" ]]; then
    path_refDB="${path_ref}/GRCm39_TrES/star"
elif [[ "${species}" != "human" ]]; then
    echo "ERROR: SPECIES must be 'human' or 'mouse' (got: ${species})" >&2
    exit 1
fi

if [[ ! -s "${USAM_IN}" ]]; then
    echo "ERROR: Input SAM missing or empty: ${USAM_IN}" >&2
    exit 1
fi

UMIlen=10

CBlen="$(
  awk -v UMIlen="${UMIlen}" '
    BEGIN { FS="\t" }
    $0 ~ /^@/ { next }
    {
      cb = ""; cr = "";
      for (i=12; i<=NF; i++) {
        if ($i ~ /^CB:Z:/) { cb = substr($i,6); break }
      }
      if (cb != "") { print length(cb); exit }

      for (i=12; i<=NF; i++) {
        if ($i ~ /^CR:Z:/) { cr = substr($i,6); break }
      }
      if (cr != "") { print length(cr) - UMIlen; exit }
    }
  ' "${USAM_IN}"
)"

if [[ -z "${CBlen}" ]] || [[ "${CBlen}" -le 0 ]]; then
  echo "ERROR: Could not detect CB length from ${USAM_IN} (no CB:Z: or CR:Z: found)." >&2
  exit 1
fi

UMIstart=$(( CBlen + 1 ))

echo "Using species=${species}"
echo "Detected CBlen=${CBlen} => UMIstart=${UMIstart} UMIlen=${UMIlen}"
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
  --soloInputSAMattrBarcodeSeq CR \
  --soloInputSAMattrBarcodeQual - \
  --soloCBwhitelist None \
  --soloCBstart 1 --soloCBlen "${CBlen}" \
  --soloUMIstart "${UMIstart}" --soloUMIlen "${UMIlen}" \
  --soloBarcodeReadLength 0 \
  --soloFeatures GeneFull \
  --soloStrand Forward \
  --soloMultiMappers EM \
  --soloUMIdedup 1MM_CR \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloCellFilter EmptyDrops_CR 15000 0.99 30 20000 90000 50 0.01 20000 0.05 10000 \
  --soloCellReadStats Standard \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes NH HI nM AS CB UB GX GN NM MD jM jI MC  \
  --outSAMunmapped None \
  --soloOutFileNames "Solo.out" "features.tsv" "barcodes.tsv" "matrix.mtx" \
  --outFileNamePrefix "${outdir}/${sample_name}."

#!/usr/bin/env bash
# Usage: SplitDuplicatesDNA.sh INPUT.bam OUTPUT.bam OUTPUT.bai MAPPED_READS.txt WARNING.tsv THREADS SAMPLE GROUP MARK SPLIT

set -euo pipefail

if [[ $# -ne 10 ]]; then
  echo "Usage: $0 INPUT.bam OUTPUT.bam OUTPUT.bai MAPPED_READS.txt WARNING.tsv THREADS SAMPLE GROUP MARK SPLIT" >&2
  exit 1
fi

input_bam="${1}"
output_bam="${2}"
output_bai="${3}"
mapped_reads_file="${4}"
warning_file="${5}"
threads="${6}"
sample_id="${7}"
group_name="${8}"
mark_name="${9}"
split_name="${10}"

if [[ ! -s "${input_bam}" ]]; then
  echo "ERROR: Duplicate-marked BAM is missing or empty: ${input_bam}" >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: Missing samtools on task PATH" >&2
  exit 1
fi

# NORMALIZE_DNA_MARKDUPLICATES is the only producer of this input and has
# already restricted it to the resolved canonical chromosomes. Preserve that
# header and coordinate order while removing only records carrying 0x400.
samtools view \
  --threads "${threads}" \
  --bam \
  --with-header \
  --exclude-flags 0x400 \
  --output "${output_bam}" \
  "${input_bam}"

samtools index \
  --threads "${threads}" \
  --bai \
  --output "${output_bai}" \
  "${output_bam}"

mapped_reads="$(samtools view --threads "${threads}" -c -F 4 "${output_bam}")"
printf '%s\n' "${mapped_reads}" > "${mapped_reads_file}"

if [[ "${mapped_reads}" -eq 0 ]]; then
  bam_path="$(readlink -f "${output_bam}")"
  cat >&2 <<'EOF'
================================================================================
WARNING: ZERO MAPPED READS IN DNA NoDup BAM
================================================================================
EOF
  echo "Sample: ${sample_id}" >&2
  echo "Group: ${group_name}" >&2
  echo "Mark: ${mark_name}" >&2
  echo "BAM: ${bam_path}" >&2
  echo "Mapped reads: ${mapped_reads}" >&2
  echo "Skipped nf-core/deeptools/bamcoverage for ${split_name}" >&2
  printf 'sample\tgroup\tmark\tbam\tmapped_reads\tskipped_output\n' > "${warning_file}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample_id}" \
    "${group_name}" \
    "${mark_name}" \
    "${bam_path}" \
    "${mapped_reads}" \
    "${split_name}_NoDup.bw" \
    >> "${warning_file}"
fi

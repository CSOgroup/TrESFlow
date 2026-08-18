#!/usr/bin/env bash
# Usage:
#   FilterCanonicalBam.sh INPUT.bam OUTPUT.bam ALLOWLIST.txt THREADS normal|uncompressed [--validation-summary PATH] [samtools-view-options...]

set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 INPUT.bam OUTPUT.bam ALLOWLIST.txt THREADS normal|uncompressed [--validation-summary PATH] [samtools-view-options...]" >&2
  exit 1
fi

input_bam="${1}"
output_bam="${2}"
allowlist="${3}"
threads="${4}"
compression_mode="${5}"
shift 5

validation_summary=""
view_options=()

while (( $# > 0 )); do
  case "${1}" in
    --validation-summary)
      if (( $# < 2 )); then
        echo "ERROR: --validation-summary requires a path" >&2
        exit 1
      fi
      validation_summary="${2}"
      shift 2
      ;;
    *)
      view_options+=("${1}")
      shift
      ;;
  esac
done

SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

validation_tmp=""
header_file=""
reheadered_bam=""

cleanup() {
  rm -f \
    "${validation_tmp:-}" \
    "${header_file:-}" \
    "${reheadered_bam:-}"
}
trap cleanup EXIT

if [[ ! -s "${input_bam}" ]]; then
  echo "ERROR: Canonical BAM filter input is missing or empty: ${input_bam}" >&2
  exit 1
fi
if [[ ! -s "${allowlist}" ]]; then
  echo "ERROR: Canonical chromosome allowlist is missing or empty: ${allowlist}" >&2
  exit 1
fi

mapfile -t canonical_contigs < <(awk 'NF && $1 !~ /^#/ { print $1 }' "${allowlist}")
if [[ ${#canonical_contigs[@]} -eq 0 ]]; then
  echo "ERROR: Canonical chromosome allowlist has no usable entries: ${allowlist}" >&2
  exit 1
fi

declare -A seen_contigs=()
for contig in "${canonical_contigs[@]}"; do
  if [[ -n "${seen_contigs[${contig}]:-}" ]]; then
    echo "ERROR: Duplicate chromosome '${contig}' in canonical allowlist ${allowlist}" >&2
    exit 1
  fi
  seen_contigs["${contig}"]=1
done

declare -A header_order=()
header_index=0
while IFS= read -r contig; do
  header_order["${contig}"]="${header_index}"
  header_index=$((header_index + 1))
done < <("${SAMTOOLS_BIN}" view -H "${input_bam}" | awk -F '\t' '
  $1 == "@SQ" {
    for (i = 2; i <= NF; i++) {
      if ($i ~ /^SN:/) {
        sub(/^SN:/, "", $i)
        print $i
        break
      }
    }
  }
')

previous_index=-1
for contig in "${canonical_contigs[@]}"; do
  if [[ -z "${header_order[${contig}]+present}" ]]; then
    echo "ERROR: Canonical chromosome '${contig}' is absent from BAM header ${input_bam}" >&2
    exit 1
  fi
  current_index="${header_order[${contig}]}"
  if (( current_index <= previous_index )); then
    echo "ERROR: Canonical allowlist order does not match BAM @SQ order at '${contig}'" >&2
    exit 1
  fi
  previous_index="${current_index}"
done

if [[ ! -s "${input_bam}.bai" && ! -s "${input_bam%.bam}.bai" && ! -s "${input_bam}.csi" ]]; then
  "${SAMTOOLS_BIN}" index --threads "${threads}" "${input_bam}"
fi

case "${compression_mode}" in
  normal)
    output_option=(--bam)
    ;;
  uncompressed)
    output_option=(--uncompressed)
    ;;
  *)
    echo "ERROR: Unsupported BAM compression mode '${compression_mode}'; expected normal or uncompressed" >&2
    exit 1
    ;;
esac

"${SAMTOOLS_BIN}" view \
  --threads "${threads}" \
  "${output_option[@]}" \
  --with-header \
  "${view_options[@]}" \
  --output "${output_bam}" \
  "${input_bam}" \
  "${canonical_contigs[@]}"

# Validate the retained records once.
#
# This single scan performs the same two validations that previously required
# two independent full BAM scans:
#   1. every retained alignment RNAME is canonical;
#   2. determine whether every retained RNEXT is canonical (or '='/'*').
#
# It also counts primary, non-supplementary R1 records. RNA_FILTERED_BAM uses
# this count for its existing retention reconciliation, avoiding a third scan.
#
# IMPORTANT: this is validation only. It does not participate in filtering and
# therefore cannot alter which records are retained.
validation_tmp="$(mktemp "${TMPDIR:-/tmp}/canonical-validation.XXXXXX.tsv")"

"${SAMTOOLS_BIN}" view "${output_bam}" | awk -v allowlist="${allowlist}" '
  BEGIN {
    while ((getline line < allowlist) > 0) {
      split(line, fields)
      if (fields[1] != "" && fields[1] !~ /^#/) allowed[fields[1]] = 1
    }

    invalid_alignment_records = 0
    invalid_mate_reference_records = 0
    primary_r1_records = 0
  }

  {
    if ($3 != "*" && !($3 in allowed)) {
      print "ERROR: Noncanonical alignment record escaped filtering: " $1 " on " $3 > "/dev/stderr"
      invalid_alignment_records++
    }

    if ($7 != "=" && $7 != "*" && !($7 in allowed)) {
      invalid_mate_reference_records++
    }

    # Equivalent to:
    #   samtools view --count --require-flags 0x40 --exclude-flags 0x900
    #
    # Use arithmetic bit tests for portable awk rather than relying on
    # implementation-specific bitwise functions.
    flag = $2 + 0
    is_r1 = int(flag / 64) % 2
    is_secondary = int(flag / 256) % 2
    is_supplementary = int(flag / 2048) % 2

    if (is_r1 && !is_secondary && !is_supplementary) {
      primary_r1_records++
    }
  }

  END {
    print "invalid_alignment_records\t" invalid_alignment_records
    print "invalid_mate_reference_records\t" invalid_mate_reference_records
    print "primary_r1_records\t" primary_r1_records
  }
' > "${validation_tmp}"

invalid_alignment_records="$(
  awk -F '\t' '$1 == "invalid_alignment_records" { print $2 }' "${validation_tmp}"
)"
invalid_mate_reference_records="$(
  awk -F '\t' '$1 == "invalid_mate_reference_records" { print $2 }' "${validation_tmp}"
)"
primary_r1_records="$(
  awk -F '\t' '$1 == "primary_r1_records" { print $2 }' "${validation_tmp}"
)"

for value_name in \
  invalid_alignment_records \
  invalid_mate_reference_records \
  primary_r1_records
do
  value="${!value_name}"
  if [[ ! "${value}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Invalid canonical-filter validation result ${value_name}='${value}'" >&2
    exit 1
  fi
done

if (( invalid_alignment_records != 0 )); then
  echo "ERROR: ${invalid_alignment_records} noncanonical alignment record(s) escaped filtering in ${output_bam}" >&2
  exit 1
fi

if [[ -n "${validation_summary}" ]]; then
  {
    printf 'metric\tcount\n'
    printf 'invalid_alignment_records\t%s\n' "${invalid_alignment_records}"
    printf 'invalid_mate_reference_records\t%s\n' "${invalid_mate_reference_records}"
    printf 'primary_r1_records\t%s\n' "${primary_r1_records}"
  } > "${validation_summary}"
fi

# A canonical-only @SQ dictionary is safe when every retained record's mate
# reference is also canonical (or uses '='/'*'). Otherwise keep the unused
# noncanonical @SQ entries so RNEXT remains valid, but never keep noncanonical
# alignment records.
if (( invalid_mate_reference_records == 0 )); then
  header_file="$(mktemp "${TMPDIR:-/tmp}/canonical-header.XXXXXX.sam")"
  reheadered_bam="$(mktemp "${TMPDIR:-/tmp}/canonical-reheader.XXXXXX.bam")"

  "${SAMTOOLS_BIN}" view -H "${output_bam}" | awk -v allowlist="${allowlist}" '
    BEGIN {
      FS = OFS = "\t"
      while ((getline line < allowlist) > 0) {
        split(line, fields)
        if (fields[1] != "" && fields[1] !~ /^#/) allowed[fields[1]] = 1
      }
    }
    $1 == "@SQ" {
      name = ""
      for (i = 2; i <= NF; i++) {
        if ($i ~ /^SN:/) {
          name = $i
          sub(/^SN:/, "", name)
          break
        }
      }
      if (name in allowed) print
      next
    }
    { print }
  ' > "${header_file}"

  "${SAMTOOLS_BIN}" reheader -P "${header_file}" "${output_bam}" > "${reheadered_bam}"
  mv "${reheadered_bam}" "${output_bam}"

  rm -f "${header_file}"
  header_file=""
  reheadered_bam=""
else
  echo "WARNING: Retaining unused noncanonical @SQ entries in ${output_bam} because a canonical alignment has a mate reference outside the allowlist" >&2
fi

"${SAMTOOLS_BIN}" quickcheck -v "${output_bam}"

rm -f "${validation_tmp}"
validation_tmp=""
trap - EXIT

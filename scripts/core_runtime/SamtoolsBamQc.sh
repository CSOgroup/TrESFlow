#!/usr/bin/env bash
# Usage: SamtoolsBamQc.sh INPUT.bam OUTPUT_PREFIX RUN_IDXSTATS THREADS

set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 INPUT.bam OUTPUT_PREFIX RUN_IDXSTATS THREADS" >&2
  exit 1
fi

input_bam="${1}"
output_prefix="${2}"
run_idxstats="${3}"
threads="${4}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

if [[ ! -e "${input_bam}" ]]; then
  echo "ERROR: Samtools QC input is missing: ${input_bam}" >&2
  exit 1
fi
if ! command -v "${SAMTOOLS_BIN}" >/dev/null 2>&1; then
  echo "ERROR: Missing configured Samtools QC executable: ${SAMTOOLS_BIN}" >&2
  exit 1
fi

set +e
"${SAMTOOLS_BIN}" quickcheck "${input_bam}"
quickcheck_exit=$?
set -e
printf 'id\tbam\texit_code\tstatus\n' > "${output_prefix}.quickcheck.tsv"
if [[ "${quickcheck_exit}" == "0" ]]; then
  quickcheck_status="pass"
else
  quickcheck_status="fail"
fi
printf '%s\t%s\t%s\t%s\n' \
  "${output_prefix}" \
  "${input_bam}" \
  "${quickcheck_exit}" \
  "${quickcheck_status}" \
  >> "${output_prefix}.quickcheck.tsv"

set +e
"${SAMTOOLS_BIN}" flagstat \
  --threads "${threads}" \
  "${input_bam}" \
  > "${output_prefix}.flagstat"
flagstat_exit=$?

"${SAMTOOLS_BIN}" stats \
  --threads "${threads}" \
  "${input_bam}" \
  > "${output_prefix}.stats"
stats_exit=$?

idxstats_exit=0
case "${run_idxstats}" in
  true)
    idxstats_threads=$((threads > 0 ? threads - 1 : 0))
    "${SAMTOOLS_BIN}" idxstats \
      --threads "${idxstats_threads}" \
      "${input_bam}" \
      > "${output_prefix}.idxstats"
    idxstats_exit=$?
    ;;
  false)
    ;;
  *)
    echo "ERROR: RUN_IDXSTATS must be true or false, received '${run_idxstats}'" >&2
    exit 1
    ;;
esac
set -e

if (( flagstat_exit != 0 || stats_exit != 0 || idxstats_exit != 0 )); then
  echo "ERROR: Samtools BAM QC command failure for ${input_bam}: flagstat=${flagstat_exit}, stats=${stats_exit}, idxstats=${idxstats_exit}" >&2
  exit 1
fi

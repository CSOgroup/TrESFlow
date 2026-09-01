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

echo "Using STAR=$(command -v STAR)"
echo "Using bedGraphToBigWig=$(command -v bedGraphToBigWig)"
echo "Using STAR index directory=${path_refDB}"
echo "Using chromosome sizes=${path_refCHROMSIZES}"

if [[ ! "${threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: THREADS must be a positive integer; received: ${threads}" >&2
    exit 1
fi

# STAR's inputAlignmentsFromBAM signal-generation path is effectively serial.
# Keep each STAR invocation at one thread and use the process-level CPU
# allocation to run the independent stranded/unstranded jobs concurrently.
STAR_COVERAGE_THREADS=1

run_stranded_star() {
    STAR \
      --runMode inputAlignmentsFromBAM \
      --runThreadN "${STAR_COVERAGE_THREADS}" \
      --genomeDir "${path_refDB}" \
      --inputBAMfile "${INBAM}" \
      --outWigType bedGraph \
      --outWigStrand Stranded \
      --outWigNorm RPM \
      --outFileNamePrefix "${outdir}/${sample_name}.stranded_"
}

run_unstranded_star() {
    STAR \
      --runMode inputAlignmentsFromBAM \
      --runThreadN "${STAR_COVERAGE_THREADS}" \
      --genomeDir "${path_refDB}" \
      --inputBAMfile "${INBAM}" \
      --outWigType bedGraph \
      --outWigStrand Unstranded \
      --outWigNorm RPM \
      --outFileNamePrefix "${outdir}/${sample_name}.unstranded_"
}

wait_for_jobs() {
    local status=0
    local pid

    for pid in "$@"; do
        if ! wait "${pid}"; then
            status=1
        fi
    done

    return "${status}"
}

echo "RNA coverage task CPUs=${threads}; STAR coverage threads per invocation=${STAR_COVERAGE_THREADS}"

# Preserve correct behavior even if a user constrains the entire pipeline to
# one CPU. With the normal >=2 CPU allocation, run the two independent STAR
# signal-generation jobs concurrently.
if (( threads >= 2 )); then
    run_stranded_star &
    star_pid_stranded=$!

    run_unstranded_star &
    star_pid_unstranded=$!

    if ! wait_for_jobs "${star_pid_stranded}" "${star_pid_unstranded}"; then
        echo "ERROR: One or more STAR RNA coverage jobs failed." >&2
        exit 1
    fi
else
    run_stranded_star
    run_unstranded_star
fi

# Exactly the same sort and bedGraphToBigWig operations as before, but the
# three independent output tracks may now be converted concurrently.
shopt -s nullglob
stranded_bg=( "${outdir}/${sample_name}.stranded_Signal.Unique.str"*.bg )
unstranded_bg=( "${outdir}/${sample_name}.unstranded_Signal.Unique.str"*.bg )
shopt -u nullglob

if (( ${#stranded_bg[@]} != 2 )); then
    echo "ERROR: Expected exactly 2 stranded unique bedGraph files; found ${#stranded_bg[@]}." >&2
    exit 1
fi

if (( ${#unstranded_bg[@]} != 1 )); then
    echo "ERROR: Expected exactly 1 unstranded unique bedGraph file; found ${#unstranded_bg[@]}." >&2
    exit 1
fi

coverage_bg=(
    "${stranded_bg[@]}"
    "${unstranded_bg[@]}"
)

convert_bedgraph() {
    local f="$1"

    sort -k1,1 -k2,2n "$f" > "${f%.bedGraph}.sorted.bg"
    bedGraphToBigWig \
        "${f%.bedGraph}.sorted.bg" \
        "${path_refCHROMSIZES}" \
        "${f%.bg}.bw"
}

if (( threads >= 3 )); then
    conversion_pids=()

    for f in "${coverage_bg[@]}"; do
        convert_bedgraph "$f" &
        conversion_pids+=( "$!" )
    done

    if ! wait_for_jobs "${conversion_pids[@]}"; then
        echo "ERROR: One or more RNA BigWig conversions failed." >&2
        exit 1
    fi

elif (( threads == 2 )); then
    convert_bedgraph "${coverage_bg[0]}" &
    conversion_pid_1=$!

    convert_bedgraph "${coverage_bg[1]}" &
    conversion_pid_2=$!

    if ! wait_for_jobs "${conversion_pid_1}" "${conversion_pid_2}"; then
        echo "ERROR: RNA BigWig conversion failed." >&2
        exit 1
    fi

    convert_bedgraph "${coverage_bg[2]}"

else
    for f in "${coverage_bg[@]}"; do
        convert_bedgraph "$f"
    done
fi

rm -f "${outdir}/${sample_name}"*.bg

#!/usr/bin/env bash

set -euo pipefail

repo_root="${TRESFLOW_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)}"
samtools_bin="${SAMTOOLS_BIN:-/home/annan/micromamba/envs/tres/bin/samtools}"
bamcoverage_bin="${BAMCOVERAGE_BIN:-/home/annan/micromamba/envs/tres/bin/bamCoverage}"
python_bin="${PYTHON3_BIN:-/home/annan/micromamba/envs/tres/bin/python3}"
star_bin="${STAR_BIN:-/home/annan/micromamba/envs/tres/bin/STAR}"
bedgraph_to_bigwig_bin="${BEDGRAPH_TO_BIGWIG_BIN:-/home/annan/micromamba/envs/tres/bin/bedGraphToBigWig}"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/tresflow-canonical-bam.XXXXXX")"

cleanup() {
  if [[ "${KEEP_CANONICAL_TEST_TMP:-0}" == "1" ]]; then
    printf 'retained canonical test directory: %s\n' "${work_dir}" >&2
  else
    rm -rf "${work_dir}"
  fi
}
trap cleanup EXIT

cat > "${work_dir}/ucsc.sam" <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chr2	LN:10000
@SQ	SN:chrX	LN:10000
@SQ	SN:chrY	LN:10000
@SQ	SN:chrM	LN:10000
@SQ	SN:chrUn_KI270442v1	LN:10000
@SQ	SN:chr1_KI270706v1_random	LN:10000
@SQ	SN:chrEBV	LN:10000
@RG	ID:cells	SM:synthetic
auto	99	chr1	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
auto	147	chr1	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
duplicate	1123	chr1	50	60	10M	=	70	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
duplicate	1171	chr1	70	60	10M	=	50	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
xread	83	chrX	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
xread	163	chrX	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
yread	99	chrY	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
yread	147	chrY	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
mito	99	chrM	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
mito	147	chrM	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
unplaced	99	chrUn_KI270442v1	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
unplaced	147	chrUn_KI270442v1	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
alternative	99	chr1_KI270706v1_random	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
alternative	147	chr1_KI270706v1_random	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
decoy	99	chrEBV	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
decoy	147	chrEBV	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
EOF

sed \
  -e 's/chrM/MT/g' \
  -e 's/chrX/X/g' \
  -e 's/chrY/Y/g' \
  -e 's/chrUn/Un/g' \
  -e 's/chrEBV/EBV/g' \
  -e 's/chr1/1/g' \
  -e 's/chr2/2/g' \
  "${work_dir}/ucsc.sam" > "${work_dir}/ensembl.sam"

run_case() {
  style="${1}"
  shift
  canonical_contigs=("$@")
  source_sam="${work_dir}/${style}.sam"
  source_bam="${work_dir}/${style}.source.bam"
  allowlist="${work_dir}/${style}.canonical.txt"
  chrom_sizes="${work_dir}/${style}.canonical.chrom.sizes"
  dna_bam="${work_dir}/${style}.NoDup.bam"
  rna_internal="${work_dir}/${style}.rna.internal.bam"
  rna_published="${work_dir}/${style}.rna.published.bam"

  awk '/^@/ { print; next } { print $0 "\tNH:i:1" }' "${source_sam}" \
    | "${samtools_bin}" view --bam - \
    | "${samtools_bin}" sort -o "${source_bam}" -
  "${samtools_bin}" index "${source_bam}"
  printf '%s\n' "${canonical_contigs[@]}" > "${allowlist}"
  for contig in "${canonical_contigs[@]}"; do
    printf '%s\t10000\n' "${contig}"
  done > "${chrom_sizes}"

  SAMTOOLS_BIN="${samtools_bin}" \
    bash "${repo_root}/scripts/core_runtime/FilterCanonicalBam.sh" \
      "${source_bam}" "${dna_bam}" "${allowlist}" 1 normal \
      --exclude-flags 0x400
  "${samtools_bin}" index "${dna_bam}"
  "${samtools_bin}" quickcheck "${dna_bam}"

  "${samtools_bin}" view --exclude-flags 0x400 "${source_bam}" \
    "${canonical_contigs[@]}" > "${work_dir}/${style}.dna.expected.sam"
  "${samtools_bin}" view "${dna_bam}" > "${work_dir}/${style}.dna.observed.sam"
  cmp "${work_dir}/${style}.dna.expected.sam" "${work_dir}/${style}.dna.observed.sam"

  SAMTOOLS_BIN="${samtools_bin}" \
    bash "${repo_root}/scripts/core_runtime/FilterCanonicalBam.sh" \
      "${source_bam}" "${rna_internal}" "${allowlist}" 1 uncompressed \
      --exclude-flags 0x100 --require-flags 0x1,0x2
  "${samtools_bin}" quickcheck "${rna_internal}"
  "${samtools_bin}" index "${rna_internal}"
  "${samtools_bin}" view --bam --with-header \
    --output "${rna_published}" "${rna_internal}"
  "${samtools_bin}" index "${rna_published}"
  "${samtools_bin}" quickcheck "${rna_published}"

  "${samtools_bin}" view "${source_bam}" "${canonical_contigs[@]}" \
    > "${work_dir}/${style}.rna.expected.sam"
  "${samtools_bin}" view "${rna_internal}" > "${work_dir}/${style}.rna.internal.sam"
  "${samtools_bin}" view "${rna_published}" > "${work_dir}/${style}.rna.published.sam"
  cmp "${work_dir}/${style}.rna.expected.sam" "${work_dir}/${style}.rna.internal.sam"
  cmp "${work_dir}/${style}.rna.internal.sam" "${work_dir}/${style}.rna.published.sam"

  for bam in "${dna_bam}" "${rna_internal}" "${rna_published}"; do
    mapfile -t header_contigs < <("${samtools_bin}" view -H "${bam}" \
      | awk -F '\t' '$1 == "@SQ" { sub(/^SN:/, "", $2); print $2 }')
    [[ "${header_contigs[*]}" == "${canonical_contigs[*]}" ]]
    for contig in "${canonical_contigs[@]}"; do
      "${samtools_bin}" view --count "${bam}" "${contig}" >/dev/null
    done
    [[ $("${samtools_bin}" view --count "${bam}" "${canonical_contigs[1]}") -eq 0 ]]
    [[ $("${samtools_bin}" view --count "${bam}" "${canonical_contigs[0]}") -gt 0 ]]
    [[ $("${samtools_bin}" view --count "${bam}" "${canonical_contigs[2]}") -gt 0 ]]
    [[ $("${samtools_bin}" view --count "${bam}" "${canonical_contigs[3]}") -gt 0 ]]
    [[ $("${samtools_bin}" view --count "${bam}" "${canonical_contigs[4]}") -gt 0 ]]
  done

  MPLCONFIGDIR="${work_dir}/mpl-${style}" "${bamcoverage_bin}" \
    --bam "${dna_bam}" \
    --binSize 10 \
    --effectiveGenomeSize 50000 \
    --numberOfProcessors 1 \
    --outFileName "${work_dir}/${style}.dna.bw" \
    >/dev/null
  MPLCONFIGDIR="${work_dir}/mpl-${style}-rna" "${bamcoverage_bin}" \
    --bam "${rna_internal}" \
    --binSize 10 \
    --effectiveGenomeSize 50000 \
    --numberOfProcessors 1 \
    --outFileName "${work_dir}/${style}.rna.bw" \
    >/dev/null

  reference_fasta="${work_dir}/${style}.fa"
  star_index="${work_dir}/${style}.star"
  star_coverage_dir="${work_dir}/${style}.star-coverage"
  mkdir -p "${star_index}" "${star_coverage_dir}"
  while read -r contig length; do
    printf '>%s\n' "${contig}" >> "${reference_fasta}"
    "${python_bin}" -c \
      'import random, sys; print("".join(random.Random(sys.argv[1]).choices("ACGT", k=int(sys.argv[2]))))' \
      "${contig}" "${length}" >> "${reference_fasta}"
  done < <("${samtools_bin}" view -H "${source_bam}" | awk -F '\t' '
    $1 == "@SQ" {
      name = $2
      seq_length = $3
      sub(/^SN:/, "", name)
      sub(/^LN:/, "", seq_length)
      print name "\t" seq_length
    }
  ')
  (
    cd "${work_dir}"
    "${star_bin}" \
      --runMode genomeGenerate \
      --runThreadN 1 \
      --genomeDir "${star_index}" \
      --genomeFastaFiles "${reference_fasta}" \
      --genomeSAindexNbases 2 \
      --genomeChrBinNbits 10 \
      >/dev/null
  )
  STAR_BIN="${star_bin}" BEDGRAPH_TO_BIGWIG_BIN="${bedgraph_to_bigwig_bin}" \
    bash "${repo_root}/scripts/core_runtime/RNA_COVERAGE.sh" \
      "${style}.rna" "${rna_internal}" "${star_index}" "${chrom_sizes}" \
      "${star_coverage_dir}" 1 \
      >/dev/null
  mapfile -t star_rna_bigwigs < <(find "${star_coverage_dir}" -type f -name '*.bw' -print | sort)
  [[ ${#star_rna_bigwigs[@]} -gt 0 ]]

  "${python_bin}" - "${allowlist}" \
    "${work_dir}/${style}.dna.bw" "${work_dir}/${style}.rna.bw" \
    "${star_rna_bigwigs[@]}" <<'PY'
import sys
from pathlib import Path

import pyBigWig

allowed = set(Path(sys.argv[1]).read_text(encoding="utf-8").split())
for bigwig_path in sys.argv[2:]:
    handle = pyBigWig.open(bigwig_path)
    observed = set(handle.chroms())
    handle.close()
    invalid = observed - allowed
    if invalid:
        raise SystemExit(f"noncanonical BigWig chromosomes in {bigwig_path}: {sorted(invalid)}")
PY
}

run_case ucsc chr1 chr2 chrX chrY chrM
run_case ensembl 1 2 X Y MT

printf 'canonical BAM and BigWig synthetic tests passed\n'

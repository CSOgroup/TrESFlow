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
@SQ	SN:chrUn_KI270442v1	LN:10000
@SQ	SN:chr2	LN:10000
@SQ	SN:chrX	LN:10000
@SQ	SN:another_noncanonical_contig	LN:10000
@SQ	SN:chrY	LN:10000
@SQ	SN:chrM	LN:10000
@SQ	SN:chrEBV	LN:10000
@RG	ID:cells	SM:synthetic
auto	99	chr1	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
auto	147	chr1	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
duplicate	1123	chr1	50	60	10M	=	70	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
duplicate	1171	chr1	70	60	10M	=	50	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
crosscanonical	99	chr2	20	60	10M	chrX	50	0	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
xread	83	chrX	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
xread	163	chrX	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
crosscanonical	147	chrX	50	60	10M	chr2	20	0	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
yread	99	chrY	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
yread	147	chrY	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
mito	99	chrM	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
mito	147	chrM	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
unplaced	99	chrUn_KI270442v1	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
unplaced	147	chrUn_KI270442v1	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
alternative	99	another_noncanonical_contig	10	60	10M	=	30	30	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
alternative	147	another_noncanonical_contig	30	60	10M	=	10	-30	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
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
  source_pg_count="$(
    "${samtools_bin}" view --no-PG -H "${source_bam}" \
      | awk '$1 == "@PG" { count++ } END { print count + 0 }'
  )"
  printf '%s\n' "${canonical_contigs[@]}" > "${allowlist}"
  for contig in "${canonical_contigs[@]}"; do
    printf '%s\t10000\n' "${contig}"
  done > "${chrom_sizes}"

  if [[ "${style}" == "ucsc" ]]; then
    unsafe_filtered_bam="${work_dir}/unsafe.full-header.bam"
    unsafe_header="${work_dir}/unsafe.reduced-header.sam"
    unsafe_bam="${work_dir}/unsafe.reheader.bam"

    "${samtools_bin}" view --no-PG --bam --with-header \
      --output "${unsafe_filtered_bam}" "${source_bam}" \
      "${canonical_contigs[@]}"
    "${samtools_bin}" view -H "${unsafe_filtered_bam}" | awk -v allowlist="${allowlist}" '
      BEGIN {
        FS = OFS = "\t"
        while ((getline line < allowlist) > 0) {
          split(line, fields)
          if (fields[1] != "" && fields[1] !~ /^#/) allowed[fields[1]] = 1
        }
      }
      $1 == "@SQ" {
        name = $2
        sub(/^SN:/, "", name)
        if (name in allowed) print
        next
      }
      { print }
    ' > "${unsafe_header}"
    "${samtools_bin}" reheader -P "${unsafe_header}" "${unsafe_filtered_bam}" \
      > "${unsafe_bam}"

    # The v1.1.0 implementation looked healthy to quickcheck, but its retained
    # numeric tid/mtid values no longer matched the shortened @SQ dictionary.
    "${samtools_bin}" quickcheck "${unsafe_bam}"
    if "${samtools_bin}" view "${unsafe_bam}" >/dev/null \
      2> "${work_dir}/unsafe.view.stderr"; then
      echo "ERROR: unsafe samtools reheader unexpectedly survived a full BAM scan" >&2
      exit 1
    fi
    grep -Fq 'Numerical result out of range' "${work_dir}/unsafe.view.stderr"
    if "${samtools_bin}" index --bai --output "${unsafe_bam}.bai" \
      "${unsafe_bam}" >/dev/null 2> "${work_dir}/unsafe.index.stderr"; then
      echo "ERROR: unsafe samtools reheader unexpectedly produced a BAI" >&2
      exit 1
    fi
    grep -Fq 'Numerical result out of range' "${work_dir}/unsafe.index.stderr"
  fi

  SAMTOOLS_BIN="${samtools_bin}" \
    bash "${repo_root}/scripts/core_runtime/FilterCanonicalBam.sh" \
      "${source_bam}" "${dna_bam}" "${allowlist}" 1 normal \
      --exclude-flags 0x400
  "${samtools_bin}" quickcheck "${dna_bam}"
  "${samtools_bin}" view "${dna_bam}" >/dev/null
  "${samtools_bin}" index --bai --output "${dna_bam}.bai" "${dna_bam}"
  [[ -s "${dna_bam}.bai" ]]
  [[ $("${samtools_bin}" view --no-PG -H "${dna_bam}" \
    | awk '$1 == "@PG" { count++ } END { print count + 0 }') \
    -eq $((source_pg_count + 1)) ]]

  "${samtools_bin}" view --exclude-flags 0x400 "${source_bam}" \
    "${canonical_contigs[@]}" > "${work_dir}/${style}.dna.expected.sam"
  "${samtools_bin}" view "${dna_bam}" > "${work_dir}/${style}.dna.observed.sam"
  cmp "${work_dir}/${style}.dna.expected.sam" "${work_dir}/${style}.dna.observed.sam"

  SAMTOOLS_BIN="${samtools_bin}" \
    bash "${repo_root}/scripts/core_runtime/FilterCanonicalBam.sh" \
      "${source_bam}" "${rna_internal}" "${allowlist}" 1 uncompressed \
      --exclude-flags 0x100 --require-flags 0x1,0x2
  "${samtools_bin}" quickcheck "${rna_internal}"
  "${samtools_bin}" view "${rna_internal}" >/dev/null
  "${samtools_bin}" index --bai --output "${rna_internal}.bai" "${rna_internal}"
  [[ -s "${rna_internal}.bai" ]]
  [[ $("${samtools_bin}" view --no-PG -H "${rna_internal}" \
    | awk '$1 == "@PG" { count++ } END { print count + 0 }') \
    -eq $((source_pg_count + 1)) ]]
  "${samtools_bin}" view --bam --with-header \
    --output "${rna_published}" "${rna_internal}"
  "${samtools_bin}" quickcheck "${rna_published}"
  "${samtools_bin}" view "${rna_published}" >/dev/null
  "${samtools_bin}" index --bai --output "${rna_published}.bai" "${rna_published}"
  [[ -s "${rna_published}.bai" ]]
  [[ $("${samtools_bin}" view --no-PG -H "${rna_published}" \
    | awk '$1 == "@PG" { count++ } END { print count + 0 }') \
    -eq $((source_pg_count + 2)) ]]

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
    cmp \
      <("${samtools_bin}" view --no-PG -H "${source_bam}" | awk '$1 == "@RG"') \
      <("${samtools_bin}" view --no-PG -H "${bam}" | awk '$1 == "@RG"')
    for contig in "${canonical_contigs[@]}"; do
      "${samtools_bin}" view --count "${bam}" "${contig}" >/dev/null
    done
    for contig in "${canonical_contigs[@]}"; do
      [[ $("${samtools_bin}" view --count "${bam}" "${contig}") -gt 0 ]]
    done
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

run_noncanonical_mate_case() {
  source_sam="${work_dir}/noncanonical-mate.sam"
  source_bam="${work_dir}/noncanonical-mate.source.bam"
  allowlist="${work_dir}/noncanonical-mate.canonical.txt"
  output_bam="${work_dir}/noncanonical-mate.filtered.bam"

  cat > "${source_sam}" <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:10000
@SQ	SN:chrUn_interspersed	LN:10000
@SQ	SN:chr2	LN:10000
@RG	ID:cells	SM:synthetic
canonical_with_noncanonical_mate	65	chr1	10	60	10M	chrUn_interspersed	30	0	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
noncanonical_mate	129	chrUn_interspersed	30	60	10M	chr1	10	0	TTACGTACGT	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
later_canonical	0	chr2	20	60	10M	*	0	0	ACGTACGTAA	IIIIIIIIII	RG:Z:cells	CB:Z:CELL1
EOF

  "${samtools_bin}" view --bam --output "${source_bam}" "${source_sam}"
  "${samtools_bin}" index "${source_bam}"
  printf 'chr1\nchr2\n' > "${allowlist}"

  SAMTOOLS_BIN="${samtools_bin}" \
    bash "${repo_root}/scripts/core_runtime/FilterCanonicalBam.sh" \
      "${source_bam}" "${output_bam}" "${allowlist}" 1 normal

  "${samtools_bin}" quickcheck "${output_bam}"
  "${samtools_bin}" view "${output_bam}" >/dev/null
  "${samtools_bin}" index --bai --output "${output_bam}.bai" "${output_bam}"
  [[ -s "${output_bam}.bai" ]]

  mapfile -t header_contigs < <("${samtools_bin}" view -H "${output_bam}" \
    | awk -F '\t' '$1 == "@SQ" { sub(/^SN:/, "", $2); print $2 }')
  [[ "${header_contigs[*]}" == "chr1 chrUn_interspersed chr2" ]]

  "${samtools_bin}" view "${source_bam}" chr1 chr2 \
    > "${work_dir}/noncanonical-mate.expected.sam"
  "${samtools_bin}" view "${output_bam}" \
    > "${work_dir}/noncanonical-mate.observed.sam"
  cmp "${work_dir}/noncanonical-mate.expected.sam" \
    "${work_dir}/noncanonical-mate.observed.sam"
}

run_case ucsc chr1 chr2 chrX chrY chrM
run_case ensembl 1 2 X Y MT
run_noncanonical_mate_case

printf 'canonical BAM and BigWig synthetic tests passed\n'

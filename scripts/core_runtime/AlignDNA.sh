#!/bin/bash

# Make sure hg19-blacklist.bed is in the assets dir (from https://github.com/nf-core/cutandrun/blob/master/assets/blacklists/hg19-blacklist.bed)

# ./AlignDNA.sh Human H3K27ac H2RD1 H2RD1_Human_H3K27ac_R1.fastq H2RD1_Human_H3K27ac_R2.fastq hg38-blacklist.bed SAM_RG_Header_H2RD1_Human_H3K27ac.tsv hg38_bwamem2_index 2913022398

set -e

threads="${ALIGN_DNA_THREADS:-8}"
view_threads="${ALIGN_DNA_VIEW_THREADS:-4}"
sort_threads="${ALIGN_DNA_SORT_THREADS:-8}"
sort_mem="${ALIGN_DNA_SORT_MEM:-1G}"

if (( view_threads > threads )); then
  view_threads="${threads}"
fi

if (( sort_threads > threads )); then
  sort_threads="${threads}"
fi

# Input Variables
modality="${1}"
sample_name="${2}"
R1="${3}"
R2="${4}"
blacklist_bed="${5}"
PathSam_header="${6}"
path_bwarefDB="${7}"
effsize="${8}"
outdir="${9}"
BWA_MEM2_BIN="${BWA_MEM2_BIN:-bwa-mem2}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

PathOutputBam=${outdir}/${sample_name}_${modality}.bam

RGID="${sample_name}_${modality}"

# Alignement
echo "Aligning"
echo "Using BWA_MEM2_BIN=${BWA_MEM2_BIN}"
echo "Using SAMTOOLS_BIN=${SAMTOOLS_BIN}"
"${BWA_MEM2_BIN}" mem -t ${threads} -C -o ${RGID}_TEMP.sam ${path_bwarefDB} ${R1} ${R2}

echo "Changing header and sorting alignment..."
# Prepare the header
"${SAMTOOLS_BIN}" view -H ${RGID}_TEMP.sam > ${RGID}_TEMPHEADER.sam
# Find the line number of the last occurrence of "@SQ"
last_sq_line=$(grep -n "@SQ" "${RGID}_TEMPHEADER.sam" | tail -n1 | cut -d':' -f1)
# Split the header based on the line number of the last "@SQ" occurrence
head -n "$last_sq_line" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER1.sam"
tail -n +"$((last_sq_line + 1))" "${RGID}_TEMPHEADER.sam" > "${RGID}_TEMPHEADER2.sam"

# Append the new header to the RG-replaced SAM file and convert to BAM
{ cat ${RGID}_TEMPHEADER1.sam ${PathSam_header} ${RGID}_TEMPHEADER2.sam; "${SAMTOOLS_BIN}" view --threads ${view_threads} ${RGID}_TEMP.sam; } | "${SAMTOOLS_BIN}" sort --threads ${sort_threads} -m ${sort_mem} -n -o ${RGID}_TEMP.bam -

echo "Removing blacklist-overlapping reads..."
# Remove reads overlapping with the blacklist
"${SAMTOOLS_BIN}" view --threads ${view_threads} --bam --with-header --output ${RGID}_TEMP_inBLRegions.bam --unoutput ${RGID}_TEMP_outBLRegions.bam -L ${blacklist_bed} ${RGID}_TEMP.bam

echo "Filtering properly paired mapped reads and sorting..."
"${SAMTOOLS_BIN}" view --threads ${view_threads} --with-header --require-flags 0x2 --output ${RGID}_TEMP_Good.bam ${RGID}_TEMP_outBLRegions.bam

## Sort and output final bam
"${SAMTOOLS_BIN}" sort -@ ${sort_threads} -m ${sort_mem} -o ${PathOutputBam} ${RGID}_TEMP_Good.bam

# Index the BAM
"${SAMTOOLS_BIN}" index --threads ${threads} --bai --output ${PathOutputBam}.bai ${PathOutputBam}

#Remove temp files
echo "Removing temp files..."
rm ${RGID}_TEMP.sam ${RGID}_TEMPHEADER.sam ${RGID}_TEMPHEADER1.sam ${RGID}_TEMPHEADER2.sam ${RGID}_TEMP.bam ${RGID}_TEMP_Good.bam ${RGID}_TEMP_inBLRegions.bam ${RGID}_TEMP_outBLRegions.bam

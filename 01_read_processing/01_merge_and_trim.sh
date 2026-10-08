#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-11              # 6 ATAC + 6 naked-DNA (gDNA) libraries
#SBATCH --mem-per-cpu 24000
#SBATCH -t 1-23:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 1: merge reads sequenced over several lanes (if any), then trim adaptors and run FastQC.
#
# Input : lane-level FASTQ files listed in ${PROJECT}/lane_files.tsv (see below)
# Output: ${RAW_READS}/<sample>_R{1,2}.fastq.merged.gz   merged reads
#         ${RESULTS}/1.trimmed/<sample>_R{1,2}.fastq.merged.gz_val_{1,2}.fq.gz   trimmed reads + FastQC reports
#
# lane_files.tsv: tab-delimited, one line per library. Column 1 is the sample name
# (e.g. 1aAc_3dpf_ATAC or 1aAc_3dpf_gDNA); remaining columns are the R1 lane files, then the matching R2 files
# separated by the literal column "R2" (any number of lanes):
#   1aAc_3dpf_ATAC   laneA_R1.fastq.gz   laneB_R1.fastq.gz   R2   laneA_R2.fastq.gz   laneB_R2.fastq.gz

source "$(dirname "$0")/../config.sh"
source trim_galore-0.5.0
ml java
source fastqc-0.11.9

SAMPLES=(${ATAC_SAMPLES} ${GDNA_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
TRIMDIR=${RESULTS}/1.trimmed
mkdir -p "${RAW_READS}" "${TRIMDIR}"

# 1a. merge lanes
line=$(awk -F'\t' -v s="${SAMPLE}" '$1==s' "${PROJECT}/lane_files.tsv")
r1=$(echo "${line}" | awk -F'\t' '{for(i=2;i<=NF && $i!="R2";i++) printf "%s ", $i}')
r2=$(echo "${line}" | awk -F'\t' '{f=0; for(i=2;i<=NF;i++){ if($i=="R2"){f=1; continue} if(f) printf "%s ", $i}}')
cat ${r1} > "${RAW_READS}/${SAMPLE}_R1.fastq.merged.gz"
cat ${r2} > "${RAW_READS}/${SAMPLE}_R2.fastq.merged.gz"

# 1b. trim adaptors (paired-end) and run FastQC on the trimmed reads
trim_galore --output_dir "${TRIMDIR}" --paired --fastqc \
  "${RAW_READS}/${SAMPLE}_R1.fastq.merged.gz" "${RAW_READS}/${SAMPLE}_R2.fastq.merged.gz"

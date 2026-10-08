#!/bin/bash -e
#SBATCH -p ei-largemem
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --array=0-11
#SBATCH --mem 512GB
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 2: build the bowtie2 index (once) and align trimmed reads to the A. calliptera genome (fAstCal1.2).
# Settings follow the ENCODE ATAC-seq pipeline: up to 4 alignments per read (-k 4), maximum fragment length 2 kb (-X2000).
#
# Run once before submitting the array:  bowtie2-build ${GENOME_FA} ${GENOME_ID}   (module: bowtie2/2.2.6)
#
# Input : ${RESULTS}/1.trimmed/*_val_{1,2}.fq.gz
# Output: ${RESULTS}/2.aligned/<sample>.bam (+ .align.log, _flagstat_qc1.txt)

source "$(dirname "$0")/../config.sh"
ml bowtie2/2.2.6
ml samtools/1.7

SAMPLES=(${ATAC_SAMPLES} ${GDNA_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
OUT=${RESULTS}/2.aligned
mkdir -p "${OUT}"

bowtie2 -k ${BOWTIE2_K} -X2000 --mm --threads ${THREADS} -x "${GENOME_ID}" \
  -1 "${RESULTS}/1.trimmed/${SAMPLE}_R1.fastq.merged.gz_val_1.fq.gz" \
  -2 "${RESULTS}/1.trimmed/${SAMPLE}_R2.fastq.merged.gz_val_2.fq.gz" \
  2> "${OUT}/${SAMPLE}.align.log" \
  | samtools view -Su /dev/stdin | samtools sort -o "${OUT}/${SAMPLE}.bam"

samtools flagstat "${OUT}/${SAMPLE}.bam" > "${OUT}/${SAMPLE}_flagstat_qc1.txt"

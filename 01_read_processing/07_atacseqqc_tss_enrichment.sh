#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH --array=0-5
#SBATCH --mem 64000
#SBATCH -t 0-12:00
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 7: run ATACseqQC (TSS enrichment, PT/NFR scores, Tn5-shifted BAM) for each ATAC library.
# Input : ${RESULTS}/4.filtered/<sample>.nochrM.nodup.filt.bam
# Output: ${RESULTS}/7.qc/<sample>/  (shifted BAM, score tables and plots)

source "$(dirname "$0")/../config.sh"
SAMPLES=(${ATAC_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
# GTF must be uncompressed for rtracklayer::import
GTF=${GTF_R101%.gz}
[ -f "${GTF}" ] || gunzip -c "${GTF_R101}" > "${GTF}"

Rscript --vanilla "$(dirname "$0")/07_atacseqqc_tss_enrichment.R" \
  "${RESULTS}/4.filtered/${SAMPLE}.nochrM.nodup.filt.bam" "${GTF}" "${RESULTS}/7.qc/${SAMPLE}" "${SAMPLE}"

# Coordinate-sort and index the Tn5-shifted BAM; these files are the input to differential accessibility (04_*) and
# TF footprinting (05_*). ATACseqQC writes the shifted BAM next to its plots.
ml samtools/1.3
mkdir -p "${RESULTS}/7.qc/shifted_bams"
samtools sort -o "${RESULTS}/7.qc/shifted_bams/${SAMPLE}.nochrM.nodup.filt.shifted.sorted.bam" \
  "${RESULTS}/7.qc/${SAMPLE}/${SAMPLE}.nochrM.nodup.filt.shifted.bam"
samtools index "${RESULTS}/7.qc/shifted_bams/${SAMPLE}.nochrM.nodup.filt.shifted.sorted.bam"

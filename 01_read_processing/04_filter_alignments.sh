#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH --array=0-5
#SBATCH --mem 120000
#SBATCH -t 0-10:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 4: post-alignment filtering (ENCODE ATAC-seq standard).
# Sort, mark duplicates (Sambamba), then keep properly paired (-f 2), primary, mapped, non-duplicate reads that pass
# platform QC with MAPQ >= 30 (-F 1804 -q 30). Fragment-length distributions are collected with Picard.
#
# Input : ${RESULTS}/3.chrM/<sample>.nochrM.bam
# Output: ${RESULTS}/4.filtered/<sample>.nochrM.nodup.filt.bam (+ .bai), .flagstat.qc, Picard insert-size metrics/histogram

source "$(dirname "$0")/../config.sh"
ml samtools/1.3
ml sambamba/0.6.5
ml python/3.5
ml zlib/1.2.8
ml glib/2.40
ml java/8.45
ml picard/1.140

SAMPLES=(${ATAC_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
IN=${RESULTS}/3.chrM/${SAMPLE}.nochrM.bam
OUT=${RESULTS}/4.filtered
mkdir -p "${OUT}/tmp"

SORTED=${OUT}/${SAMPLE}.nochrM.sorted.bam
DUP=${OUT}/${SAMPLE}.nochrM.sorted.dup.bam
FILT=${OUT}/${SAMPLE}.nochrM.nodup.filt.bam

sambamba sort -m 88G -t 1 --tmpdir "${OUT}/tmp" -o "${SORTED}" -u "${IN}"
sambamba markdup -l 0 -t 1 "${SORTED}" "${DUP}"
samtools view -F 1804 -f 2 -q 30 -b "${DUP}" > "${FILT}"
samtools index "${FILT}" "${FILT}.bai"
samtools flagstat "${FILT}" > "${OUT}/${SAMPLE}.nochrM.nodup.filt.flagstat.qc"

# fragment-length distribution
java -Xmx2g -jar "${PICARD_JAR:-CollectInsertSizeMetrics.jar}" \
  R="${GENOME_FA}" I="${FILT}" \
  O="${FILT}_PicardInsertMetrics.txt" H="${FILT}_insert_size_histogram.pdf" M=0.5

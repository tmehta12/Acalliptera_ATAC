#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-2                 # 3dpf, 7dpf, 12dpf
#SBATCH --mem-per-cpu 32000
#SBATCH -t 0-03:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# TF footprinting is run on the two replicates of each stage pooled, to increase footprint depth.
#   a. merge the Tn5-shifted replicate BAMs (${RESULTS}/7.qc/shifted_bams, from 01_read_processing/07)
#   b. sort, mark duplicates (Sambamba), keep properly paired MAPQ >= 30 reads, sort and index
#   c. union of the two replicate peak sets (concatenated, not merged, so per-peak statistics are retained)
#
# Output: ${RESULTS}/10.footprinting/Ac_<stage>_ATAC.nochrM.REPmerged.nodup.filt.sorted.bam (+ .bai)
#         ${RESULTS}/10.footprinting/Ac_<stage>_ATAC_peaks.final.REPcat.narrowPeak

source "$(dirname "$0")/../config.sh"
ml samtools/1.3
ml sambamba/0.6.5

STAGE_ARR=(${STAGES})
STAGE=${STAGE_ARR[${SLURM_ARRAY_TASK_ID}]}
case ${STAGE} in 3dpf) n=1;; 7dpf) n=2;; 12dpf) n=3;; esac
BAMS=${RESULTS}/7.qc/shifted_bams
OUT=${RESULTS}/10.footprinting
mkdir -p "${OUT}" && cd "${OUT}"

R1=${n}aAc_${STAGE}_ATAC; R2=${n}bAc_${STAGE}_ATAC
P=Ac_${STAGE}_ATAC.nochrM.REPmerged

samtools merge -f ${P}.bam ${BAMS}/${R1}.nochrM.nodup.filt.shifted.sorted.bam ${BAMS}/${R2}.nochrM.nodup.filt.shifted.sorted.bam
samtools sort ${P}.bam -o ${P}.sorted.bam
sambamba markdup -l 0 -t 2 ${P}.sorted.bam ${P}.sorted.dup.bam
samtools view -F 1804 -f 2 -q 30 -b ${P}.sorted.dup.bam > ${P}.nodup.filt.bam
samtools sort ${P}.nodup.filt.bam -o ${P}.nodup.filt.sorted.bam
samtools index ${P}.nodup.filt.sorted.bam ${P}.nodup.filt.sorted.bam.bai
samtools flagstat ${P}.nodup.filt.bam > ${P}.flagstat.qc
rm -f ${P}.bam ${P}.sorted.bam ${P}.sorted.dup.bam ${P}.nodup.filt.bam

cat ${RESULTS}/6.idr/${R1}_peaks.final.narrowPeak ${RESULTS}/6.idr/${R2}_peaks.final.narrowPeak \
  | sort -k1,1 -k2,2n | sed -e 's|        |\t|g' > Ac_${STAGE}_ATAC_peaks.final.REPcat.narrowPeak

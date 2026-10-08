#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-2                 # 3dpf, 7dpf, 12dpf
#SBATCH --mem-per-cpu 32000
#SBATCH -t 0-23:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# TF footprints and motif matching (HINT-ATAC, Regulatory Genomics Toolbox v0.13.0).
#   a. rgt-hint footprinting: footprints within the stage's union of replicate peaks (bias-corrected, paired-end ATAC mode)
#   b. rgt-motifanalysis matching (FPR 0.0001): motif matches within footprints, using species-specific and cichlid-wide
#      PWMs plus vertebrate PWMs (JASPAR, HOCOMOCO) - see rgt_data_config.md
#
# Input : ${RESULTS}/10.footprinting/Ac_<stage>_ATAC.nochrM.REPmerged.nodup.filt.sorted.bam and Ac_<stage>_ATAC_peaks.final.REPcat.narrowPeak
# Output: ${RESULTS}/10.footprinting/ac_fp/Ac_<stage>/Ac_<stage>.bed (footprints, column 5 = tag count)
#         ${RESULTS}/10.footprinting/ac_fp/Ac_<stage>/Ac_<stage>_mpbs.bed (motif-predicted binding sites: motif, bit-score, strand)

source "$(dirname "$0")/../config.sh"

STAGE_ARR=(${STAGES})
STAGE=${STAGE_ARR[${SLURM_ARRAY_TASK_ID}]}
FP=${RESULTS}/10.footprinting
OUTDIR=${FP}/ac_fp/Ac_${STAGE}
mkdir -p "${OUTDIR}"

BAM=${FP}/Ac_${STAGE}_ATAC.nochrM.REPmerged.nodup.filt.sorted.bam
PEAKS=${FP}/Ac_${STAGE}_ATAC_peaks.final.REPcat.narrowPeak

# a. footprints
rgt-hint footprinting --atac-seq --paired-end --organism=${GENOME_ID} \
  --output-location="${OUTDIR}" --output-prefix="Ac_${STAGE}" "${BAM}" "${PEAKS}"

# b. motif matching within footprints
rgt-motifanalysis matching --filter "database:cichlidacCSsp,cichlidCW,cichlidJASPAR,jaspar_vertebrates,hocomoco" \
  --organism ${GENOME_ID} --output-location "${OUTDIR}" --input-files "${OUTDIR}/Ac_${STAGE}.bed"

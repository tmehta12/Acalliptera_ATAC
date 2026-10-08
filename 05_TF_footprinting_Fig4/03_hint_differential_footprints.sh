#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --array=0-2                 # 3 vs 7, 3 vs 12, 7 vs 12 dpf
#SBATCH --mem-per-cpu 32000
#SBATCH -t 0-23:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Bias-corrected TF footprint profiles and differential footprint statistics between stages
# (rgt-hint differential, HINT-ATAC). Used for the TF footprint line plots.
#
# Input : ac_fp/Ac_<stage>/Ac_<stage>_mpbs.bed and the pooled BAMs (02_hint_footprinting_motif_matching.sh)
# Output: ${RESULTS}/10.footprinting/differential/<stage1>_vs_<stage2>/ (line plots, differential_statistics.txt)

source "$(dirname "$0")/../config.sh"

PAIRS=("3dpf 7dpf" "3dpf 12dpf" "7dpf 12dpf")
set -- ${PAIRS[${SLURM_ARRAY_TASK_ID}]}
A=$1; B=$2
FP=${RESULTS}/10.footprinting
OUTDIR=${FP}/differential/${A}_vs_${B}
mkdir -p "${OUTDIR}"

rgt-hint differential --organism ${GENOME_ID} --bc --nc 32 \
  --mpbs-files ${FP}/ac_fp/Ac_${A}/Ac_${A}_mpbs.bed,${FP}/ac_fp/Ac_${B}/Ac_${B}_mpbs.bed \
  --reads-files ${FP}/Ac_${A}_ATAC.nochrM.REPmerged.nodup.filt.sorted.bam,${FP}/Ac_${B}_ATAC.nochrM.REPmerged.nodup.filt.sorted.bam \
  --conditions ${A},${B} --output-location "${OUTDIR}"

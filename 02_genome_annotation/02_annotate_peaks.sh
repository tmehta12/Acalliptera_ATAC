#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 12000
#SBATCH -t 0-03:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Assign every peak to the genomic feature(s) it overlaps (annotation from 01_build_feature_annotation.sh).
# A peak overlapping several features (e.g. an aCNE inside a promoter) appears once per feature.
#
# Input : ${RESULTS}/6.idr/*_peaks.final.narrowPeak                       replicate peak sets (column 11 = IDR T/F)
#         ${RESULTS}/6.idr/Ac_<stage>_ATAC_peaks.final.allIDR.narrowPeak  both replicates per stage
# Output: ${RESULTS}/8.peak_features/<peakfile>.features.gff
#         columns: chrom, peak_id, feature, start, end, score, strand, frame, attributes where
#         attributes = signalValue;pValue;qValue;summit;IDR(T/F);gene_id;gene_name

source "$(dirname "$0")/../config.sh"
source bedtools-2.30.0

OUT=${RESULTS}/8.peak_features
mkdir -p "${OUT}"

for peaks in ${RESULTS}/6.idr/*_peaks.final.narrowPeak ${RESULTS}/6.idr/Ac_*_ATAC_peaks.final.allIDR.narrowPeak; do
  name=$(basename "${peaks}")
  awk '{gsub(/.*\//,"",$4);print}' OFS='\t' "${peaks}" | sort -V -k1,1 -k2,2n > "${OUT}/${name}"
  bedtools intersect -a "${ANNOT_BED}" -b "${OUT}/${name}" -wb \
    | awk '{print $11,$14,$4,$12,$13,$15,$6,".",$17";"$18";"$19";"$20";"$21";"$7";"$10}' OFS='\t' \
    > "${OUT}/${name}.features.gff"
done

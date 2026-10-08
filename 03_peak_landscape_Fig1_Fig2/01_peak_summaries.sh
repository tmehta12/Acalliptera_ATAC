#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 8000
#SBATCH -t 0-01:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Fig. 1B and 1D inputs: peak totals per stage, and the distance of each peak summit to the nearest TSS.
#
# Input : ${RESULTS}/6.idr/Ac_<stage>_ATAC_peaks.final.allIDR.narrowPeak, ${RESULTS}/6.idr/*_peaks.final.narrowPeak,
#         ${PROJECT}/annotation/Astatotilapia_calliptera.fAstCal1.2.100_gene.bed (made by 02_genome_annotation/01_build_feature_annotation.sh)
# Output: ${RESULTS}/9.peak_summaries/peak_counts.tsv            total peaks (IDR T and F) per stage
#         ${RESULTS}/9.peak_summaries/tss_distance_summary.tsv   per replicate: peaks with summit <=100 bp, 100 bp-5 kb, >5 kb from the nearest TSS

source "$(dirname "$0")/../config.sh"
source bedtools-2.30.0

OUT=${RESULTS}/9.peak_summaries
mkdir -p "${OUT}" && cd "${OUT}"
GENES=${PROJECT}/annotation/Astatotilapia_calliptera.fAstCal1.2.100_gene.bed   # chr_start, chr_end, chrom, start, end, gene, ., strand, gene_id, ...

# total peaks per stage
printf "stage\tpeaks\n" > peak_counts.tsv
for stage in ${STAGES}; do
  printf "%s\t%s\n" "${stage}" "$(wc -l < ${RESULTS}/6.idr/Ac_${stage}_ATAC_peaks.final.allIDR.narrowPeak | tr -d ' ')" >> peak_counts.tsv
done

# TSS positions (single-base BED): gene start on the + strand, gene end on the - strand
awk -F'\t' '{if($8=="+") tss=$4; else tss=$5; print $3"\t"tss"\t"tss+1"\t"$9"\t.\t"$8}' "${GENES}" | sort -k1,1 -k2,2n > tss.bed

# distance of each peak summit (start + narrowPeak column 10) to the nearest TSS
printf "stage\treplicate\tproximal_le100bp\tpromoter_distal_100bp_5kb\tfar_distal_gt5kb\n" > tss_distance_summary.tsv
for rep in ${ATAC_SAMPLES}; do
  stage=$(echo "${rep}" | awk -F'_' '{print $2}')
  awk -F'\t' '{summit=$2+$10; print $1"\t"summit"\t"summit+1"\t"$4}' "${RESULTS}/6.idr/${rep}_peaks.final.narrowPeak" \
    | sort -k1,1 -k2,2n > summit_${rep}.bed
  bedtools closest -a summit_${rep}.bed -b tss.bed -d -t first > closest_${rep}.bed
  prox=$(awk -F'\t' '$NF<=100' closest_${rep}.bed | wc -l | tr -d ' ')
  promdist=$(awk -F'\t' '$NF>100 && $NF<=5000' closest_${rep}.bed | wc -l | tr -d ' ')
  fardist=$(awk -F'\t' '$NF>5000' closest_${rep}.bed | wc -l | tr -d ' ')
  printf "%s\t%s\t%s\t%s\t%s\n" "${stage}" "${rep}" "${prox}" "${promdist}" "${fardist}" >> tss_distance_summary.tsv
  rm -f summit_${rep}.bed closest_${rep}.bed
done

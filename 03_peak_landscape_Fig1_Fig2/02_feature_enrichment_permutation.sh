#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 16000
#SBATCH -t 1-23:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Fig. 1C: are ATAC-seq peaks enriched in, or depleted from, each genomic feature?
# Peaks (merged, primary linkage groups only) are randomly re-positioned along the genome 200 times (bedtools shuffle,
# keeping each peak on its chromosome); the number of peaks overlapping each feature in the shuffled sets gives the
# expected count. Fold enrichment = observed / mean expected. See 02b_summarise_permutation_enrichment.R for p/q values.
#
# Input : ${RESULTS}/6.idr/Ac_<stage>_ATAC_peaks.final.allIDR.narrowPeak, ${ANNOT_BED}, ${CHROM_SIZES}
# Output: ${RESULTS}/9.peak_summaries/feature_enrichment/observed_summary.txt, null_<stage>_<feature>.txt

source "$(dirname "$0")/../config.sh"
source bedtools-2.30.0

NPERM=200
FEATS="five_prime_utr three_prime_utr exon intron intergenic 5kb_gene_promoter hCNE aCNE"
OUT=${RESULTS}/9.peak_summaries/feature_enrichment
mkdir -p "${OUT}" && cd "${OUT}"

# workspace: the numbered linkage groups (same sequences as the differential-accessibility analysis)
awk -F'\t' '$1 ~ /^[0-9]+$/' "${CHROM_SIZES}" > genome.txt

for feat in ${FEATS}; do
  awk -F'\t' -v f="${feat}" '$4==f' "${ANNOT_BED}" | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > annot_${feat}.bed
done

: > observed_summary.txt
for stage in ${STAGES}; do
  awk -F'\t' '$1 ~ /^[0-9]+$/' "${RESULTS}/6.idr/Ac_${stage}_ATAC_peaks.final.allIDR.narrowPeak" \
    | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -i - > segments_${stage}.bed

  for feat in ${FEATS}; do
    obs=$(bedtools intersect -u -a segments_${stage}.bed -b annot_${feat}.bed | wc -l | tr -d ' ')
    echo "${stage} ${feat} ${obs}" >> observed_summary.txt
    : > null_${stage}_${feat}.txt
  done

  for i in $(seq 1 ${NPERM}); do
    bedtools shuffle -i segments_${stage}.bed -g genome.txt -chrom -seed ${i} | sort -k1,1 -k2,2n > shuffled_${stage}.bed
    for feat in ${FEATS}; do
      bedtools intersect -u -a shuffled_${stage}.bed -b annot_${feat}.bed | wc -l | tr -d ' ' >> null_${stage}_${feat}.txt
    done
  done
  rm -f shuffled_${stage}.bed
done

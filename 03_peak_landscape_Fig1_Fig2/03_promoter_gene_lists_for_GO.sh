#!/bin/bash -e

# Fig. 2: genes whose 5 kb promoter overlaps an ATAC-seq peak at each stage (replicate peaks pooled, IDR T and F).
# The three lists are the query sets for the GO biological-process enrichment, which was run on the g:Profiler web server
# (g:GOst, https://biit.cs.ut.ee/gprofiler/gost):
#   organism: Astatotilapia calliptera, data source: GO biological process only, significance threshold: Benjamini-Hochberg FDR,
#   adjusted p < 0.05, statistical background: annotated genes (default). The GO:BP rows of the exported CSV files are
#   supplied in data/GO_BP_promoter_peaks_<stage>.tsv and are the input of 05_plot_Fig2.R.
#
# Input : ${RESULTS}/8.peak_features/Ac_<stage>_ATAC_peaks.final.allIDR.narrowPeak.features.gff  (02_genome_annotation/02_annotate_peaks.sh)
# Output: ${RESULTS}/9.peak_summaries/promoter_genes_<stage>.txt   one Ensembl gene ID per line

source "$(dirname "$0")/../config.sh"

OUT=${RESULTS}/9.peak_summaries
mkdir -p "${OUT}"
for stage in ${STAGES}; do
  # attribute field 6 of the features GFF is the Ensembl gene ID
  awk -F'\t' '$3=="5kb_gene_promoter"{split($9,a,";"); print a[6]}' \
    "${RESULTS}/8.peak_features/Ac_${stage}_ATAC_peaks.final.allIDR.narrowPeak.features.gff" | sort -u > "${OUT}/promoter_genes_${stage}.txt"
done

#!/usr/bin/env Rscript
# Merge the tiled 150 bp windows into differentially accessible regions (csaw::mergeWindows, tolerance 100 bp, max width
# 5 kb) and compute a region-level FDR with csaw::combineTests (Simes' method, BH-adjusted). A region is called when
# its FDR < 0.05; its direction is that of the window with the strongest evidence ('up' = gain, 'down' = loss,
# 'mixed' = windows change in both directions; mixed regions are not assigned to either direction downstream).
#
# Usage: Rscript 02_merge_windows_to_regions.R <outdir from 01_csaw_differential_windows.R> <promoter_bed>
# Outputs: DA_<comp>_sig_REGIONS.tsv, DA_<comp>_sig_REGIONS_promoter_genes.txt, DA_REGIONS_summary.tsv

suppressPackageStartupMessages({ library(csaw); library(GenomicRanges) })
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1]); promoter.bed <- args[2]

prom <- read.table(promoter.bed, header = FALSE, stringsAsFactors = FALSE, sep = "\t",
                   col.names = c("chrom", "start0", "end", "name", "score", "strand", "geneID", "c8", "c9", "symbol"))
prom.gr <- GRanges(prom$chrom, IRanges(prom$start0 + 1, prom$end), geneID = prom$geneID)

merge_one <- function(comp_name, tol = 100) {
  all.tab <- read.table(sprintf("DA_%s_all_windows.tsv", comp_name), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  all.gr  <- GRanges(all.tab$seqnames, IRanges(all.tab$start, all.tab$end))

  merged <- mergeWindows(all.gr, tol = tol, max.width = 5000)
  comb   <- combineTests(merged$ids, data.frame(logFC = all.tab$logFC, PValue = all.tab$PValue))
  comb$FDR <- p.adjust(comb$PValue, method = "BH")

  region.gr <- merged$regions
  sig.idx <- which(comb$FDR < 0.05)
  sig.gr  <- region.gr[sig.idx]
  sig.tab <- data.frame(as.data.frame(sig.gr)[, c("seqnames", "start", "end", "width")],
                        rep.logFC = comb$rep.logFC[sig.idx], direction = comb$direction[sig.idx],
                        PValue = comb$PValue[sig.idx], FDR = comb$FDR[sig.idx])

  ov <- findOverlaps(sig.gr, prom.gr, ignore.strand = TRUE)
  genes <- unique(prom.gr$geneID[subjectHits(ov)])
  write.table(sig.tab, sprintf("DA_%s_sig_REGIONS.tsv", comp_name), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(genes, sprintf("DA_%s_sig_REGIONS_promoter_genes.txt", comp_name))

  data.frame(comparison = comp_name, n_regions_total = length(region.gr), n_sig_regions_FDR05 = length(sig.idx),
             n_higher_late = sum(sig.tab$rep.logFC > 0), n_higher_early = sum(sig.tab$rep.logFC < 0),
             n_promoter_overlap = length(unique(queryHits(ov))), n_unique_genes = length(genes))
}

summary.df <- rbind(merge_one("3v7dpf"), merge_one("3v12dpf"), merge_one("7v12dpf"))
write.table(summary.df, "DA_REGIONS_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summary.df)

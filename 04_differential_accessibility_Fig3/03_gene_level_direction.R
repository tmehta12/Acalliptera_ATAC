#!/usr/bin/env Rscript
# Gene-level differential accessibility.
#   1. Genes gaining / losing promoter accessibility: promoters overlapping a significant region whose direction is
#      'up' (gain) or 'down' (loss) -> DA_<comp>_gain_genes.txt / DA_<comp>_loss_genes.txt
#   2. gene_DA_matrix.tsv: per gene and comparison, mean and strongest region log2FC, minimum region FDR and number of
#      regions in the promoter (input to the Fig. 3B heatmap)
#
# Usage: Rscript 03_gene_level_direction.R <outdir from 01/02> <promoter_bed>

suppressPackageStartupMessages(library(GenomicRanges))
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1]); promoter.bed <- args[2]
comps <- c("3v7dpf", "3v12dpf", "7v12dpf")

prom <- read.table(promoter.bed, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                   col.names = c("chrom", "start0", "end", "name", "score", "strand", "geneID", "c8", "c9", "symbol"))
prom.gr <- GRanges(prom$chrom, IRanges(prom$start0 + 1, prom$end), geneID = prom$geneID, symbol = prom$symbol)

# 1. direction-specific gene lists
for (comp in comps) {
  sig <- read.table(sprintf("DA_%s_sig_REGIONS.tsv", comp), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  for (dir in c("up", "down")) {
    sub <- sig[sig$direction == dir, ]
    ov  <- findOverlaps(GRanges(sub$seqnames, IRanges(sub$start, sub$end)), prom.gr, ignore.strand = TRUE)
    genes <- unique(prom.gr$geneID[subjectHits(ov)])
    writeLines(genes, sprintf("DA_%s_%s_genes.txt", comp, if (dir == "up") "gain" else "loss"))
  }
}

# 2. gene x comparison matrix
out <- list()
for (cn in comps) {
  sig <- read.table(sprintf("DA_%s_sig_REGIONS.tsv", cn), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ov  <- findOverlaps(GRanges(sig$seqnames, IRanges(sig$start, sig$end)), prom.gr, ignore.strand = TRUE)
  d   <- data.frame(geneID = prom.gr$geneID[subjectHits(ov)], logFC = sig$rep.logFC[queryHits(ov)], FDR = sig$FDR[queryHits(ov)])
  agg <- aggregate(cbind(logFC, FDR) ~ geneID, d, mean)
  agg <- merge(agg, setNames(aggregate(FDR ~ geneID, d, min), c("geneID", "minFDR")), by = "geneID")
  agg$n_regions <- as.integer(table(d$geneID)[agg$geneID])
  agg <- merge(agg, setNames(aggregate(logFC ~ geneID, d, function(x) x[which.max(abs(x))]), c("geneID", "maxLogFC")), by = "geneID")
  colnames(agg)[-1] <- paste0(colnames(agg)[-1], "_", cn)
  out[[cn]] <- agg
}
m <- Reduce(function(a, b) merge(a, b, by = "geneID", all = TRUE), out)
sym <- prom[!duplicated(prom$geneID), c("geneID", "symbol")]
m <- merge(sym, m, by = "geneID", all.y = TRUE)
write.table(m, "gene_DA_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

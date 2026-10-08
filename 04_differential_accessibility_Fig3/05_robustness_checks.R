#!/usr/bin/env Rscript
# Robustness of the accessibility-gain bias and replicate concordance (Supplementary Table S3).
#   A. for each comparison, the direction of change of the significant windows in every individual replicate
#      (late replicate vs mean of the early stage, and late mean vs each early replicate)
#   B. replicate-vs-replicate agreement within stage (median log2FC and correlation of log CPM across windows tested)
#   C. median log2FC of genes gaining / losing accessibility, and overlap of the 7 vs 12 dpf gain/loss genes with the
#      earlier comparisons
#
# Usage: Rscript 05_robustness_checks.R <outdir from 01-03> <promoter_bed>

suppressPackageStartupMessages({ library(csaw); library(edgeR); library(GenomicRanges) })
args <- commandArgs(trailingOnly = TRUE)
setwd(args[1]); promoter.bed <- args[2]

win  <- readRDS("win.data.all.rds"); bins <- readRDS("bins.all.rds")
colnames(win) <- c("3dpf_A", "3dpf_B", "7dpf_A", "7dpf_B", "12dpf_A", "12dpf_B")
nf   <- normFactors(bins, se.out = FALSE)                      # background-bin normalisation, all six libraries
lcpm <- log2(t(t(assay(win) + 0.5) / (win$totals * nf / 1e6)))

# A. per-replicate direction of the significant windows
for (cn in c("3v7dpf", "3v12dpf", "7v12dpf")) {
  sig <- read.table(sprintf("DA_%s_sig_windows.tsv", cn), header = TRUE, sep = "\t")
  idx <- queryHits(findOverlaps(rowRanges(win), GRanges(sig$seqnames, IRanges(sig$start, sig$end)), type = "equal"))
  L <- lcpm[idx, , drop = FALSE]
  st <- strsplit(cn, "v")[[1]]; early <- sub("dpf", "", st[1]); late <- sub("dpf", "", st[2])
  e <- paste0(early, "dpf_", c("A", "B")); l <- paste0(late, "dpf_", c("A", "B"))
  cat(sprintf("\n== %s: %d significant windows matched (of %d)\n", cn, length(idx), nrow(sig)))
  for (ee in e) { f <- rowMeans(L[, l, drop = FALSE]) - L[, ee]
    cat(sprintf("  late (mean) vs %-8s: median log2FC = %5.2f; %% windows higher in late = %.1f\n", ee, median(f), 100 * mean(f > 0))) }
  for (ll in l) { f <- L[, ll] - rowMeans(L[, e, drop = FALSE])
    cat(sprintf("  %-8s vs early (mean): median log2FC = %5.2f; %% windows higher in late = %.1f\n", ll, median(f), 100 * mean(f > 0))) }
}

# B. replicate agreement across all windows tested in 3 vs 7 dpf
tested <- read.table("DA_3v7dpf_all_windows.tsv", header = TRUE, sep = "\t")
idx <- queryHits(findOverlaps(rowRanges(win), GRanges(tested$seqnames, IRanges(tested$start, tested$end)), type = "equal"))
L <- lcpm[idx, ]
cat(sprintf("\nWindows tested in 3 vs 7 dpf (n = %d)\n", length(idx)))
cat(sprintf("  median log2FC A-B: 3 dpf %.2f; 7 dpf %.2f; 12 dpf %.2f\n", median(L[, "3dpf_A"] - L[, "3dpf_B"]),
            median(L[, "7dpf_A"] - L[, "7dpf_B"]), median(L[, "12dpf_A"] - L[, "12dpf_B"])))
cat(sprintf("  correlation (log CPM) A vs B: 3 dpf %.3f; 7 dpf %.3f; 12 dpf %.3f; 3 dpf A vs 7 dpf A %.3f\n",
            cor(L[, "3dpf_A"], L[, "3dpf_B"]), cor(L[, "7dpf_A"], L[, "7dpf_B"]), cor(L[, "12dpf_A"], L[, "12dpf_B"]), cor(L[, "3dpf_A"], L[, "7dpf_A"])))

# C. gene-level effect sizes and overlap between comparisons
prom <- read.table(promoter.bed, sep = "\t", stringsAsFactors = FALSE, quote = "",
                   col.names = c("chrom", "start0", "end", "name", "score", "strand", "geneID", "c8", "c9", "symbol"))
pg <- GRanges(prom$chrom, IRanges(prom$start0 + 1, prom$end), geneID = prom$geneID)
gene_lfc <- function(cn, dir) {
  s <- read.table(sprintf("DA_%s_sig_REGIONS.tsv", cn), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  s <- s[s$direction == dir, ]
  ov <- findOverlaps(GRanges(s$seqnames, IRanges(s$start, s$end)), pg, ignore.strand = TRUE)
  tapply(s$rep.logFC[queryHits(ov)], pg$geneID[subjectHits(ov)], mean)
}
cat("\n")
for (cn in c("3v7dpf", "3v12dpf", "7v12dpf")) {
  up <- gene_lfc(cn, "up"); dn <- gene_lfc(cn, "down")
  cat(sprintf("%s: gain genes n = %d, median log2FC %.2f (IQR %.2f-%.2f) | loss genes n = %d, median log2FC %.2f (IQR %.2f to %.2f)\n",
              cn, length(up), median(up), quantile(up, .25), quantile(up, .75), length(dn), median(dn), quantile(dn, .25), quantile(dn, .75)))
}
gl <- function(cn, lab) readLines(sprintf("DA_%s_%s_genes.txt", cn, lab))
g7g <- gl("7v12dpf", "gain"); g7l <- gl("7v12dpf", "loss")
g37g <- gl("3v7dpf", "gain"); g37l <- gl("3v7dpf", "loss"); g312g <- gl("3v12dpf", "gain"); g312l <- gl("3v12dpf", "loss")
cat(sprintf("\n7 vs 12 dpf gain genes (n = %d): also gain 3 vs 7 = %d; also gain 3 vs 12 = %d (%.0f%%); not gaining in 3 vs 7 = %d\n",
            length(g7g), sum(g7g %in% g37g), sum(g7g %in% g312g), 100 * mean(g7g %in% g312g), sum(!g7g %in% g37g)))
cat(sprintf("7 vs 12 dpf loss genes (n = %d): gain 3 vs 7 = %d (%.0f%%); gain 3 vs 12 = %d; lose 3 vs 12 = %d; lose 3 vs 7 = %d\n",
            length(g7l), sum(g7l %in% g37g), 100 * mean(g7l %in% g37g), sum(g7l %in% g312g), sum(g7l %in% g312l), sum(g7l %in% g37l)))
cat(sprintf("genes gaining in 3 vs 7 and losing in 7 vs 12 but not significant in 3 vs 12: %d\n", sum(g7l %in% g37g & !(g7l %in% c(g312g, g312l)))))

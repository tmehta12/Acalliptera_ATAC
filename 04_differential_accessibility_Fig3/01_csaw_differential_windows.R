#!/usr/bin/env Rscript
# Differential chromatin accessibility between stages with csaw + edgeR (Fig. 3).
#
#   - reads are counted in 150 bp windows (50 bp spacing) across the numbered linkage groups (LG1-LG23; unplaced
#     scaffolds are excluded) from the Tn5-shifted, filtered BAM files of the six libraries
#   - windows are retained if their signal is > 3-fold above a background estimate from 2 kb bins
#     (csaw::filterWindowsGlobal), and libraries are normalised for composition bias with the 2 kb background bins
#   - pairwise edgeR quasi-likelihood F-tests (3 vs 7, 3 vs 12, 7 vs 12 dpf), Benjamini-Hochberg FDR < 0.05;
#     a positive log2FC means higher accessibility at the later stage ('gain')
#   - significant windows are overlapped with 5 kb gene promoters
#
# Usage:
#   Rscript 01_csaw_differential_windows.R <bam_dir> <promoter_bed> <chrom_sizes> <outdir>
#     bam_dir       directory with <sample>.nochrM.nodup.filt.shifted.sorted.bam (${RESULTS}/7.qc/shifted_bams from 01_read_processing/07)
#     promoter_bed  10-column 5 kb promoter BED (the '5kb_gene_promoter' rows of the feature annotation)
#     chrom_sizes   two-column chromosome sizes file
# Outputs (in outdir): DA_<comp>_all_windows.tsv, DA_<comp>_sig_windows.tsv, DA_<comp>_sig_promoter_genes.txt,
#   DA_summary.tsv, bam_read_counts.tsv, csaw_library_totals.tsv, win.data.all.rds, bins.all.rds
# Tested with csaw 1.42.0, edgeR 4.6.3, R 4.5.

suppressPackageStartupMessages({ library(csaw); library(edgeR); library(rtracklayer); library(GenomicRanges); library(Rsamtools) })

args <- commandArgs(trailingOnly = TRUE)
bamdir <- args[1]; promoter.bed <- args[2]; chrom.sizes.file <- args[3]; outdir <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE); setwd(outdir)
set.seed(1)

bam.files <- c(
  "3dpf_A"  = file.path(bamdir, "1aAc_3dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"),
  "3dpf_B"  = file.path(bamdir, "1bAc_3dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"),
  "7dpf_A"  = file.path(bamdir, "2aAc_7dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"),
  "7dpf_B"  = file.path(bamdir, "2bAc_7dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"),
  "12dpf_A" = file.path(bamdir, "3aAc_12dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"),
  "12dpf_B" = file.path(bamdir, "3bAc_12dpf_ATAC.nochrM.nodup.filt.shifted.sorted.bam"))
stopifnot(all(file.exists(bam.files)))
stage <- setNames(c("3dpf", "3dpf", "7dpf", "7dpf", "12dpf", "12dpf"), names(bam.files))

# numbered linkage groups only
cs <- read.table(chrom.sizes.file, header = FALSE, stringsAsFactors = FALSE, col.names = c("chrom", "length"))
primary.chroms <- cs$chrom[grepl("^[0-9]+$", cs$chrom)]

# mapped reads per library
read.counts <- do.call(rbind, lapply(names(bam.files), function(nm) {
  idx <- idxstatsBam(bam.files[nm])
  data.frame(sample = nm, stage = stage[nm], mapped_reads_idx = sum(idx$mapped), unmapped_reads_idx = sum(idx$unmapped),
             approx_pairs = round(sum(idx$mapped) / 2))
}))
write.table(read.counts, "bam_read_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# window counts (150 bp, 50 bp spacing) and 2 kb background bins
param <- readParam(pe = "both", max.frag = 2000, restrict = primary.chroms, dedup = FALSE, minq = NA)
win.data.all <- windowCounts(bam.files, width = 150, spacing = 50, ext = NA, param = param, filter = 10)
colnames(win.data.all) <- names(bam.files)
write.table(data.frame(sample = names(bam.files), stage = stage, csaw_library_total = win.data.all$totals),
            "csaw_library_totals.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
bins.all <- windowCounts(bam.files, bin = TRUE, width = 2000, param = param)
colnames(bins.all) <- names(bam.files)
saveRDS(win.data.all, "win.data.all.rds"); saveRDS(bins.all, "bins.all.rds")

# 5 kb promoters
prom <- read.table(promoter.bed, header = FALSE, stringsAsFactors = FALSE, sep = "\t",
                   col.names = c("chrom", "start0", "end", "name", "score", "strand", "geneID", "c8", "c9", "symbol"))
prom.gr <- GRanges(prom$chrom, IRanges(prom$start0 + 1, prom$end), strand = ifelse(prom$strand %in% c("+", "-"), prom$strand, "*"),
                   geneID = prom$geneID, symbol = prom$symbol)

run_comparison <- function(early_label, late_label, comp_name) {
  samp.idx <- which(stage %in% c(early_label, late_label))
  samp.idx <- samp.idx[order(match(stage[samp.idx], c(early_label, late_label)))]
  samp.names <- names(stage)[samp.idx]

  win.sub  <- win.data.all[, samp.idx]
  bins.sub <- bins.all[, samp.idx]

  filter.stat <- filterWindowsGlobal(win.sub, bins.sub)
  filtered <- win.sub[filter.stat$filter > log2(3), ]
  filtered <- normFactors(bins.sub, se.out = filtered)

  group  <- factor(stage[samp.names], levels = c(early_label, late_label))
  design <- model.matrix(~group); colnames(design) <- c("Intercept", "late_vs_early")
  y   <- estimateDisp(asDGEList(filtered), design, robust = TRUE)
  res <- glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = "late_vs_early")

  tab <- res$table; tab$FDR <- p.adjust(tab$PValue, method = "BH")
  full <- cbind(as.data.frame(rowRanges(filtered))[, c("seqnames", "start", "end")],
                logFC = tab$logFC, logCPM = tab$logCPM, Fstat = tab$F, PValue = tab$PValue, FDR = tab$FDR)
  sig <- full[full$FDR < 0.05, ]

  sig.gr <- GRanges(sig$seqnames, IRanges(sig$start, sig$end))
  ov <- findOverlaps(sig.gr, prom.gr, ignore.strand = TRUE)
  overlap.genes <- unique(prom.gr$geneID[subjectHits(ov)])
  sig$prom_overlap <- seq_len(nrow(sig)) %in% unique(queryHits(ov))
  sig$direction <- ifelse(sig$logFC > 0, paste0("higher_", late_label), paste0("higher_", early_label))

  write.table(full, sprintf("DA_%s_all_windows.tsv", comp_name), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(sig,  sprintf("DA_%s_sig_windows.tsv", comp_name), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(overlap.genes, sprintf("DA_%s_sig_promoter_genes.txt", comp_name))

  data.frame(comparison = comp_name, n_windows_tested = nrow(full), n_sig_DA_FDR05 = nrow(sig),
             n_higher_late = sum(sig$logFC > 0), n_higher_early = sum(sig$logFC < 0),
             n_promoter_overlap = length(unique(queryHits(ov))), n_unique_genes = length(overlap.genes))
}

summary.df <- rbind(run_comparison("3dpf", "7dpf", "3v7dpf"),
                    run_comparison("3dpf", "12dpf", "3v12dpf"),
                    run_comparison("7dpf", "12dpf", "7v12dpf"))
write.table(summary.df, "DA_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(summary.df)

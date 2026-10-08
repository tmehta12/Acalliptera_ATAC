# Fig. 1: ATAC-seq peak landscape of A. calliptera embryos at 3, 7 and 12 dpf.
#   A  experimental overview (schematic, drawn separately; supplied as a PNG)
#   B  total accessible-chromatin peaks per stage
#   C  enrichment of peaks in genomic features (permutation test, 02_feature_enrichment_permutation.sh)
#   D  distance of peak summits to the nearest TSS (mean +/- SD over the two replicates)
#
# Usage: Rscript 04_plot_Fig1.R <9.peak_summaries dir> <schematic.png> <output dir>

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(cowplot); library(png) })

args <- commandArgs(trailingOnly = TRUE)
summ_dir <- args[1]; schematic <- args[2]; outdir <- args[3]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

OK <- "#009E73"; MID <- "#E69F00"; VM <- "#D55E00"
stage_cols <- c("3 dpf" = OK, "7 dpf" = MID, "12 dpf" = VM)
stage_lab  <- function(x) factor(x, levels = c("3dpf", "7dpf", "12dpf"), labels = c("3 dpf", "7 dpf", "12 dpf"))

common_theme <- theme(
  axis.ticks.length = unit(-0.1, "cm"),
  axis.text   = element_text(size = 15, color = "black"),
  axis.title  = element_text(size = 19, face = "bold"),
  panel.grid.minor = element_blank(),
  strip.text  = element_text(face = "italic", size = 16),
  legend.text = element_text(size = 14), legend.title = element_text(size = 15))

# A. schematic
p_a <- ggdraw() + draw_image(readPNG(schematic))

# B. peaks per stage
peaks <- read.table(file.path(summ_dir, "peak_counts.tsv"), header = TRUE)
peaks$stage_lab <- stage_lab(peaks$stage)
p_b <- ggplot(peaks, aes(stage_lab, peaks, fill = stage_lab)) +
  geom_col() + scale_fill_manual(values = stage_cols, guide = "none") +
  scale_y_continuous(labels = scales::comma) +
  xlab("Developmental stage") + ylab("Total accessible\nchromatin peaks") + theme_bw() + common_theme

# C. feature enrichment
enr <- read.table(file.path(summ_dir, "feature_enrichment", "feature_enrichment_results.tsv"), header = TRUE, stringsAsFactors = FALSE)
feat_labels <- c(five_prime_utr = "5' UTR", three_prime_utr = "3' UTR", exon = "Exon", intron = "Intron", intergenic = "Intergenic",
                 `5kb_gene_promoter` = "5kb gene promoter", hCNE = "hCNE", aCNE = "aCNE")
enr$feature_lab <- factor(feat_labels[enr$feature], levels = c("aCNE", "hCNE", "Intergenic", "5kb gene promoter", "5' UTR", "Exon", "Intron", "3' UTR"))
enr$stage_lab   <- stage_lab(enr$stage)
p_c <- ggplot(enr, aes(feature_lab, fold_enrichment)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = observed, colour = qvalue)) +
  scale_color_gradient(low = OK, high = VM, limits = c(0, NA)) +
  coord_flip() + facet_grid(~stage_lab) + theme_bw() +
  xlab("Annotated region") + ylab("Fold enrichment") +
  labs(color = "q-value\n(FDR<0.05)", size = "Number of\npeaks overlapping\nannotation") +
  common_theme +
  theme(axis.text.x = element_text(margin = margin(5, 5, 0, 5, "pt")),
        axis.text.y = element_text(margin = margin(5, 5, 5, 5, "pt")), legend.title.align = 0.5)

# D. distance to TSS
tss <- read.table(file.path(summ_dir, "tss_distance_summary.tsv"), header = TRUE)
long <- tss %>% pivot_longer(c(proximal_le100bp, promoter_distal_100bp_5kb, far_distal_gt5kb), names_to = "category", values_to = "n")
long$category <- factor(long$category, levels = c("proximal_le100bp", "promoter_distal_100bp_5kb", "far_distal_gt5kb"),
                        labels = c("Proximal (≤100 bp)", "Promoter-distal (100 bp–5 kb)", "Far-distal (>5 kb)"))
d <- long %>% group_by(stage, category) %>% summarise(mean_n = mean(n), sd_n = sd(n), .groups = "drop")
d$stage_lab <- stage_lab(d$stage)
p_d <- ggplot(d, aes(category, mean_n, fill = stage_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = stage_cols) + theme_bw() +
  xlab("Peak summit distance to nearest TSS") + ylab("Mean number of peaks\n(± SD across replicates)") + labs(fill = "Stage") +
  common_theme + theme(axis.text.x = element_text(angle = 15, hjust = 1))

fig1 <- (p_a + p_b + plot_layout(widths = c(1.4, 1))) / p_c / p_d +
  plot_layout(heights = c(1.1, 1, 1)) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 26, face = "bold"))

ggsave(file.path(outdir, "Fig1_multipanel.png"), fig1, units = "in", width = 13, height = 15, dpi = 300, bg = "white")

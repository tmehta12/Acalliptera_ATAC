# Fig. 2: GO biological processes enriched among genes with a peak in their 5 kb promoter at each stage.
#   A  the ten most significant GO:BP terms enriched at all three stages
#   B  stage-specific GO:BP terms (enriched at one stage only; top 15 per stage by adjusted p-value)
#   C  number of stage-unique GO:BP terms per stage
#
# Usage: Rscript 05_plot_Fig2.R <data dir with GO_BP_promoter_peaks_<stage>.tsv> <output dir>
#   GO_BP_promoter_peaks_<stage>.tsv: GO:BP rows of the g:Profiler (g:GOst) export for the genes with a promoter peak at
#   that stage (gene lists from 03_promoter_gene_lists_for_GO.sh); columns term_id, term_name, adjusted_p_value,
#   term_size, query_size, intersection_size, effective_domain_size.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork); library(scales) })

args <- commandArgs(trailingOnly = TRUE)
godir <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

OK <- "#009E73"; VM <- "#D55E00"; MID <- "#E69F00"
stages <- c("3dpf", "7dpf", "12dpf"); slab <- c("3 dpf", "7 dpf", "12 dpf")

bp <- bind_rows(lapply(stages, function(s) {
  d <- read.table(file.path(godir, sprintf("GO_BP_promoter_peaks_%s.tsv", s)), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
  d$stage <- s; d
}))
bp$fold_enrichment <- (bp$intersection_size / bp$query_size) / (bp$term_size / bp$effective_domain_size)
bp$stage <- factor(bp$stage, levels = stages, labels = slab)
bp$nlp <- -log10(bp$adjusted_p_value)
bp$term_lab <- gsub("behavior", "behaviour", bp$term_name)

AXT <- 18; AXL <- 22; STR <- 20; LGT <- 17; LGL <- 18
common <- theme_bw() + theme(
  axis.text = element_text(size = AXT, colour = "black"),
  axis.title = element_text(size = AXL, face = "bold"),
  strip.text = element_text(size = STR, face = "italic"), strip.text.y = element_text(size = STR, face = "italic", angle = 0),
  legend.text = element_text(size = LGT), legend.title = element_text(size = LGL),
  panel.grid.minor = element_blank(), axis.ticks.length = unit(-0.1, "cm"),
  legend.key.height = unit(0.9, "cm"))
wrapl <- label_wrap(46)

# ---------- A: top 10 terms shared by all three stages (smallest adj. p across stages)
nst <- bp %>% group_by(term_name) %>% summarise(n = n_distinct(stage), minp = min(adjusted_p_value), .groups = "drop")
top_shared <- nst %>% filter(n == 3) %>% arrange(minp) %>% slice_head(n = 10) %>% pull(term_name)
dA <- bp %>% filter(term_name %in% top_shared)
dA$term_lab <- factor(dA$term_lab, levels = rev(gsub("behavior", "behaviour", top_shared)))
pA <- ggplot(dA, aes(x = fold_enrichment, y = term_lab)) +
  geom_point(aes(size = intersection_size, colour = nlp)) +
  facet_grid(. ~ stage) +
  scale_colour_gradient(low = OK, high = VM, name = expression(-log[10]~"(adj. P)")) +
  scale_size_continuous(range = c(5, 14), name = "Number\nof genes", breaks = pretty_breaks(3)) +
  scale_y_discrete(labels = wrapl) +
  scale_x_continuous(breaks = c(1.07, 1.13), limits = c(1.05, 1.16), expand = expansion(mult = 0.02)) +
  guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
  xlab("Fold enrichment") + ylab("GO Biological Process\n(shared across all stages)") + common

# ---------- B: stage-unique terms (top 15 per stage by adj. p)
nst_all <- bp %>% group_by(term_name) %>% summarise(n = n_distinct(stage), .groups = "drop")
uniq <- bp %>% inner_join(nst_all %>% filter(n == 1), by = "term_name")
dB <- uniq %>% group_by(stage) %>% arrange(adjusted_p_value, .by_group = TRUE) %>% slice_head(n = 15) %>% ungroup()
dB <- dB %>% arrange(stage, adjusted_p_value)
dB$term_lab <- factor(dB$term_lab, levels = rev(unique(dB$term_lab)))
pB <- ggplot(dB, aes(x = fold_enrichment, y = term_lab)) +
  geom_point(aes(size = intersection_size, colour = nlp)) +
  facet_grid(stage ~ ., scales = "free_y", space = "free_y") +
  scale_colour_gradient(low = OK, high = VM, name = expression(-log[10]~"(adj. P)")) +
  scale_size_continuous(range = c(4, 12), name = "Number\nof genes", breaks = pretty_breaks(3)) +
  scale_y_discrete(labels = wrapl) +
  scale_x_continuous(breaks = c(1.05, 1.10, 1.15), expand = expansion(mult = 0.08)) +
  guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) +
  xlab("Fold enrichment") + ylab("GO Biological Process\n(unique to one stage)") + common

# ---------- C: number of stage-unique terms
cnt <- uniq %>% count(stage)
pC <- ggplot(cnt, aes(x = stage, y = n, fill = stage)) + geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 8) +
  scale_fill_manual(values = c("3 dpf" = OK, "7 dpf" = MID, "12 dpf" = VM), guide = "none") +
  scale_y_continuous(limits = c(0, 88), expand = c(0, 0)) +
  xlab("Developmental stage") + ylab("Number of\nstage-unique terms") + common

ggsave(file.path(outdir, "Fig2a.GO_shared_dotplot.png"), pA, width = 14, height = 5.4, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig2b.GO_specific_dotplot.png"), pB, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(file.path(outdir, "Fig2c.unique_term_counts.png"), pC, width = 6.5, height = 4.6, dpi = 300, bg = "white")

fig2 <- pA / pB / wrap_plots(pC, plot_spacer(), widths = c(1, 1.1)) +
  plot_layout(heights = c(5.4, 12, 4.6)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 28, face = "bold"))
ggsave(file.path(outdir, "Fig2_multipanel.png"), fig2, width = 14, height = 22.8, dpi = 300, bg = "white")
cat("Shared top10:\n"); print(top_shared); cat("\nUnique counts:\n"); print(cnt)
cat("\nFE range panel A:", range(dA$fold_enrichment), " B:", range(dB$fold_enrichment), "\n")

#!/usr/bin/env Rscript
# Fig. 3: differential promoter accessibility between 3, 7 and 12 dpf.
#   A  volcano plots of the 150 bp windows for the three comparisons (log2FC = later vs earlier stage); labels give the
#      numbers of significant merged regions gaining / losing accessibility
#   B  log2FC of genes representing the four temporal patterns of promoter opening across the three comparisons
#      (blocks; genes ranked by region-level FDR, named genes always shown; grey = no significant change)
#   C  GO:BP enrichment of genes gaining promoter accessibility
#   D  GO:BP enrichment of genes losing promoter accessibility
#
# Usage: Rscript 06_plot_Fig3.R <workdir with outputs of 01-04> <output dir>
#   workdir needs: DA_<comp>_all_windows.tsv, DA_REGIONS_summary.tsv, gene_DA_matrix.tsv, GO_DA_<comp>_<gain|loss>.tsv
# Environment variable NPB sets the number of top-ranked genes per heatmap block (6 in the paper).

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ComplexHeatmap); library(circlize); library(magick); library(grid) })

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1]); OUT <- args[2]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

OK <- "#009E73"; VM <- "#D55E00"
comps <- c("3v7dpf", "3v12dpf", "7v12dpf"); comp_labels <- c("3 vs 7 dpf", "3 vs 12 dpf", "7 vs 12 dpf")
N_PER_BLOCK <- as.integer(Sys.getenv("NPB", "6"))

############################################################
# Panel A: volcano plots
############################################################
# numbers of significant merged regions per comparison (loss, gain), from DA_REGIONS_summary.tsv
da <- read.table("DA_REGIONS_summary.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sig_counts <- lapply(comps, function(cn) c(loss = da$n_higher_early[da$comparison == cn], gain = da$n_higher_late[da$comparison == cn]))
names(sig_counts) <- comps

vl <- list()
for (i in seq_along(comps)) {
  df <- read.table(sprintf("DA_%s_all_windows.tsv", comps[i]), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df$sig <- "NS"; df$sig[df$FDR < 0.05 & df$logFC > 0] <- "Gain"; df$sig[df$FDR < 0.05 & df$logFC < 0] <- "Loss"
  df$sig <- factor(df$sig, levels = c("Loss", "NS", "Gain"))
  set.seed(1)                                                       # non-significant windows are down-sampled for plotting
  ns <- which(df$sig == "NS"); keep <- if (length(ns) > 30000) sample(ns, 30000) else ns
  d <- df[c(which(df$sig != "NS"), keep), ]; d$comparison <- comp_labels[i]; vl[[i]] <- d
}
vdf <- bind_rows(vl); vdf$comparison <- factor(vdf$comparison, levels = comp_labels)
cl <- data.frame(comparison = factor(comp_labels, levels = comp_labels),
                 label = sapply(comps, function(cn) sprintf("Loss: %s\nGain: %s", format(sig_counts[[cn]]["loss"], big.mark = ","), format(sig_counts[[cn]]["gain"], big.mark = ","))),
                 x = -7.5, y = 27)
pv <- ggplot(vdf, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Loss" = OK, "NS" = "grey80", "Gain" = VM), breaks = c("Loss", "Gain")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text(data = cl, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 7.2, fontface = "bold", lineheight = 0.9) +
  theme_bw() + facet_grid(~comparison) +
  xlab(expression(log[2]~"Fold Change (late vs early)")) + ylab(expression(-log[10]~"(p-value)")) +
  labs(color = "Direction\n(FDR<0.05)") +
  theme(axis.ticks.length = unit(-0.1, "cm"),
        axis.text.x = element_text(margin = margin(5, 5, 0, 5, "pt"), size = 24), axis.text.y = element_text(margin = margin(5, 5, 5, 5, "pt"), size = 24),
        axis.title = element_text(size = 36, face = "bold"), axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "italic", size = 26),
        legend.title.align = 0.5, legend.text = element_text(size = 23), legend.title = element_text(size = 24))
ggsave(file.path(OUT, "Fig3a.volcano.png"), pv, units = "in", width = 19, height = 6.5, dpi = 300)

############################################################
# Panel B: heatmap of temporal patterns of promoter opening
############################################################
m <- read.table("gene_DA_matrix.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
a <- "logFC_3v7dpf"; b <- "logFC_3v12dpf"; cc <- "logFC_7v12dpf"
m$has_symbol <- !(is.na(m$symbol) | m$symbol %in% c("ensembl", ".") | grepl(":", m$symbol))

blocks <- list(
  "Opens 3→7 and 3→12 dpf,\nstable 7→12 dpf" = m[!is.na(m[[a]]) & m[[a]] > 0 & !is.na(m[[b]]) & m[[b]] > 0 & is.na(m[[cc]]), ],
  "Opens 3→7 dpf,\ncontinues to open\n7→12 dpf"       = m[!is.na(m[[a]]) & m[[a]] > 0 & !is.na(m[[cc]]) & m[[cc]] > 0, ],
  "Opens 3→7 dpf,\ncloses 7→12 dpf"                = m[!is.na(m[[a]]) & m[[a]] > 0 & !is.na(m[[cc]]) & m[[cc]] < 0, ],
  "Opens only\nafter 7 dpf"                                = m[is.na(m[[a]]) & !is.na(m[[cc]]) & m[[cc]] > 0, ])
rank_fdr <- function(df, i) if (i == 1) pmax(df$minFDR_3v7dpf, df$minFDR_3v12dpf) else df$minFDR_7v12dpf
named <- list(NULL, c("pawr", "ddit3"), c("neurod4", "irx7", "notch1b", "isl2b"), c("myo5b", "scinlb"))   # genes discussed in the text

sel <- list()
for (i in seq_along(blocks)) {
  df <- blocks[[i]][blocks[[i]]$has_symbol, ]
  df <- df[order(rank_fdr(df, i)), ]
  df <- df[!duplicated(paste(round(df[[a]], 3), round(df[[b]], 3), round(df[[cc]], 3))), ]   # drop genes with identical profiles (shared promoters)
  top <- head(df, N_PER_BLOCK)
  extra <- if (!is.null(named[[i]])) blocks[[i]][blocks[[i]]$symbol %in% named[[i]] & !(blocks[[i]]$symbol %in% top$symbol), ] else df[0, ]
  top$block <- names(blocks)[i]; top$named <- FALSE
  if (nrow(extra) > 0) { extra$block <- names(blocks)[i]; extra$named <- TRUE }
  sel[[i]] <- rbind(top, extra)
}
sel <- do.call(rbind, sel)
sel$block <- factor(sel$block, levels = names(blocks))
sel_out <- sel; sel_out$block <- gsub("\n", " ", as.character(sel_out$block), fixed = TRUE)
write.table(sel_out[, c("block", "geneID", "symbol", a, b, cc, "minFDR_3v7dpf", "minFDR_3v12dpf", "minFDR_7v12dpf", "named")],
            file.path(OUT, "Fig3B_gene_matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mat <- as.matrix(sel[, c(a, b, cc)]); rownames(mat) <- sel$symbol
colnames(mat) <- comp_labels
ht <- Heatmap(mat, name = "log2FC", na_col = "grey85", col = colorRamp2(c(-6, 0, 6), c(OK, "white", VM)),
              cluster_rows = FALSE, cluster_columns = FALSE,
              row_split = sel$block, cluster_row_slices = FALSE, row_gap = unit(5, "mm"),
              row_title_rot = 0, row_title_gp = gpar(fontsize = 17, fontface = "bold"),
              show_row_names = TRUE, row_names_side = "right", row_names_gp = gpar(fontsize = 19, fontface = "italic"),
              column_names_side = "top", column_names_rot = 45, column_names_gp = gpar(fontsize = 19),
              rect_gp = gpar(col = "white", lwd = 1.5),
              heatmap_legend_param = list(title_gp = gpar(fontsize = 18, fontface = "bold"), labels_gp = gpar(fontsize = 16),
                                          legend_height = unit(4, "cm"), at = c(-6, -3, 0, 3, 6), labels = c("≤-6", "-3", "0", "3", "≥6")),
              width = unit(7.5, "cm"), height = unit(nrow(sel) * 0.78, "cm"))
png(file.path(OUT, "Fig3b.heatmap.png"), width = 2900, height = 3900, res = 300)
draw(ht, heatmap_legend_side = "right", padding = unit(c(3, 3, 3, 3), "mm"))
invisible(dev.off())

############################################################
# Panels C and D: GO enrichment of gaining / losing genes (top 10 terms per comparison)
############################################################
go_plot <- function(label, ylab_txt, text_size) {
  lst <- list()
  for (i in seq_along(comps)) {
    f <- sprintf("GO_DA_%s_%s.tsv", comps[i], label)
    if (!file.exists(f)) next
    df <- read.table(f, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    if (nrow(df) == 0) next
    df <- df[order(df$p_value), ][1:min(10, nrow(df)), ]
    df$comparison <- comp_labels[i]; lst[[comps[i]]] <- df
  }
  d <- bind_rows(lst)
  d$comparison <- factor(d$comparison, levels = comp_labels)
  d$term_name <- factor(d$term_name, levels = rev(unique(d$term_name)))
  ggplot(d, aes(x = term_name, y = fold_enrichment)) +
    geom_point(aes(size = intersection_size, colour = -log10(p_value))) +
    scale_color_gradient(low = OK, high = VM, limits = c(1, NA)) +
    coord_flip() + theme_bw() + facet_grid(~comparison, scales = "free_y", space = "free_y") +
    xlab(ylab_txt) + ylab("Fold enrichment") + labs(color = "-log10\nadjusted\np-value", size = "Number\nof genes") +
    theme(axis.ticks.length = unit(-0.1, "cm"), axis.text.x = element_text(size = 22), axis.text.y = element_text(size = text_size),
          axis.title = element_text(size = 30, face = "bold"), axis.text = element_text(color = "black"),
          panel.grid.minor = element_blank(), strip.text = element_text(face = "italic", size = 23),
          legend.title.align = 0.5, legend.text = element_text(size = 21), legend.title = element_text(size = 22))
}
ggsave(file.path(OUT, "Fig3c.GO_gain.png"), go_plot("gain", "GO Biological Process\n(genes gaining accessibility)", 19), units = "in", width = 15.5, height = 7.5, dpi = 300)
ggsave(file.path(OUT, "Fig3d.GO_loss.png"), go_plot("loss", "GO Biological Process\n(genes losing accessibility)", 21), units = "in", width = 13, height = 7, dpi = 300)

############################################################
# Multipanel: A on top; B (tall) left; C over D on the right
############################################################
lab <- function(img, l, w, strip = 220, size = 180) {
  img <- image_border(image_scale(img, as.character(w)), "white", paste0("0x", strip))
  image_annotate(img, l, gravity = "northwest", location = "+20+20", size = size, font = "Arial", weight = 700, color = "black")
}
W <- 4300; wB <- 2050; gap <- 60
imgA <- lab(image_read(file.path(OUT, "Fig3a.volcano.png")), "A", W)
imgB <- lab(image_border(image_trim(image_read(file.path(OUT, "Fig3b.heatmap.png"))), "white", "30x0"), "B", wB)
wR <- W - wB - gap
imgC <- lab(image_read(file.path(OUT, "Fig3c.GO_gain.png")), "C", wR)
imgD <- lab(image_read(file.path(OUT, "Fig3d.GO_loss.png")), "D", wR)
right <- image_append(c(imgC, imgD), stack = TRUE)
hmax <- max(image_info(imgB)$height, image_info(right)$height)
pad <- function(img, h) { i <- image_info(img); if (i$height < h) image_extent(img, paste0(i$width, "x", h), gravity = "north", color = "white") else img }
row2 <- image_append(c(pad(imgB, hmax), image_blank(gap, hmax, "white"), pad(right, hmax)), stack = FALSE)
image_write(image_append(c(imgA, row2), stack = TRUE), file.path(OUT, "Fig3_multipanel.png"), format = "png")

# Fig. 4: TF footprint-supported motif profiling of chromatin accessibility in A. calliptera embryos.
#   A  enrichment of motif-predicted binding sites in genomic features (GAT; 05_motif_enrichment_in_annotations_GAT.sh)
#   B  enrichment of candidate TF motifs in promoters of genes with accessible promoters (06_motif_enrichment_in_promoters.sh)
#   C  k-means clustering of mean motif bit-scores (row Z-scores) across the three stages
#   D  the ten most stage-invariant ('conserved') and ten most stage-variable ('diverged') TFs
# Also writes the elbow (within-cluster sum of squares) and average silhouette-width plots used to choose k.
#
# Usage: Rscript 07_plot_Fig4.R <indir> <outdir>
#   indir must contain: gatnormed_Collated_motif_AllAnnot_ChrBgd.tsv, Acal_MOTIFS_OUTPUT_details.simp.txt,
#   All_MOTIFS_OUTPUT_details.simp.cand.txt, Ac_{3dpf,7dpf,12dpf}_mpbs_TC.geneprom.genomeavgbitscore.bed

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(reshape2); library(ggh4x)
  library(ComplexHeatmap); library(circlize); library(cluster); library(factoextra); library(magick)
})

args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]; outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(indir)

K_SILHOUETTE <- 10     # clusters shown in the silhouette plot (supplementary figure)
K_HEATMAP    <- 7      # clusters shown in Fig. 4C (elbow plot and silhouette widths)

############################################################
# Fig. 4A: motif-site enrichment in genomic features
############################################################
gat <- read.table("gatnormed_Collated_motif_AllAnnot_ChrBgd.tsv", header = TRUE)
gat$annotation <- factor(gat$annotation,
  levels = c("aCNE", "hCNE", "intergenic", "5kb_gene_promoter", "5kb_gene_promoter_aCNE", "5kb_gene_promoter_hCNE", "five_prime_utr", "exon", "intron", "three_prime_utr"),
  labels = c("aCNE", "hCNE", "Intergenic", "5kb gene promoter", "aCNE in 5kb gene promoter", "hCNE in 5kb gene promoter", "5' UTR", "Exon", "Intron", "3' UTR"))
gat$tissue <- factor(gat$tissue, levels = c("3dpf", "7dpf", "12dpf"))
gat <- subset(gat, qvalue < 0.05 & annotation != "Exon")

p4a <- ggplot(gat, aes(x = annotation, y = fold)) +
  geom_point(aes(size = overlap_nsegments, colour = qvalue)) +
  scale_color_gradient(low = "#009E73", high = "#D55E00", limits = c(0, NA)) +
  coord_flip() + theme_bw() + facet_grid(~tissue, margins = FALSE) +
  xlab("Annotated region") + ylab("Fold enrichment") +
  labs(color = "q-value\n(FDR<0.05)", size = "Number\nof TF footprints\noverlapping\nannotation") +
  theme(axis.ticks.length = unit(-0.1, "cm"),
        axis.text.x = element_text(margin = margin(5, 5, 0, 5, "pt"), size = 16), axis.text.y = element_text(margin = margin(5, 5, 5, 5, "pt"), size = 14),
        axis.title = element_text(size = 30, face = "bold"), axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "italic", size = 16),
        legend.title.align = 0.5, legend.text = element_text(size = 17), legend.title = element_text(size = 18))
ggsave(file.path(outdir, "Fig4a.motif_annotation_enrichment.png"), p4a, units = "in", width = 10, height = 5, dpi = 300)

############################################################
# Fig. 4B: candidate TF motifs enriched in promoters
############################################################
cand <- read.table("All_MOTIFS_OUTPUT_details.simp.cand.txt", header = FALSE, stringsAsFactors = TRUE, sep = "\t")
# columns: V2 stage, V3 TF motif, V5 FDR, V9 target promoters with the motif, V10 fold enrichment
cand$V3 <- as.factor(cand$V3)
cand$V2 <- factor(cand$V2, levels = c("3dpf", "7dpf", "12dpf"))
cand$V3 <- factor(cand$V3, levels = rev(levels(cand$V3)))

p4b <- ggplot(cand, aes(x = V3, y = V10)) +
  geom_point(aes(size = V9, colour = -log10(V5))) +
  scale_color_gradient(low = "#009E73", high = "#D55E00", limits = c(1, NA)) +
  coord_flip() + theme_bw() + facet_nested_wrap(~V2) +
  xlab("Transcription Factor (TF)\nmotif") + ylab("Fold enrichment") +
  labs(color = "-log10\nadjusted\np-value\n(FDR<0.05)", size = "Number\nof genes") +
  theme(axis.ticks.length = unit(-0.1, "cm"),
        axis.text.x = element_text(margin = margin(5, 5, 0, 5, "pt"), angle = 90, vjust = 0.5, size = 14), axis.text.y = element_text(margin = margin(5, 5, 5, 5, "pt"), size = 12),
        axis.title = element_text(size = 25, face = "bold"), axis.text = element_text(color = "black"),
        panel.grid.minor = element_blank(), legend.title.align = 0.5,
        strip.text = element_text(face = "italic", size = 14.5), legend.text = element_text(size = 16), legend.title = element_text(size = 17))
ggsave(file.path(outdir, "Fig4b.TFenrichmentcands.png"), p4b, units = "in", width = 10, height = 5, dpi = 300)

############################################################
# Fig. 4C: k-means clustering of motif activity (mean bit-score across gene promoters)
############################################################
read_bits <- function(f, stage) { d <- read.table(f, header = FALSE, stringsAsFactors = FALSE, sep = "\t"); colnames(d) <- c("motifID", "value", "species"); d$stage <- stage; d }
all_motifbit <- bind_rows(read_bits("Ac_3dpf_mpbs_TC.geneprom.genomeavgbitscore.bed", "3dpf"),
                          read_bits("Ac_7dpf_mpbs_TC.geneprom.genomeavgbitscore.bed", "7dpf"),
                          read_bits("Ac_12dpf_mpbs_TC.geneprom.genomeavgbitscore.bed", "12dpf"))

# elbow plot: total within-cluster sum of squares against the number of clusters
wssplot <- function(data, nc = 20, seed = 1234) {
  data <- as.matrix(data); mode(data) <- "numeric"
  data <- data[complete.cases(data), , drop = FALSE]
  nc <- min(nc, nrow(data) - 1)
  wss <- numeric(nc); wss[1] <- (nrow(data) - 1) * sum(apply(data, 2, var))
  for (i in 2:nc) { set.seed(seed); wss[i] <- sum(kmeans(data, centers = i, nstart = 25)$withinss) }
  plot(1:nc, wss, type = "b", xlab = "Number of Clusters", ylab = "Within groups sum of squares")
}

make_stage_matrix <- function(df) {
  mat <- dcast(df, motifID ~ stage, value.var = "value", fun.aggregate = mean)
  rownames(mat) <- mat$motifID; mat$motifID <- NULL
  mat <- as.matrix(mat); mode(mat) <- "numeric"
  mat <- mat[, c("3dpf", "7dpf", "12dpf"), drop = FALSE]
  mat <- t(scale(t(mat)))                                            # row Z-score
  mat[complete.cases(mat), , drop = FALSE]
}
motif_mat <- make_stage_matrix(all_motifbit)

png(file.path(outdir, "Ac_wsskmeans.png"), width = 1200, height = 1000)
wssplot(motif_mat, nc = 15)
dev.off()

set.seed(1234)
km_res <- kmeans(motif_mat, centers = K_SILHOUETTE, iter.max = 1000, nstart = 20, algorithm = "MacQueen")
sil <- silhouette(km_res$cluster, dist(motif_mat))
png(file.path(outdir, "Ac_avgsilhouettewidth.png"), width = 1200, height = 1000)
print(fviz_silhouette(sil, label = FALSE, print.summary = TRUE))
dev.off()

ht <- Heatmap(motif_mat, name = "Z-score", km = K_HEATMAP,
              col = colorRamp2(c(-2, 0, 2), c("#009E73", "white", "#D55E00")),
              show_row_names = FALSE, show_column_names = TRUE, cluster_columns = FALSE,
              row_gap = unit(2.5, "mm"), column_order = c("3dpf", "7dpf", "12dpf"), column_labels = c("3 dpf", "7 dpf", "12 dpf"),
              row_dend_reorder = TRUE, cluster_row_slices = TRUE,
              row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 11), row_title_gp = gpar(fontsize = 11),
              width = unit(4.5, "cm"))
png(file.path(outdir, "Fig4c.Ac_TFheatmap.png"), width = 1000, height = 2100, res = 300)
draw(ht)
dev.off()

############################################################
# Fig. 4D: the ten most conserved and ten most diverged TFs across stages
############################################################
extract_tf_name <- function(x) ifelse(grepl("^[^_]+_(HUMAN|MOUSE)\\.", x), sub("_.*", "", x), sub(".*\\.", "", x))
all_motifbit$TF_name <- extract_tf_name(all_motifbit$motifID)
collapsed <- all_motifbit %>% group_by(TF_name, stage) %>% dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>% dplyr::select(TF_name, value, stage)

tf_mat <- make_stage_matrix(collapsed %>% dplyr::mutate(motifID = TF_name) %>% dplyr::select(motifID, value, stage))
tf_stats <- data.frame(TF_name = rownames(tf_mat), sd = apply(tf_mat, 1, sd), range = apply(tf_mat, 1, function(x) max(x) - min(x)), stringsAsFactors = FALSE)
top_conserved <- tf_stats %>% arrange(sd, range) %>% slice_head(n = 10) %>% pull(TF_name)
top_diverged  <- tf_stats %>% arrange(desc(range), desc(sd)) %>% filter(!TF_name %in% top_conserved) %>% slice_head(n = 10) %>% pull(TF_name)

selected_mat <- tf_mat[c(top_conserved, top_diverged), c("3dpf", "7dpf", "12dpf"), drop = FALSE]
row_group <- factor(c(rep("Conserved", length(top_conserved)), rep("Diverged", length(top_diverged))), levels = c("Conserved", "Diverged"))

png(file.path(outdir, "Fig4d.Ac_selectedTFs_heatmap.png"), width = 1800, height = 2800, res = 300)
ht_top <- Heatmap(selected_mat, name = "Z-score", col = colorRamp2(c(-2, 0, 2), c("#009E73", "white", "#D55E00")),
                  split = row_group, cluster_row_slices = FALSE, cluster_columns = FALSE, column_order = c("3dpf", "7dpf", "12dpf"),
                  show_row_names = TRUE, row_labels = rownames(selected_mat), show_column_names = TRUE, column_labels = c("3 dpf", "7 dpf", "12 dpf"),
                  row_names_gp = gpar(fontsize = 16), column_names_gp = gpar(fontsize = 19),
                  row_title_gp = gpar(fontsize = 19, fontface = "bold"), row_title_rot = 90,
                  width = unit(4.5, "cm"), row_gap = unit(2.5, "mm"),
                  heatmap_legend_param = list(title = "Z-score", title_gp = gpar(fontsize = 18, fontface = "bold"),
                                              labels_gp = gpar(fontsize = 16), legend_height = unit(4, "cm")))
draw(ht_top)
dev.off()
write.table(data.frame(TF_name = rownames(selected_mat), group = row_group), file.path(outdir, "Fig4d.Ac_selectedTFs_list.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

############################################################
# Multipanel
############################################################
add_label <- function(img, lab, width_px, top_strip = 170, label_size = 130) {
  img <- image_border(image_scale(img, as.character(width_px)), color = "white", geometry = paste0("0x", top_strip))
  image_annotate(img, text = lab, gravity = "northwest", location = "+20+20", size = label_size, font = "Arial", weight = 700, color = "black")
}
rd <- function(f) image_read(file.path(outdir, f))
bottom_row <- image_append(c(add_label(rd("Fig4c.Ac_TFheatmap.png"), "C", 1180), add_label(rd("Fig4d.Ac_selectedTFs_heatmap.png"), "D", 1180)), stack = FALSE)
fig4 <- image_append(c(add_label(rd("Fig4a.motif_annotation_enrichment.png"), "A", 2400),
                       add_label(rd("Fig4b.TFenrichmentcands.png"), "B", 2400), bottom_row), stack = TRUE)
image_write(fig4, file.path(outdir, "Fig4_multipanel.png"), format = "png")

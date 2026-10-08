# Fig. 5: candidate TF regulators and TF-target networks across embryonic stages.
#   A  hub TFs: the 15 motifs with the highest fold enrichment in gene promoters at each stage (06_motif_enrichment_in_promoters.sh)
#   B  stage-specific TF -> target-gene networks: footprints inside promoter aCNEs (04_promoter_footprints_and_edges.sh).
#      Only TFs linked to two or more targets are drawn (when >= 3 such TFs exist).
#
# Usage: Rscript 01_plot_Fig5.R <tf_dir> <promoter_bed> <outdir>
#   tf_dir        hubTFs_<stage>_full.tsv and Ac_<stage>_mpbs_TC.geneprom.aCNE.TF_TG_edges.txt
#   promoter_bed  10-column 5 kb promoter BED (gene_id in column 7, gene symbol in column 10), used to label targets

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(igraph); library(ggraph); library(ggrepel); library(patchwork); library(magick)
})
args <- commandArgs(trailingOnly = TRUE)
TFD <- args[1]; promoter_bed <- args[2]; OUT <- args[3]
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

OK <- "#009E73"; MID <- "#E69F00"; VM <- "#D55E00"
stages <- c("3dpf", "7dpf", "12dpf"); stage_labels <- c("3 dpf", "7 dpf", "12 dpf")

LBL   <- 4.3                # TF / gene label size (mm, ~12 pt)
TITLE <- 30                 # stage titles (pt)
KEY   <- LBL / 0.3528       # legend text, same size as the labels
HB    <- 6.4                # height of panel B (in)
AW    <- 2700               # width of panel A in the multipanel (px)

############################################################
# Panel A
############################################################
hub_df <- bind_rows(lapply(seq_along(stages), function(i) {
  df <- read.table(file.path(TFD, sprintf("hubTFs_%s_full.tsv", stages[i])), header = FALSE, sep = "\t",
                   col.names = c("TF", "target_count", "fold_enrichment"))
  df <- df[order(-df$fold_enrichment), ][1:15, ]
  df$stage <- stage_labels[i]; df
}))
hub_df$stage <- factor(hub_df$stage, levels = stage_labels)
tf_order <- hub_df %>% group_by(TF) %>% summarise(m = mean(fold_enrichment)) %>% arrange(m) %>% pull(TF)
hub_df$TF <- factor(hub_df$TF, levels = tf_order)

p_a <- ggplot(hub_df, aes(x = TF, y = fold_enrichment, fill = stage)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  coord_flip() +
  scale_fill_manual(values = c("3 dpf" = OK, "7 dpf" = MID, "12 dpf" = VM)) +
  theme_bw() + xlab("Transcription factor (TF) motif") + ylab("Fold enrichment (promoter target count)") + labs(fill = "Stage") +
  theme(axis.ticks.length = unit(-0.1, "cm"), axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 18, face = "bold"), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14), legend.title = element_text(size = 15))
ggsave(file.path(OUT, "Fig5a.hubTF_foldenrichment.png"), p_a, units = "in", width = 8, height = 7, dpi = 300)

############################################################
# Panel B
############################################################
# gene ID -> symbol from the promoter annotation; two genes without a symbol were named manually
prom <- read.table(promoter_bed, sep = "\t", quote = "", stringsAsFactors = FALSE,
                   col.names = c("chrom", "start0", "end", "name", "score", "strand", "geneID", "c8", "c9", "symbol"))
lookup <- setNames(prom$symbol, prom$geneID)
lookup["ENSACLG00000019367"] <- "Tcb1"
lookup["ENSACLG00000017959"] <- "IGKC-like"

# ring layout for the 7 dpf network: TFs with a single target sit on a ring around that target,
# TFs linking several targets sit at the centroid of their targets
radial_xy <- function(g) {
  el <- as_data_frame(g, what = "edges"); vn <- V(g)$name
  tg <- unique(el$to); tfs <- setdiff(vn, tg)
  single <- function(f) sum(el$from == f) == 1
  n_single <- sapply(tg, function(t) sum(el$to == t & sapply(el$from, single)))
  tg <- tg[order(-n_single)]
  Rg <- pmax(1.6, 0.115 * n_single[tg])
  cx <- numeric(length(tg)); cx[1] <- 0
  if (length(tg) > 1) for (i in 2:length(tg)) cx[i] <- cx[i - 1] + Rg[i - 1] + Rg[i] + 3.4
  xy <- matrix(0, nrow = length(vn), ncol = 2, dimnames = list(vn, c("x", "y")))
  xy[tg, ] <- cbind(cx, 0)
  for (i in seq_along(tg)) {
    mine <- el$from[el$to == tg[i] & sapply(el$from, single)]
    mine <- mine[order(tolower(mine))]
    ang <- pi / 2 - 2 * pi * (seq_along(mine) - 1) / max(1, length(mine))
    xy[mine, ] <- cbind(cx[i] + Rg[i] * cos(ang), Rg[i] * sin(ang))
  }
  for (f in tfs[sapply(tfs, function(f) sum(el$from == f) > 1)]) xy[f, ] <- colMeans(xy[el$to[el$from == f], , drop = FALSE])
  xy[vn, ]
}

make_network <- function(stage, label, radial = FALSE) {
  edges <- readLines(file.path(TFD, sprintf("Ac_%s_mpbs_TC.geneprom.aCNE.TF_TG_edges.txt", stage)))
  edges <- edges[edges != ""]
  tf     <- sub("-ENSACLG.*$", "", edges)
  rest   <- sub("^.*?-(ENSACLG.*)$", "\\1", edges)
  geneid <- sub(":.*", "", rest)
  symbol <- sub("^[^:]*:", "", rest)
  target <- ifelse(symbol %in% c("ensembl", ""), ifelse(!is.na(lookup[geneid]) & lookup[geneid] != "", lookup[geneid], geneid), symbol)
  df <- unique(data.frame(from = tf, to = target, stringsAsFactors = FALSE))
  g <- graph_from_data_frame(df, directed = TRUE)
  V(g)$type <- ifelse(V(g)$name %in% df$from, "TF", "Target gene")
  deg <- degree(g, mode = "out"); hub <- names(deg[deg >= 2])
  if (length(hub) >= 3) {
    keep <- df[df$from %in% hub, ]
    g <- graph_from_data_frame(keep, directed = TRUE)
    V(g)$type <- ifelse(V(g)$name %in% keep$from, "TF", "Target gene")
  }
  set.seed(42)
  lay <- if (radial) { xy <- radial_xy(g); ggraph(g, layout = "manual", x = xy[, 1], y = xy[, 2]) } else ggraph(g, layout = "fr")
  lay +
    geom_edge_link(alpha = if (radial) 0.3 else 0.4, colour = "grey50", arrow = arrow(length = unit(2.2, "mm")), end_cap = circle(2.5, "mm")) +
    geom_node_point(aes(colour = type, size = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = LBL, fontface = "bold", max.overlaps = Inf,
                   box.padding = 0.25, point.padding = 0.2, segment.size = 0.25, segment.color = "grey40", min.segment.length = 0) +
    scale_colour_manual(values = c("TF" = VM, "Target gene" = OK), name = "Node type") +
    scale_size_manual(values = if (radial) c("TF" = 4, "Target gene" = 6) else c("TF" = 7, "Target gene" = 5), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.15)) +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    labs(title = label) + theme_void() +
    theme(plot.title = element_text(face = "bold", size = TITLE, hjust = 0.5, margin = margin(b = 4)),
          legend.text = element_text(size = KEY), legend.title = element_text(size = KEY, face = "bold"))
}

pB <- (make_network("3dpf", "3 dpf") | make_network("7dpf", "7 dpf", radial = TRUE) | make_network("12dpf", "12 dpf")) +
  plot_layout(widths = c(1, 2.6, 0.9), guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.size = unit(7, "mm"))
ggsave(file.path(OUT, "Fig5b.TF_target_networks.png"), pB, width = 13, height = HB, units = "in", dpi = 300, bg = "white")

############################################################
# Multipanel
############################################################
W <- 3900; strip <- 190; lsz <- 150
lab <- function(img, l) {
  img <- image_border(img, "white", paste0("0x", strip))
  image_annotate(img, l, gravity = "northwest", location = "+20+10", size = lsz, font = "Arial", weight = 700, color = "black")
}
A0 <- image_scale(image_read(file.path(OUT, "Fig5a.hubTF_foldenrichment.png")), as.character(AW))
A0 <- image_extent(A0, paste0(W, "x", image_info(A0)$height), gravity = "center", color = "white")
fig <- image_append(c(lab(A0, "A"), lab(image_read(file.path(OUT, "Fig5b.TF_target_networks.png")), "B")), stack = TRUE)
image_write(fig, file.path(OUT, "Fig5_multipanel.png"), format = "png")

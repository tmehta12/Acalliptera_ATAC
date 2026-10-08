#!/usr/bin/env Rscript
# GO biological-process enrichment of genes gaining or losing promoter accessibility in each comparison (Fig. 3C, D).
# g:Profiler (gprofiler2::gost), organism 'acalliptera', GO:BP, FDR-corrected, default annotated background.
#
# Usage: Rscript 04_go_enrichment_gain_loss.R <outdir from 03_gene_level_direction.R>
# Output: GO_DA_<comp>_<gain|loss>.tsv (including the 'intersection' column of query genes in each term)

suppressPackageStartupMessages(library(gprofiler2))
setwd(commandArgs(trailingOnly = TRUE)[1])

for (comp in c("3v7dpf", "3v12dpf", "7v12dpf")) {
  for (label in c("gain", "loss")) {
    genes <- readLines(sprintf("DA_%s_%s_genes.txt", comp, label))
    if (length(genes) < 5) next
    res <- gost(query = genes, organism = "acalliptera", ordered_query = FALSE, significant = TRUE, user_threshold = 0.05,
                correction_method = "fdr", sources = "GO:BP", evcodes = TRUE)
    if (is.null(res)) next
    df <- res$result
    df$fold_enrichment <- (df$intersection_size / df$query_size) / (df$term_size / df$effective_domain_size)
    df <- df[order(df$p_value), c("term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size",
                                  "effective_domain_size", "fold_enrichment", "intersection")]
    write.table(df, sprintf("GO_DA_%s_%s.tsv", comp, label), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

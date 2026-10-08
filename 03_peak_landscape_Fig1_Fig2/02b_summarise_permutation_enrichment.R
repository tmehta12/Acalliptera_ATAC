# Summarise the permutation test (02_feature_enrichment_permutation.sh) into fold enrichment, p- and q-values.
#   fold enrichment = observed / mean(null)
#   empirical p     = (number of null counts at least as extreme as observed + 1) / (permutations + 1)
#                     ('at least as extreme' = >= observed for enrichment, <= observed for depletion)
#   q-value         = Benjamini-Hochberg across all stage x feature tests
# Usage: Rscript 02b_summarise_permutation_enrichment.R <feature_enrichment_dir>

args <- commandArgs(trailingOnly = TRUE)
dir  <- args[1]
obs  <- read.table(file.path(dir, "observed_summary.txt"), col.names = c("stage", "feature", "observed"))

res <- do.call(rbind, lapply(seq_len(nrow(obs)), function(i) {
  null_counts <- scan(file.path(dir, sprintf("null_%s_%s.txt", obs$stage[i], obs$feature[i])), quiet = TRUE)
  expected <- mean(null_counts)
  fold     <- obs$observed[i] / expected
  extreme  <- if (fold >= 1) sum(null_counts >= obs$observed[i]) else sum(null_counts <= obs$observed[i])
  data.frame(stage = obs$stage[i], feature = obs$feature[i], observed = obs$observed[i],
             expected = expected, fold_enrichment = fold, pvalue = (extreme + 1) / (length(null_counts) + 1))
}))
res$qvalue <- p.adjust(res$pvalue, method = "BH")
write.table(res, file.path(dir, "feature_enrichment_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

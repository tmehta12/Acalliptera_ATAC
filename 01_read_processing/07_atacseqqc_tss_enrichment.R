#!/usr/bin/env Rscript
# Step 7: ATAC-seq quality control with ATACseqQC (v1.18.0) and Tn5-shifted BAM files.
#   - shifts read start sites (+4 / -5 bp) -> <sample>.nochrM.nodup.filt.shifted.bam (input for TF footprinting)
#   - promoter/transcript-body (PT) score, nucleosome-free-region (NFR) score and TSS enrichment (TSSE) score
#   - splits reads into nucleosome-free / mono- / di- / tri-nucleosome fractions and plots the signal around TSSs
#
# Usage:
#   Rscript 07_atacseqqc_tss_enrichment.R <filtered.bam> <annotation.gtf> <outdir> <sample_id>
# e.g.
#   Rscript 07_atacseqqc_tss_enrichment.R 4.filtered/1aAc_3dpf_ATAC.nochrM.nodup.filt.bam \
#           Astatotilapia_calliptera.fAstCal1.2.101.gtf 7.qc/1aAc_3dpf_ATAC 1aAc_3dpf_ATAC
#
# The genome is supplied to ATACseqQC as a BSgenome package forged from the fAstCal1.2 FASTA
# (BSgenome.Acal.Ensembl.fAstCal1.2; BSgenome::forgeBSgenomeDataPkg with a seed file naming the primary-assembly sequences).
# Requires R >= 4.0 and: rtracklayer, GenomicFeatures, ensembldb, GenomicRanges, ATACseqQC, ChIPpeakAnno,
# GenomicAlignments, Rsamtools, BSgenome, Biostrings, BSgenome.Acal.Ensembl.fAstCal1.2

suppressPackageStartupMessages({
  library(rtracklayer); library(GenomicFeatures); library(GenomicRanges); library(ATACseqQC)
  library(ChIPpeakAnno); library(GenomicAlignments); library(Rsamtools); library(BSgenome)
  library(Biostrings); library(BSgenome.Acal.Ensembl.fAstCal1.2)
})

args       <- commandArgs(trailingOnly = TRUE)
bamfile    <- args[1]
gtf        <- args[2]
outPath    <- args[3]
sample_id  <- args[4]
genome     <- getBSgenome("BSgenome.Acal.Ensembl.fAstCal1.2")
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)

txs <- transcripts(makeTxDbFromGRanges(import(gtf)))

# --- shift reads (input BAM must be sorted, indexed and not already shifted) ---
possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]), paste0(possibleTag[2, ], possibleTag[1, ]))
bamTop100   <- scanBam(BamFile(bamfile, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]$tag
tags        <- names(bamTop100)[lengths(bamTop100) == 100]
gal         <- readBamFile(bamfile, tag = tags, asMates = TRUE, bigFile = TRUE)
shifted_bam <- file.path(outPath, sub("\\.bam$", ".shifted.bam", basename(bamfile)))
gal1        <- shiftGAlignmentsList(gal, outbam = shifted_bam)

# --- PT score and NFR score ---
pt <- PTscore(gal1, txs)
tiff(file.path(outPath, paste0(sample_id, "_1-PTscore.tiff")), units = "in", width = 5, height = 5, res = 100)
plot(pt$log2meanCoverage, pt$PT_score, xlab = "log2 mean coverage", ylab = "Promoter vs Transcript",
     main = "Promoter/Transcript body (PT) score")
dev.off()

nfr <- NFRscore(gal1, txs)
tiff(file.path(outPath, paste0(sample_id, "_2-NFRscore.tiff")), units = "in", width = 5, height = 5, res = 100)
plot(nfr$log2meanCoverage, nfr$NFR_score, xlab = "log2 mean coverage", ylab = "Nucleosome Free Regions score",
     main = "NFRscore for 200bp flanking TSSs", xlim = c(-10, 0), ylim = c(-5, 5))
dev.off()

# --- TSS enrichment: mean TSSE < 5 'concerning', 5-7 'acceptable', >= 7 'ideal' ---
tsse    <- TSSEscore(gal1, txs)
summ    <- data.frame(lapply(summary(tsse$TSSEscore), function(x) t(data.frame(x))))
rownames(summ) <- sample_id
summ$cutoff <- ifelse(summ$Mean < 5, "concerning", ifelse(summ$Mean < 7, "acceptable", "ideal"))
summ <- cbind(sample = rownames(summ), data.frame(summ, row.names = NULL))
write.table(summ, file.path(outPath, paste0(sample_id, "_3-TSSscore.txt")), quote = FALSE, sep = "\t", row.names = FALSE)

# --- nucleosome-free / mono- / di- / tri-nucleosome fractions and signal around TSSs ---
objs <- splitGAlignmentsByCut(gal1, txs = txs, genome = genome, outPath = outPath)
null <- writeListOfGAlignments(objs, outPath)
bamfiles <- file.path(outPath, c("NucleosomeFree.bam", "mononucleosome.bam", "dinucleosome.bam", "trinucleosome.bam"))

try({
  tiff(file.path(outPath, paste0(sample_id, "_4-cumulativepercscore.tiff")), units = "in", width = 5, height = 5, res = 100)
  cumulativePercentage(bamfiles[1:2], as(seqinfo(genome), "GRanges"))
  dev.off()
})

TSS         <- unique(promoters(txs, upstream = 0, downstream = 1))
librarySize <- estLibSize(bamfiles)
NTILE <- 101; ups <- dws <- 1010
sigs <- enrichedFragments(gal = objs[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")],
                          TSS = TSS, librarySize = librarySize, seqlev = seqlevels(TSS), TSS.filter = 0.5,
                          n.tile = NTILE, upstream = ups, downstream = dws)
sigs.log2 <- lapply(sigs, function(x) log2(x + 1))
tiff(file.path(outPath, paste0(sample_id, "_5-logtransformedTSSsignalheatmap.tiff")), units = "in", width = 5, height = 7, res = 100)
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width = ups + dws), zeroAt = .5, n.tile = NTILE)
dev.off()

out   <- featureAlignedDistribution(sigs, reCenterPeaks(TSS, width = ups + dws), zeroAt = .5, n.tile = NTILE, type = "l",
                                    ylab = "Averaged coverage")
range01 <- function(x) (x - min(x)) / (max(x) - min(x))
out   <- apply(out, 2, range01)
tiff(file.path(outPath, paste0(sample_id, "_6-rescaledTSSsignal.tiff")), units = "in", width = 8, height = 6, res = 100)
matplot(out, type = "l", xaxt = "n", xlab = "Position (bp)", ylab = "Fraction of signal")
axis(1, at = seq(0, 100, by = 10) + 1, labels = c("-1K", seq(-800, 800, by = 200), "1K"), las = 2)
abline(v = seq(0, 100, by = 10) + 1, lty = 2, col = "gray")
dev.off()

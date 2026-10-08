# Chromatin accessibility during early embryogenesis in *Astatotilapia calliptera*

Code accompanying **"Genome-wide chromatin accessibility during early embryogenesis in a generalist cichlid fish"**
(Midhukrishna, Smith, Joyce, Di-Palma, Haerty and Mehta).

ATAC-seq was carried out on *A. calliptera* ('Ruvuma') embryos at **3, 7 and 12 days post fertilisation (dpf)**, two
biological replicates per stage, each with a matched naked-DNA (gDNA) control, aligned to the fAstCal1.2 genome.
This repository contains the scripts used for each part of the Results, from raw reads to the figures.

* Steps that need a cluster are SLURM shell scripts (`.sh`); downstream statistics and plotting are R scripts (`.R`).
* The scripts are numbered in the order they are run; each one lists its inputs and outputs in its header.
* This is a record of the analysis for the paper rather than a packaged pipeline.

## Layout

| Folder | What it does | Paper |
|---|---|---|
| [`01_read_processing`](01_read_processing) | Trimming, alignment, mitochondrial-read removal, filtering, MACS2 peak calling, IDR / FRiP, TSS enrichment (ATACseqQC) and Tn5-shifted BAMs | Methods: ATAC-seq processing; Fig. 1 |
| [`02_genome_annotation`](02_genome_annotation) | Genomic feature annotation (exons, introns, UTRs, intergenic, 5 kb promoters, conserved non-coding elements) and assignment of peaks to features | Methods: peak annotation |
| [`03_peak_landscape_Fig1_Fig2`](03_peak_landscape_Fig1_Fig2) | Peak totals, distance to TSS, enrichment in genomic features, GO enrichment of promoter peaks | Fig. 1, Fig. 2 |
| [`04_differential_accessibility_Fig3`](04_differential_accessibility_Fig3) | csaw / edgeR differential accessibility between stages, region merging, gene-level direction, GO enrichment, robustness checks | Fig. 3 |
| [`05_TF_footprinting_Fig4`](05_TF_footprinting_Fig4) | HINT-ATAC footprinting, motif matching, promoter footprints, motif enrichment (GAT, enrichAnalyzer), clustering of TF activity | Fig. 4 |
| [`06_TF_networks_Fig5`](06_TF_networks_Fig5) | Hub TFs and stage-specific TF-target networks within promoter aCNEs | Fig. 5 |
| [`data`](data) | Small input tables supplied with the repository | |

Shared settings (paths, sample names, key parameters) are in [`config.sh`](config.sh); the shell scripts source it.
Software versions are listed in [`environment.md`](environment.md).

## Running the analysis

1. **Set up.** Edit `config.sh` (paths to the genome FASTA, Ensembl GTFs, the conserved-non-coding-element BED, and the
   working directory). The `#SBATCH` partition names in the shell scripts (`ei-medium`, `ei-largemem`) are those of the
   Earlham Institute cluster; change them for your system. Software is loaded with `module`/`source` commands as used on
   that cluster.
2. **Reads to peaks** - `01_read_processing/01` ... `07`, run in order (each is a SLURM array over the libraries, except
   `06_idr_and_frip.sh`). Step 3 also needs the mitochondrial genome FASTA (NCBI NC_018560.1). `07` produces the
   Tn5-shifted, sorted BAM files used by `04_*` and `05_*`.
3. **Annotation** - `02_genome_annotation/01_build_feature_annotation.sh`, then `02_annotate_peaks.sh`.
   The conserved non-coding elements are an input (see [`CNE_methods.md`](02_genome_annotation/CNE_methods.md)).
4. **Fig. 1 and 2** - `03_*`: `01` and `02` (cluster) then `02b` (R); `04_plot_Fig1.R`. For Fig. 2, `03_promoter_gene_lists_for_GO.sh`
   writes the gene lists that were submitted to the g:Profiler web server (g:GOst; settings in the script header);
   the GO:BP rows of its exports are in `data/GO_BP_promoter_peaks_<stage>.tsv` and are read by `05_plot_Fig2.R`.
   Panel A of Fig. 1 is a separately drawn schematic supplied as a PNG.
5. **Fig. 3** - `04_*`: `01` (csaw, a few hours on a workstation or cluster node) to `05`, then `06_plot_Fig3.R`.
   The promoter BED needed by `04_*` is the `5kb_gene_promoter` rows of the annotation:
   `awk '$4=="5kb_gene_promoter"' ${ANNOT_BED} > promoters.5kb.bed`.
6. **Fig. 4 and 5** - `05_*`: `01` to `06` (cluster), `07_plot_Fig4.R`; then `06_TF_networks_Fig5/01_plot_Fig5.R`.
   RGT configuration is described in [`rgt_data_config.md`](05_TF_footprinting_Fig4/rgt_data_config.md).

## Data

* Raw ATAC-seq and naked-DNA sequencing reads: *[accession to be added]*.
* Reference genome and annotation: Ensembl, *Astatotilapia_calliptera* fAstCal1.2 (releases 100 and 101).
* Processed peaks, differential-accessibility tables and supplementary tables: *[Zenodo DOI to be added]*.

## Main tools

Trim Galore, Bowtie2, SAMtools, Sambamba, bedtools, MACS2, IDR, Picard, ATACseqQC, UCSC tools, csaw, edgeR, g:Profiler
(`gprofiler2`), GAT, HINT-ATAC (Regulatory Genomics Toolbox), enrichAnalyzer, R (ggplot2, ComplexHeatmap, ggraph,
patchwork, magick). See [`environment.md`](environment.md).

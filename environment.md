# Software

## Cluster steps (`01_read_processing`, `02_genome_annotation`, `05_TF_footprinting_Fig4`)

Versions as loaded in the scripts (`module` / `source` lines).

| Tool | Version | Used in |
|---|---|---|
| Trim Galore! / FastQC | 0.5.0 / 0.11.9 | 01_read_processing/01 |
| Bowtie2 | 2.2.6 | 01_read_processing/02 |
| SAMtools | 1.9 (alignment, all later steps) | 01_read_processing/02-07, 05_* |
| BLAST | 2.3.0 | 01_read_processing/03 |
| Sambamba | 0.6.5 | 01_read_processing/04, 05_*/01 |
| Picard | 1.140 | 01_read_processing/04 |
| bedtools | 2.30.0 (peak calling, IDR/FRiP, annotation, footprints) | 01_read_processing/05-06, 02_*, 05_* |
| MACS2 | 2.1.1 | 01_read_processing/05 |
| UCSC tools (bedClip, bedToBigBed) | v333 | 01_read_processing/05 |
| IDR | 2.0.4 | 01_read_processing/06 |
| R / ATACseqQC / BSgenome | 4.5.2 / 1.18.0 | 01_read_processing/07 |
| Python | 3.5 (HPC), 3.x | 01_read_processing/03, 02_*/promoters_from_gene_bed.py |
| Regulatory Genomics Toolbox (HINT-ATAC, `rgt-motifanalysis`) | 0.13.0 | 05_*/02-03 |
| GAT | conda `gat` (python 3.10, numpy < 1.24) | 05_*/05 |
| enrichAnalyzer | https://github.com/Roy-lab/enrichAnalyzer_Nongraph | 05_*/06 |
| PHAST (phyloFit, phastCons, phyloP), Cactus | 1.5; 2.0.3 | CNE_methods.md |

## R analysis and plotting (`03` - `06`)

R 4.5 and the following package versions.

| Package | Version | Package | Version |
|---|---|---|---|
| csaw | 1.42.0 | ggplot2 | 4.0.3 |
| edgeR | 4.6.3 | dplyr / tidyr / reshape2 | 1.2.1 / 1.3.2 / 1.4.5 |
| rtracklayer | 1.68.0 | ComplexHeatmap | 2.24.1 |
| GenomicRanges | 1.60.0 | circlize | 0.4.18 |
| Rsamtools | 2.24.1 | ggraph / igraph / ggrepel | 2.2.2 / 2.3.2 / 0.9.8 |
| gprofiler2 | 0.2.4 | patchwork / cowplot / ggh4x | 1.3.2 / 1.2.0 / 0.3.1 |
| cluster / factoextra | 2.1.8.2 / 2.0.0 | magick / png | 2.9.1 / 0.1.9 |


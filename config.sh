#!/bin/bash
# Shared settings for the HPC (SLURM) scripts. Edit the paths below once, then every script sources this file.

# ---- Working directories (edit) --------------------------------------------------------------
export PROJECT=/path/to/AcalATAC                    # top-level working directory on the HPC
export RAW_READS=${PROJECT}/0.rawreads              # merged R1/R2 FASTQ files (see 01_read_processing/01_merge_and_trim.sh)
export RESULTS=${PROJECT}/results                   # per-step outputs are written below this directory

# ---- Reference files (edit) ------------------------------------------------------------------
export GENOME_ID=fAstCal1.2
export GENOME_FA=/path/to/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa
export GTF_R100=/path/to/Astatotilapia_calliptera.fAstCal1.2.100.gtf.gz   # Ensembl release 100: exon/intron/UTR/intergenic features
export GTF_R101=/path/to/Astatotilapia_calliptera.fAstCal1.2.101.gtf.gz   # Ensembl release 101: 5 kb promoters, ATACseqQC, footprint annotation
export CNE_BED=${PROJECT}/annotation/Astatotilapia_calliptera.fAstCal1.2.CNEsannot.bed  # hCNE/aCNE coordinates in fAstCal1.2 (see 02_genome_annotation/CNE_methods.md)
export CHROM_SIZES=${GENOME_FA}.chrom.sizes         # created in 02_genome_annotation/01_build_feature_annotation.sh
export MT_ACCESSION=NC_018560.1                     # A. calliptera mitochondrial genome (NCBI)
export ANNOT_BED=${PROJECT}/annotation/Astatotilapia_calliptera.fAstCal1.2.annot.bed   # built in 02_genome_annotation

# ---- Samples ---------------------------------------------------------------------------------
# Two biological replicates (a, b) per stage; each ATAC library has a matched naked-DNA (gDNA) control.
export STAGES="3dpf 7dpf 12dpf"
export ATAC_SAMPLES="1aAc_3dpf_ATAC 1bAc_3dpf_ATAC 2aAc_7dpf_ATAC 2bAc_7dpf_ATAC 3aAc_12dpf_ATAC 3bAc_12dpf_ATAC"
export GDNA_SAMPLES="1aAc_3dpf_gDNA 1bAc_3dpf_gDNA 2aAc_7dpf_gDNA 2bAc_7dpf_gDNA 3aAc_12dpf_gDNA 3bAc_12dpf_gDNA"

# ---- Parameters used in the paper ------------------------------------------------------------
export BOWTIE2_K=4                 # bowtie2 -k 4 -X2000 --mm (ENCODE ATAC-seq settings)
export THREADS=32
export MACS2_P=0.05                # MACS2 -p 0.05 --nomodel --shift -75 --extsize 150
export IDR_THRESH=0.1              # IDR threshold for replicate peak sets
export FPR=0.0001                  # HINT-ATAC / motif matching false-positive rate

# ---- Software paths (edit) -------------------------------------------------------------------
export ENRICH_ANALYZER=/path/to/enrichAnalyzer/src/enrichAnalyzer     # https://github.com/Roy-lab/enrichAnalyzer_Nongraph
export RGT_MOTIFS=${HOME}/rgtdata/motifs                               # *.mtf motif repositories used by RGT (see 05_TF_footprinting_Fig4/rgt_data_config.md)

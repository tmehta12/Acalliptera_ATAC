# RGT data configuration for A. calliptera

HINT-ATAC and `rgt-motifanalysis` read genome and motif settings from `~/rgtdata/data.config.user`. Entries used for
fAstCal1.2 (replace paths with your own; the gene-alias and gene-list files are the Ensembl-release-101 annotation in
RGT format):

```ini
[fAstCal1.2]
genome: /path/to/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa
chromosome_sizes: /path/to/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes
genes_Gencode: /path/to/fAstCal1.2.genes_Gencode.bed      # Ensembl gene IDs
genes_RefSeq: /path/to/fAstCal1.2.genes_RefSeq.bed        # gene symbols
annotation: /path/to/Astatotilapia_calliptera.fAstCal1.2.101.gtf
gene_alias: /path/to/fAstCal1.2.alias.txt

[MotifData]
pwm_dataset: motifs
logo_dataset: logos
repositories: cichlidacCSsp, cichlidCW, cichlidJASPAR, jaspar_vertebrates, hocomoco, jaspar_plants, uniprobe_primary, uniprobe_secondary
```

* `cichlidacCSsp` - *A. calliptera*-specific PWMs and `cichlidCW` - cichlid-wide PWMs, both from the earlier multi-tissue
  cichlid ATAC-seq study; `cichlidJASPAR` - JASPAR vertebrate PWMs re-scored for cichlids; `jaspar_vertebrates`,
  `hocomoco` - vertebrate PWM repositories (the Methods additionally list GTRD and UniPROBE). Each repository is a
  folder of `.pwm` files plus a `.mtf` file in `~/rgtdata/motifs/`.
* In the custom `cichlid*.mtf` files, column 10 holds the score thresholds for false-positive rates
  0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001 as `4,8,12,16,20,24`; the analysis uses FPR 0.0001
  (`rgt-motifanalysis matching` default). Pseudocounts: 1.0.
* 2,042 motifs were loaded in the footprinting runs; 4,408 motif IDs map to 1,442 unique TF names
  (`allTFmotifs_motifname.txt`, built with `cat *.mtf | awk '{print $2,toupper($4)}' OFS='\t' | sort -u -k1,1`).

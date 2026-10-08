# Input tables supplied with the repository

| File | Description | Used by |
|---|---|---|
| `allTFmotifs_motifname.txt` | motif ID -> TF name (upper case) for the 4,408 motifs loaded in the footprinting runs (1,442 unique TF names) | `05_TF_footprinting_Fig4/06_motif_enrichment_in_promoters.sh` (regenerated there from the RGT `.mtf` files) |

Larger inputs are not stored here: raw reads, the genome FASTA and Ensembl annotation (see the main README), the
conserved-non-coding-element BED (see `02_genome_annotation/CNE_methods.md`) and the processed peak and differential
accessibility tables (Zenodo).

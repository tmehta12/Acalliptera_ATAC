# Conserved non-coding elements (CNEs) used in the annotation

The hCNE / aCNE coordinates used in `01_build_feature_annotation.sh` (`${CNE_BED}`) were generated in a companion
cichlid CNE analysis and are supplied to this repository as a finished BED file `Astatotilapia_calliptera.fAstCal1.2.CNEsannot.bed`; the full annotation file is also provided `Astatotilapia_calliptera.fAstCal1.2.annot.bed`; the full CNE-calling scripts are not duplicated here. The procedure, as described in the Methods, was:

1. **Multiple genome alignment** of six cichlid genomes - *A. calliptera* (fAstCal1.2), *M. zebra* (UMD2a),
   *P. nyererei* (v1), *A. burtoni* (v1), *N. brichardi* (v1) and *O. niloticus* (UMD_NMBU) - with Cactus (v2.0.3),
   written as MAF with *O. niloticus* as reference. Tree:
   `(((((Astatotilapia_calliptera,Metriaclima_zebra),Pundamilia_nyererei),Astatotilapia_burtoni),Neolamprologus_brichardi),Oreochromis_niloticus)`

2. **Neutral model** (PHAST v1.5):
   ```bash
   phyloFit --tree "(((((Metriaclima_zebra,Astatotilapia_calliptera),Pundamilia_nyererei),Astatotilapia_burtoni),Neolamprologus_brichardi),Oreochromis_niloticus)" \
            --subst-mod REV --out-root neutral alignment.maf
   ```

3. **Split the MAF** by *O. niloticus* chromosome/scaffold (UCSC `mafSplit`), then predict CNEs per chromosome:
   ```bash
   phastCons --target-coverage 0.3 --expected-length 30 --most-conserved cnes.bed --estimate-trees trees \
             --msa-format MAF chromosome.maf neutral.mod
   ```
   Elements with sequence identity >= 90% over >= 30 bp are **hCNEs**; those with identity < 90% over >= 30 bp are
   **pseudo-aCNEs**.

4. **Acceleration test** of the pseudo-aCNEs with phyloP (likelihood-ratio test, CONACC mode):
   ```bash
   phyloP --method LRT --mode CONACC --features pseudo_aCNEs.bed neutral.mod chromosome.maf
   ```
   Pseudo-aCNEs with a negative CONACC score, likelihood ratio < 0.05 and `altsubscale` > 1 are **aCNEs**.

5. **Coordinates in *A. calliptera*.** CNEs are defined on *O. niloticus*; their *A. calliptera* coordinates were
   obtained from the pairwise *O. niloticus*-*A. calliptera* alignment blocks with `mafsInRegion` (UCSC tools) and
   bedtools (merging blocks closer than 5 bp, discarding mapped elements < 30 bp). pseudo-aCNEs were then dropped
   and each CNE was assigned its closest gene (`bedtools closest`).

The resulting file has ten columns: chrom, start, end, CNE type (`hCNE` or `aCNE`), score, strand, closest gene ID,
transcript (`.`), CNE details (total conservation score; length; average conservation; distance to gene) and
closest gene symbol.

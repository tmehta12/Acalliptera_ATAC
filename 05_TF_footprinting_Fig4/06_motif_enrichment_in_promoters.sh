#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 48000
#SBATCH -t 0-05:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Fig. 4B and Fig. 5A: enrichment of TF motifs (footprint-supported) in the promoters of genes with an accessible
# promoter footprint at each stage, using enrichAnalyzer (hypergeometric test, Benjamini-Hochberg FDR).
#
#   1. motif ID -> TF name table (allTFmotifs_motifname.txt) from the RGT motif repositories
#   2. footprint-supported motif matches in gene promoters, keeping gene, stage and motif
#   3. enrichAnalyzer input files: gene sets per stage (clusters) and gene -> TF motif annotations
#   4. enrichAnalyzer
#   5. tidy output: stage names, species column -> Acal_MOTIFS_OUTPUT_details.simp.txt
#
# Input : ${RESULTS}/10.footprinting/ac_fp/Ac_<stage>/Ac_<stage>_mpbs.bed, ${ANNOT_BED}
# Output: ${RESULTS}/10.footprinting/motif_enrichment/enrichment/output/Acal_MOTIFS_OUTPUT_details.simp.txt
#   columns: species, stage, TF motif, p-value, adjusted p-value (FDR), n genes with any term in the stage set,
#            n genes with the term, n target promoters, n target promoters with the term, fold enrichment
#            (the original 10th column listing the target gene IDs is dropped)

source "$(dirname "$0")/../config.sh"
ml gcc
ml zlib
source bedtools-2.30.0

FP=${RESULTS}/10.footprinting
WD=${FP}/motif_enrichment
mkdir -p "${WD}/enrichment"/{input1,input2,output} && cd "${WD}"

# 1. motif ID -> TF name (upper case); 4,408 motifs / 1,442 unique TF names
cat ${RGT_MOTIFS}/*.mtf | awk '{print $2,toupper($4)}' OFS='\t' | sort -u -k1,1 > allTFmotifs_motifname.txt

# 2. promoter footprints: chrom, start, end, motifID, bit-score, strand, gene_id, gene_name, feature, strand
for stage in ${STAGES}; do
  f=${FP}/ac_fp/Ac_${stage}/Ac_${stage}_mpbs.bed
  sed 's/\t$//' ${f} > Ac_${stage}_mpbs.notab.bed
  bedtools intersect -a "${ANNOT_BED}" -b Ac_${stage}_mpbs.notab.bed -wb | awk '$4=="5kb_gene_promoter"' \
    | awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > Ac_${stage}_mpbs.geneprom.bed
  rm Ac_${stage}_mpbs.notab.bed
done

# 3a. gene <-> stage (one line per promoter footprint)
: > Ac_speciesspecnames_tissueassign.txt
for stage in ${STAGES}; do awk -v s=${stage} '{print $7,s}' OFS='\t' Ac_${stage}_mpbs.geneprom.bed >> Ac_speciesspecnames_tissueassign.txt; done
# 3b. gene <-> TF name (via motif ID); motifs without a TF name are dropped
: > all_motifnames_regnet.tmp
for stage in ${STAGES}; do awk '{print $7,$4,"Ac"}' OFS='\t' Ac_${stage}_mpbs.geneprom.bed >> all_motifnames_regnet.tmp; done
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$2]){print $0,a[$2]}else{print $0,"NULL"}}' allTFmotifs_motifname.txt all_motifnames_regnet.tmp \
  | awk '{print $1,$4,$2,$3}' OFS='\t' > all_motifnames_regnet.txt
cut -f1,2 all_motifnames_regnet.txt > Acal_motifnames_regnet.txt
rm all_motifnames_regnet.tmp
# 3c. enrichAnalyzer annotation files: gene -> TF ('gotermap') and number of genes per TF ('genecnt')
grep -v NULL Acal_motifnames_regnet.txt | awk 'NF==2{printf("%s\t%s\t1\n",$1,$2)}' > enrichment/input1/Acal_gotermap.txt
grep -v NULL Acal_motifnames_regnet.txt | cut -f2 | sort | uniq -c | awk '{printf("%s\t%s\t1\n",$2,$1)}' > enrichment/input1/Acal_genecnt.txt
# 3d. gene sets per stage ('Cluster<stage>' <TAB> gene#gene#...; one entry per footprint)
cp Ac_speciesspecnames_tissueassign.txt enrichment/input2/Acal_speciesspecnames_tissueassign.txt
awk -F'\t' '{a[$2]=(a[$2]==""?$1:a[$2]"#"$1)} END{for(k in a) print "Cluster"k"\t"a[k]}' Ac_speciesspecnames_tissueassign.txt | sort -k1,1 \
  > enrichment/output/Acal_MOTIFS_INPUT.txt

# 4. enrichAnalyzer: <gene sets> <gene-set names> <annotation file prefix> <min overlap> <output prefix> persg
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${ENRICH_ANALYZER})/../x86_64/lib:$(dirname ${ENRICH_ANALYZER})/../x86_64/lib2
${ENRICH_ANALYZER} enrichment/output/Acal_MOTIFS_INPUT.txt enrichment/input2/Acal_speciesspecnames_tissueassign.txt \
  enrichment/input1/Acal_ 1 enrichment/output/Acal_MOTIFS_OUTPUT persg

# 5. tidy: replace 'Cluster<stage>' by the stage and add a species column; drop the target-gene list (column 10)
awk -v sp="A. calliptera" '{print sp,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' enrichment/output/Acal_MOTIFS_OUTPUT_details.txt \
  | sed 's|Cluster12dpf|12dpf|g' | sed 's|Cluster7dpf|7dpf|g' | sed 's|Cluster3dpf|3dpf|g' \
  > enrichment/output/Acal_MOTIFS_OUTPUT_details.simp.txt

# number of TF motifs (of the 1,442 input) enriched at any stage
cut -f3 enrichment/output/Acal_MOTIFS_OUTPUT_details.simp.txt | awk '!v[$0]++' | wc -l

# candidate TFs highlighted in Fig. 4B: factors with established stage-specific developmental roles
# (EGR1 neural activity; SOX10 neural crest; CRX and RXRA retinal / retinoic-acid development; HNF4A, FOXA1, FOXA2 tissue differentiation)
for tf in EGR1 CRX SOX10 RXRA HNF4A FOXA1 FOXA2 hoxc6 elf4; do
  grep -wiF ${tf} enrichment/output/Acal_MOTIFS_OUTPUT_details.simp.txt | grep -v '::' >> All_MOTIFS_OUTPUT_details.simp.cand.txt || true
done

# hub TFs for Fig. 5A: per stage, TFs ranked by fold enrichment among significant (FDR < 0.05) motifs
for stage in ${STAGES}; do
  awk -F'\t' -v s="${stage}" '$2==s && $5<0.05 {print $3"\t"$9"\t"$10}' enrichment/output/Acal_MOTIFS_OUTPUT_details.simp.txt \
    | sort -t$'\t' -k2,2nr > hubTFs_${stage}_full.tsv
done

#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 12000
#SBATCH -t 0-03:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Build the genomic feature annotation used to annotate ATAC-seq peaks and TF footprints in A. calliptera (fAstCal1.2).
#
# Features: exon, intron, intergenic, five_prime_utr, three_prime_utr, 5kb_gene_promoter, hCNE, aCNE
# Output  : ${ANNOT_BED}   10 columns: chrom, start, end, feature, score, strand, gene_id, transcript_id, exon_number, gene_name
#                          (0-based BED; a 1-based GFF-style copy is written alongside)
#
# Inputs  : ${GTF_R100}  Ensembl release 100 gene models -> exon, intron, intergenic, UTRs and 5 kb promoters
#           ${CNE_BED}   hCNE / aCNE coordinates mapped to fAstCal1.2 (see CNE_methods.md)

source "$(dirname "$0")/../config.sh"
source bedtools-2.30.0

SCRIPTS=$(cd "$(dirname "$0")" && pwd)
OUT=$(dirname "${ANNOT_BED}")
mkdir -p "${OUT}" && cd "${OUT}"
PFX=Astatotilapia_calliptera.fAstCal1.2

# 0. chromosome sizes (name, length)
samtools faidx "${GENOME_FA}"
awk '{print $1"\t"$2}' "${GENOME_FA}.fai" | sort -V -k1,1 > "${CHROM_SIZES}"

# 1. gene BED (chrom, start, end, feature, ., strand, gene_id, ., ., gene_name) from release 100
gzip -dc "${GTF_R100}" \
  | awk 'BEGIN{OFS="\t"} $3=="gene" {print $1,$4-1,$5,$3,".",$7,$10,".",".",$14}' \
  | sed 's|"||g' | sed 's|;||g' | sort -V -k1,1 -k2,2n > ${PFX}.100_gene.sorted.bed
# helper with a chr_start / chr_end key in front of each gene, used to name intergenic regions by their flanking genes
awk '{print $1"_"$2,$1"_"$3,$0}' OFS='\t' ${PFX}.100_gene.sorted.bed > ${PFX}.100_gene.bed

# 2. exons: merge overlapping exons, then map gene / transcript / exon / name back using the chr_start key
gzip -dc "${GTF_R100}" \
  | awk 'BEGIN{OFS="\t"} $3=="exon" {print $1,$4-1,$5,$3,".",$7,$10,$14,$18,$20}' \
  | sed 's|"||g' | sed 's|;||g' | sort -k1,1 -k2,2n | awk '{print $1"_"$2,$0}' OFS='\t' > ${PFX}.100_exon_pre-merged.bed
gzip -dc "${GTF_R100}" \
  | awk 'BEGIN{OFS="\t"} $3=="exon" {print $1,$4-1,$5,$3,".",$7,$10,$14,$18,$20}' \
  | sed 's|"||g' | sed 's|;||g' | sort -k1,1 -k2,2n | bedtools merge -i - | awk '{print $1"_"$2,$0}' OFS='\t' > ${PFX}.100_exon_merged.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$1]){print $0,a[$1]}else{print $0,"NULL"}}' \
  ${PFX}.100_exon_pre-merged.bed ${PFX}.100_exon_merged.bed | cut -f6-15 > ${PFX}.100_exon_merged.map.bed

# 3. introns = gene bodies minus exons
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,".",".",$10}' ${PFX}.100_gene.sorted.bed \
  | bedtools subtract -a stdin -b ${PFX}.100_exon_merged.map.bed | sed 's|gene|intron|g' > ${PFX}.100_intron.map.bed

# 4. 5' and 3' UTRs
for utr in five_prime_utr:5UTR three_prime_utr:3UTR; do
  gzip -dc "${GTF_R100}" \
    | awk -v f="${utr%%:*}" 'BEGIN{OFS="\t"} $3==f {print $1,$4-1,$5,$3,".",$7,$10,$14,".",$18}' \
    | sed 's|"||g' | sed 's|;||g' | sort -V -k1,1 -k2,2n > ${PFX}.100_${utr##*:}.map.bed
done

# 5. intergenic = complement of genes; each region is named after the genes on either side
gzip -dc "${GTF_R100}" \
  | awk 'BEGIN{OFS="\t"} $3=="gene" {print $1,$4-1,$5,$3,".",$7,$10,".",".",$14}' \
  | sed 's|"||g' | sed 's|;||g' | sort -V -k1,1 -k2,2n \
  | bedtools complement -i stdin -g "${CHROM_SIZES}" | awk '{print $1"_"$2,$1"_"$3,$0}' OFS='\t' > ${PFX}.100_intergenic.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$2]=$0;next}{if(a[$1]){print $0,a[$1]}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"}}' \
  ${PFX}.100_gene.bed ${PFX}.100_intergenic.bed \
  | awk '{print $1,$2,$3,$4,$5,"intergenic",".",$13,$14"_intergenic_",$15,$16,$17"_intergenic_"}' OFS='\t' > ${PFX}.100_intergenic.firstmap.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$0;next}{if(a[$2]){print $0,a[$2]}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"}}' \
  ${PFX}.100_gene.bed ${PFX}.100_intergenic.firstmap.bed \
  | awk '{print $3,$4,$5,$6,$7,$8,$9$21,$10,$11,$12$24}' OFS='\t' > ${PFX}.100_intergenic.map.bed
rm -f ${PFX}.100_intergenic.firstmap.bed

# 6. 5 kb promoters (strand-aware; shortened at neighbouring genes, see promoters_from_gene_bed.py); gene names are mapped back from the gene BED
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$7,$6}' ${PFX}.100_gene.sorted.bed > ${PFX}.100_genes_for_promoters.bed   # chrom, start, end, gene_id, strand
python3 "${SCRIPTS}/promoters_from_gene_bed.py" ${PFX}.100_genes_for_promoters.bed
# script output columns: chrom, start, end, strand, gene_id, gene_id -> chrom, start, end, gene_id, gene_id, strand (start < end)
awk 'BEGIN{OFS="\t"}{s=($2<$3)?$2:$3; e=($2<$3)?$3:$2; print $1,s,e,$5,$6,$4}' \
  ${PFX}.100_genes_for_promoters5kb_promoters.stranded.GENEBED.bed > ${PFX}.5kb_promoters.bed
awk 'BEGIN{OFS="\t"}NR==FNR{a[$9]=$0;next}{if(a[$4]){print $0,a[$4]}else{print $0,"NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL"}}' \
  ${PFX}.100_gene.bed ${PFX}.5kb_promoters.bed \
  | awk '{print $1,$2,$3,"5kb_gene_promoter",".",$6,$4,".",".",$18}' OFS='\t' > ${PFX}.5kb_promoters.map.bed

# 7. promoters lie within intergenic space, so remove them from the intergenic regions
bedtools subtract -a ${PFX}.100_intergenic.map.bed -b ${PFX}.5kb_promoters.map.bed > ${PFX}.100_intergenic.NEWminus5kbProm.map.bed

# 8. combine all features (CNEs overlap the other features by design) and sort
cat ${PFX}.100_intergenic.NEWminus5kbProm.map.bed "${CNE_BED}" ${PFX}.5kb_promoters.map.bed \
    ${PFX}.100_5UTR.map.bed ${PFX}.100_exon_merged.map.bed ${PFX}.100_intron.map.bed ${PFX}.100_3UTR.map.bed \
  | sort -V -k1,1 -k2,2n > "${ANNOT_BED}"

# 9. 1-based GFF-style copy: chrom, gene_name, feature, start+1, end, score, strand, exon_number, gene_id;transcript_id
awk '{print $1,$10,$4,$2+1,$3,$5,$6,$9,$7";"$8}' OFS='\t' "${ANNOT_BED}" > "${ANNOT_BED%.bed}.gff"

#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 18000
#SBATCH -t 0-02:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Footprints and predicted TF binding sites (TFBS) in gene promoters, and TF -> target-gene edges (Figs 4 and 5).
#   a. attach the footprint tag count to every motif-predicted binding site (bedtools intersect with the footprint BED)
#   b. keep sites in 5 kb gene promoters (annotation from 02_genome_annotation) -> Ac_<stage>_mpbs_TC.geneprom.bed
#   c. summary statistics per stage (footprints, TFBS, bit-score, tag count, TF-target edges)
#   d. candidate TF-target edges whose footprint lies in an aCNE inside the promoter -> *.aCNE.TF_TG_edges.txt
#      (one 'TF-geneID:symbol' per line), and pairwise shared / stage-unique edges between stages
#   e. mean motif bit-score per motif across all gene promoters, per stage (input to the clustering in 07_plot_Fig4.R)
#
# Input : ${RESULTS}/10.footprinting/ac_fp/Ac_<stage>/Ac_<stage>.bed and Ac_<stage>_mpbs.bed, ${ANNOT_BED}
# Output: ${RESULTS}/10.footprinting/promoter/

source "$(dirname "$0")/../config.sh"
source bedtools-2.30.0

FP=${RESULTS}/10.footprinting
OUT=${FP}/promoter
mkdir -p "${OUT}" && cd "${OUT}"

# a. + b. tag counts and promoter subset
for stage in ${STAGES}; do
  id=Ac_${stage}
  sed -i 's/[ \t]\+$//' ${FP}/ac_fp/${id}/${id}.bed ${FP}/ac_fp/${id}/${id}_mpbs.bed       # remove trailing tabs
  # column 5 = motif bit-score, column 11 = footprint tag count (number of reads)
  bedtools intersect -a ${FP}/ac_fp/${id}/${id}_mpbs.bed -b ${FP}/ac_fp/${id}/${id}.bed -wao > ${id}_mpbs_TC.bed
  bedtools intersect -a ${id}_mpbs_TC.bed -b "${ANNOT_BED}" -wb | awk '$17=="5kb_gene_promoter"' > ${id}_mpbs_TC.geneprom.bed
done

# c. summaries
: > Ac_total_footprints.tsv; : > Ac_promoter_footprints.tsv; : > Ac_promoter_tfbs_unique.tsv
: > Ac_promoter_mean_bitscore.tsv; : > Ac_promoter_mean_tagcount.tsv; : > Ac_promoter_TF_TG_edges.tsv
for stage in ${STAGES}; do
  id=Ac_${stage}
  echo -e "Ac\t${stage}\t$(cut -f10 ${id}_mpbs_TC.bed | awk '!v[$0]++' | wc -l)" >> Ac_total_footprints.tsv            # footprints in the genome
  echo -e "Ac\t${stage}\t$(cut -f10 ${id}_mpbs_TC.geneprom.bed | awk '!v[$0]++' | wc -l)" >> Ac_promoter_footprints.tsv   # footprints in promoters
  echo -e "Ac\t${stage}\t$(awk '{print $1"_"$2"_"$3}' ${id}_mpbs_TC.geneprom.bed | awk '!v[$0]++' | wc -l)" >> Ac_promoter_tfbs_unique.tsv   # unique TFBS positions
  echo -e "Ac\t${stage}\t$(cut -f5,10 ${id}_mpbs_TC.geneprom.bed | awk '!v[$2]++' | awk '{x+=$1} END{print x/NR}')" >> Ac_promoter_mean_bitscore.tsv
  echo -e "Ac\t${stage}\t$(awk '{print $1"_"$2"_"$3,$11}' OFS='\t' ${id}_mpbs_TC.geneprom.bed | awk '!v[$1]++' | awk '{x+=$2} END{print x/NR}')" >> Ac_promoter_mean_tagcount.tsv
  echo -e "Ac\t${stage}\t$(awk '{print $4"-"$20}' ${id}_mpbs_TC.geneprom.bed | awk '!v[$0]++' | wc -l)" >> Ac_promoter_TF_TG_edges.tsv    # TF-target gene relationships
done

# d. candidate TF-target edges with a footprint in a promoter aCNE
for stage in ${STAGES}; do
  id=Ac_${stage}_mpbs_TC.geneprom
  bedtools intersect -a ${id}.bed -b "${ANNOT_BED}" -wb | awk '$27=="aCNE"' > ${id}.aCNE.bed
  awk '{print $4"-"$20":"$23}' ${id}.aCNE.bed \
    | sed -E 's/^([A-Z0-9]+)_([A-Z]+)\..*-/\1-/; s/^MA[0-9]+\.[0-9]+\.([^-]+)-/\1-/' | sort -u > ${id}.aCNE.TF_TG_edges.txt
done

pairs=("3dpf 7dpf" "3dpf 12dpf" "7dpf 12dpf")
for p in "${pairs[@]}"; do
  set -- ${p}
  a=Ac_${1}_mpbs_TC.geneprom.aCNE.TF_TG_edges; b=Ac_${2}_mpbs_TC.geneprom.aCNE.TF_TG_edges
  comm -12 <(sort ${a}.txt) <(sort ${b}.txt) > ${a}_vs_${b}.shared.txt
  comm -23 <(sort ${a}.txt) <(sort ${b}.txt) > ${a}_unique_vs_${b}.txt
  comm -13 <(sort ${a}.txt) <(sort ${b}.txt) > ${b}_unique_vs_${a}.txt
done

# e. mean bit-score of every motif across all gene promoters, per stage
for stage in ${STAGES}; do
  id=Ac_${stage}_mpbs_TC.geneprom
  awk '{print $4,$5}' OFS='\t' ${id}.bed > ${id}.tmp
  awk '{a[$1]+=$2} END{for (i in a) print i,a[i]}' OFS='\t' ${id}.tmp > ${id}.tmp1          # summed bit-score per motif
  awk '{c[$1]++} END{for (k in c) print k,c[k]}' OFS='\t' ${id}.tmp > ${id}.tmp2             # number of sites per motif
  awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$1]){print $0,a[$1]}else{print $0,"NULL"}}' ${id}.tmp2 ${id}.tmp1 \
    | awk -v s="A. calliptera" '{print $1,$2/$3,s}' OFS='\t' > Ac_${stage}_mpbs_TC.geneprom.genomeavgbitscore.bed
  rm -f ${id}.tmp ${id}.tmp1 ${id}.tmp2
done

#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 48000
#SBATCH -t 0-23:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Fig. 4A: enrichment of TF footprint-predicted binding sites (motif matches within footprints) in genomic features,
# with the Genomic Association Tester (GAT, https://gat.readthedocs.io), BH-corrected q-values, genome (chromosome
# lengths) as workspace.
#
# GAT installation (conda): conda create -n gatenv python=3.10 gat -c conda-forge -c bioconda; conda install "numpy<1.24"
#
# Input : ${RESULTS}/10.footprinting/ac_fp/Ac_<stage>/Ac_<stage>_mpbs.bed, ${ANNOT_BED}, ${CHROM_SIZES}
# Output: ${RESULTS}/10.footprinting/GAT/gatnormed_Collated_motif_AllAnnot_ChrBgd.tsv

source "$(dirname "$0")/../config.sh"
eval "$(conda shell.bash hook)"; conda activate gatenv

FP=${RESULTS}/10.footprinting
OUT=${FP}/GAT
mkdir -p "${OUT}" && cd "${OUT}"

# segments: motif-predicted binding sites; annotations: genomic features (first six BED columns); workspace: chromosome lengths
cut -f1-6 "${ANNOT_BED}" > ac.annot.basic.bed
awk '{print $1,"0",$2}' OFS='\t' "${CHROM_SIZES}" > workspace.bed

for stage in ${STAGES}; do
  sort -k1,1 -k2,2n ${FP}/ac_fp/Ac_${stage}/Ac_${stage}_mpbs.bed > ac.${stage}_collated.ATAC_mpbs.bed
  gat-run.py --verbose=5 --log=gatnormed_ac_${stage}_motif_allAnnot_ChrBgd.tsv.log \
    --segments=ac.${stage}_collated.ATAC_mpbs.bed --annotations=ac.annot.basic.bed --workspace=workspace.bed \
    --counter=segment-overlap --ignore-segment-tracks --qvalue-method=BH --pvalue-method=norm \
    >& gatnormed_ac_${stage}_motif_AllAnnot_ChrBgd.tsv
done

# collate the three stages, adding species and stage
printf 'species\ttissue\ttrack\tannotation\tobserved\texpected\tCI95low\tCI95high\tstddev\tfold\tl2fold\tpvalue\tqvalue\ttrack_nsegments\ttrack_size\ttrack_density\tannotation_nsegments\tannotation_size\tannotation_density\toverlap_nsegments\toverlap_size\toverlap_density\tpercent_overlap_nsegments_track\tpercent_overlap_size_track\tpercent_overlap_nsegments_annotation\tpercent_overlap_size_annotation\n' > gatnormed_Collated_motif_AllAnnot_ChrBgd.tsv
for f in gatnormed_*_*_motif_AllAnnot_ChrBgd.tsv; do
  awk '{print FILENAME,$0}' OFS='\t' ${f} | sed '1,1d' | sed 's|_motif_AllAnnot_ChrBgd.tsv||g' | sed 's|gatnormed_||g' | sed 's|_|\t|' \
    >> gatnormed_Collated_motif_AllAnnot_ChrBgd.tsv
done

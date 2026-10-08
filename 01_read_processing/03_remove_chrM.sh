#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-5               # ATAC libraries only
#SBATCH --mem 48000
#SBATCH -t 0-23:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 3: remove reads mapping to mitochondrial sequence and record fragment-length counts.
# The A. calliptera mitochondrial genome (NC_018560.1) is BLASTed against the assembly; scaffolds that match
# (pident >= 93, e-value <= 1e-10, >= 75% of the mitochondrial genome covered) are removed from the BAM.
#
# Download the mitochondrial FASTA once on a node with internet access:
#   wget -O ${RESULTS}/3.chrM/${MT_ACCESSION}.fasta "https://www.ncbi.nlm.nih.gov/search/api/sequence/${MT_ACCESSION}/?report=fasta"
#
# Input : ${RESULTS}/2.aligned/<sample>.bam
# Output: ${RESULTS}/3.chrM/<sample>.nochrM.bam (+ index), <sample>.nochrM_frag_length_count.txt

source "$(dirname "$0")/../config.sh"
ml samtools/1.3
ml blast/2.3.0
ml python/3.5

SAMPLES=(${ATAC_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
OUT=${RESULTS}/3.chrM
mkdir -p "${OUT}"

# 3a. BLAST the mitochondrial genome against the assembly (database is built once)
[ -f "${GENOME_FA}.nsq" ] || makeblastdb -in "${GENOME_FA}" -parse_seqids -dbtype nucl
blastn -db "${GENOME_FA}" -outfmt 6 -evalue 1e-3 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 \
  -out "${OUT}/${GENOME_ID}.genome_mt.blast" -query "${OUT}/${MT_ACCESSION}.fasta"

# 3b. keep hits passing the identity / e-value / coverage thresholds (writes <prefix>.filtered.blast)
python "$(dirname "$0")/filter_mito_blast_hits.py" "${OUT}/${GENOME_ID}.genome_mt.blast" "${OUT}/${MT_ACCESSION}.fasta"
FILT=${OUT}/${GENOME_ID}.genome_mt.filtered.blast

# 3c. drop every read on a scaffold that matches the mitochondrial genome
IN=${RESULTS}/2.aligned/${SAMPLE}.bam
NOCHRM=${OUT}/${SAMPLE}.nochrM.bam
if [ -s "${FILT}" ] && ! grep -q "no hits" "${FILT}"; then
  scaffolds=$(cut -f2 "${FILT}" | awk '!x[$0]++')
  samtools idxstats "${IN}" | cut -f1 | grep -v -x -F "${scaffolds}" | xargs samtools view -b "${IN}" > "${NOCHRM}"
else
  cp "${IN}" "${NOCHRM}"
fi
samtools index "${NOCHRM}"

# 3d. fragment length counts of proper pairs
samtools view "${NOCHRM}" | awk '$9>0' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' \
  > "${OUT}/${SAMPLE}.nochrM_frag_length_count.txt"

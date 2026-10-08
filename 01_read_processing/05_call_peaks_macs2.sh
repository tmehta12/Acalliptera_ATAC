#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH --array=0-5
#SBATCH --mem 48000
#SBATCH -t 0-12:59
#SBATCH -o slurm.%A.%a.out
#SBATCH -e slurm.%A.%a.err

# Step 5: ATAC-seq peak calling with MACS2, using the matched naked-DNA library as control.
#   a. convert the filtered ATAC BAM and the (aligned) gDNA control BAM to tagAlign
#   b. Tn5-shift the ATAC tagAlign (+4 bp on the + strand, -5 bp on the - strand)
#   c. MACS2: -p 0.05 --nomodel --shift -75 --extsize 150, with pile-up signal tracks
#   d. fold-enrichment and log-likelihood-ratio tracks; narrowPeak -> bigBed
#
# Input : ${RESULTS}/4.filtered/<sample>.nochrM.nodup.filt.bam ; ${RESULTS}/2.aligned/<sample with _gDNA>.bam
# Output: ${RESULTS}/5.peaks/<sample>_peaks.narrowPeak(.gz/.bb), <sample>.tn5.tagAlign.gz, *_FE.bdg, *_logLR.bdg

source "$(dirname "$0")/../config.sh"
ml MACS
ml bedtools/2.25.0
ml GCC
ml zlib
ml samtools/1.3
source ucsc_utils-v333
ulimit -Sn 10000

SAMPLES=(${ATAC_SAMPLES})
SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
OUT=${RESULTS}/5.peaks
mkdir -p "${OUT}"
cd "${OUT}"

TEST_BAM=${RESULTS}/4.filtered/${SAMPLE}.nochrM.nodup.filt.bam
CONTROL_BAM=${RESULTS}/2.aligned/${SAMPLE/_ATAC/_gDNA}.bam
TAG_TEST=${SAMPLE}.nochrM.nodup.filt.tagAlign.gz
TAG_CONTROL=${SAMPLE/_ATAC/_gDNA}.tagAlign.gz
TAG_SHIFT=${SAMPLE}.tn5.tagAlign.gz

# genome size and chromosome boundaries (MACS2 -g, bedClip)
samtools faidx "${GENOME_FA}"
GENOME_SIZE=$(awk '{x+=$2} END{print x}' "${GENOME_FA}.fai")
awk '{print $1"\t"$2}' "${GENOME_FA}.fai" > scaffold_lengths.txt

# a. tagAlign
bedtools bamtobed -i "${TEST_BAM}"    | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > "${TAG_TEST}"
bedtools bamtobed -i "${CONTROL_BAM}" | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > "${TAG_CONTROL}"

# b. Tn5 shift
zcat "${TAG_TEST}" | awk -F $'\t' 'BEGIN{OFS=FS}{ if ($6=="+") {$2=$2+4} else if ($6=="-") {$3=$3-5} print $0}' | gzip -c > "${TAG_SHIFT}"

# c. MACS2 (--shift -75 --extsize 150: cut sites are extended to nucleosome-sized fragments)
macs2 callpeak -t "${TAG_SHIFT}" -c "${TAG_CONTROL}" -f BED -n "${SAMPLE}" -g "${GENOME_SIZE}" \
  -p ${MACS2_P} --nomodel --shift -75 --extsize 150 -B --SPMR --keep-dup all --call-summits

# d. signal tracks: fold enrichment and log-likelihood ratio over the control
macs2 bdgcmp -t "${SAMPLE}_treat_pileup.bdg" -c "${SAMPLE}_control_lambda.bdg" -o "${SAMPLE}_FE.bdg"    -m FE
macs2 bdgcmp -t "${SAMPLE}_treat_pileup.bdg" -c "${SAMPLE}_control_lambda.bdg" -o "${SAMPLE}_logLR.bdg" -m logLR -p 0.00001

# narrowPeak -> bigBed (UCSC tools v333)
cat > narrowPeak.as <<'AS'
table narrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
	string chrom;        "Reference sequence chromosome or scaffold"
	uint   chromStart;   "Start position in chromosome"
	uint   chromEnd;     "End position in chromosome"
	string name;         "Name given to a region (preferably unique). Use . if no name is assigned"
	uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
	char[1]  strand;     "+ or - or . for unknown"
	float  signalValue;  "Measurement of average enrichment for the region"
	float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
	float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
	int   peak;          "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
AS
sort -k1,1 -k2,2n "${SAMPLE}_peaks.narrowPeak" > "${SAMPLE}_peaks.tmp"
bedClip "${SAMPLE}_peaks.tmp" scaffold_lengths.txt "${SAMPLE}_peaks.tmp2"
bedToBigBed -type=bed6+4 -as=narrowPeak.as "${SAMPLE}_peaks.tmp2" scaffold_lengths.txt "${SAMPLE}_peaks.narrowPeak.bb"
rm -f "${SAMPLE}_peaks.tmp" "${SAMPLE}_peaks.tmp2"
gzip -f "${SAMPLE}_peaks.narrowPeak"

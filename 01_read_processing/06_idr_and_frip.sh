#!/bin/bash -e
#SBATCH -p ei-medium
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 16000
#SBATCH -t 0-04:59
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err

# Step 6: reproducibility of the two biological replicates per stage (ENCODE ATAC-seq pipeline).
#   a. IDR (v2.0.4) between replicate a and b, ranked by p-value, on the pooled-replicate peak list
#   b. every replicate peak is flagged T (rank within the number of peaks passing IDR < 0.1) or F;
#      the flagged files are *_peaks.final.narrowPeak (column 11 = T/F)
#   c. replicate peak sets per stage are concatenated (Ac_<stage>_ATAC_peaks.final.allIDR.narrowPeak) - these 'all IDR'
#      peaks (T and F) are the peak sets used throughout the paper
#   d. fraction of reads in peaks (FRiP)
#
# IDR installation: https://github.com/kundajelab/idr (v2.0.4)
# Input : ${RESULTS}/5.peaks/<sample>_peaks.narrowPeak.gz, <sample>.tn5.tagAlign.gz
# Output: ${RESULTS}/6.idr/

source "$(dirname "$0")/../config.sh"
ml bedtools/2.25.0

IN=${RESULTS}/5.peaks
OUT=${RESULTS}/6.idr
mkdir -p "${OUT}" && cd "${OUT}"

IDR_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')   # global IDR column is -log10(IDR) scaled; 0.1 -> 1

printf "sample\ttotal_peaks\tpeaks_passing_IDR-%s\tpeaks_failing_IDR-%s\tpercentage_peaks_pass_IDR-%s\n" ${IDR_THRESH} ${IDR_THRESH} ${IDR_THRESH} > IDR_summary.txt
: > FRiP.txt

for stage in ${STAGES}; do
  case ${stage} in 3dpf) n=1;; 7dpf) n=2;; 12dpf) n=3;; esac
  r1=${n}aAc_${stage}_ATAC
  r2=${n}bAc_${stage}_ATAC
  pair=${r1}_${r2}

  gunzip -c "${IN}/${r1}_peaks.narrowPeak.gz" > "${r1}_peaks.narrowPeak"
  gunzip -c "${IN}/${r2}_peaks.narrowPeak.gz" > "${r2}_peaks.narrowPeak"
  cat "${r1}_peaks.narrowPeak" "${r2}_peaks.narrowPeak" > "${pair}_pooled_peaks.narrowPeak"

  # a. IDR
  idr --samples "${r1}_peaks.narrowPeak" "${r2}_peaks.narrowPeak" --peak-list "${pair}_pooled_peaks.narrowPeak" \
      --input-file-type narrowPeak --output-file "${pair}.IDR${IDR_THRESH}output" --rank p.value \
      --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR

  # b. number of pooled peaks passing the threshold, then T/F flag on each replicate (ranked by signal value, column 7)
  npass=$(awk -v t="${IDR_TRANSFORMED}" '$12>=t' "${pair}.IDR${IDR_THRESH}output" | sort -u | wc -l)
  for rep in ${r1} ${r2}; do
    sort -k7,7rn "${rep}_peaks.narrowPeak" \
      | awk -v n="${npass}" 'BEGIN{OFS="\t"}{if(NR<=n) print $0,"T"; else print $0,"F"}' > "${rep}_peaks.final.narrowPeak"

    total=$(wc -l < "${rep}_peaks.final.narrowPeak")
    pass=$(awk '$11=="T"' "${rep}_peaks.final.narrowPeak" | wc -l)
    fail=$(awk '$11=="F"' "${rep}_peaks.final.narrowPeak" | wc -l)
    perc=$(awk -v p=${pass} -v t=${total} 'BEGIN{print p/t*100}')
    printf "%s\t%s\t%s\t%s\t%s\n" "${rep}" "${total}" "${pass}" "${fail}" "${perc}" >> IDR_summary.txt

    # d. FRiP: >0.3 pass, >0.2 acceptable
    inpeaks=$(bedtools intersect -a <(zcat -f "${IN}/${rep}.tn5.tagAlign.gz") -b "${rep}_peaks.final.narrowPeak" -wa -u | wc -l)
    reads=$(zcat "${IN}/${rep}.tn5.tagAlign.gz" | wc -l)
    awk -v s=${rep} -v a=${inpeaks} -v b=${reads} 'BEGIN{f=a/b; print s"\t"a"/"b"\t"f"\t"(f>=0.3?"FRiP PASS":(f>=0.2?"FRiP PASS - acceptable":"FRiP FAIL"))}' >> FRiP.txt
  done

  # c. both replicates (T and F) per stage
  cat "${r1}_peaks.final.narrowPeak" "${r2}_peaks.final.narrowPeak" | sort -k1,1 -k2,2n > "Ac_${stage}_ATAC_peaks.final.allIDR.narrowPeak"
done

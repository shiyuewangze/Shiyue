#call peak
macs2 callpeak -c INPUT_hisat2_MAPQ20.sorted.bam -t IP_hisat2_MAPQ20.sorted.bam  -n sample-1 -f BAM  -g hs  -q 0.01  --nomodel --keep-dup all --outdir /peak
grep -v "#"  sample-1_peaks.xls | grep -v start | awk '{if($7>5) print $1"\t"$2"\t"$3"\t"$8}'  > sample-1_peaks_high-confidence.bed
intersectBed -a sample-1_peaks_high-confidence.bed -b sample-2_peaks_high-confidence.bed > sample12_overlappeak.bed 

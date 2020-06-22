#!/usr/bin/env bash

export PERL5LIB=/root/perl5/lib/perl5/x86_64-linux-thread-multi:/usr/share/perl5/czplib-1.0.7:/usr/local/lib/czplib/czplib:$PERL5LIB



IP1="samp1IP_hisat2_MAPQ20.sorted.bam"
INPUT1="samp1INPUT_hisat2_MAPQ20.sorted.bam"

#bedgraph
/data/software/deepTools-2.5.7/bin/bamCompare -b1 $IP1 -b2 $INPUT1 -o subtractsample1.bw --operation subtract  --outFileFormat bigwig --binSize 10 

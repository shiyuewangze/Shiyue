#!/usr/bin/env bash

export PERL5LIB=/root/perl5/lib/perl5/x86_64-linux-thread-multi:/usr/share/perl5/czplib-1.0.7:/usr/local/lib/czplib/czplib:$PERL5LIB

bin3=/data/software/samtools-1.8
hisat2=/data/software/hisat2-2.1.0

java=/data/software/jdk1.8.0_181/bin/java

/data/software/fastx_toolkit-0.0.14/src/fastx_clipper/fastx_clipper -a AGATCGGAAGAGCACACG -l 24 -n -i RIP_raw1.fastq -Q 33 -o RIP_trim_1.fastq
/data/software/FastQC/fastqc RIP_trim_1.fastq  -t 4 -o ./fastqc

cutadapt -a "A{8}"  -o RIP_qua_1.fastq  RIP_trim_1.fastq
/data/software/FastQC/fastqc RIP_qua_1.fastq  -t 4 -o ./fastqc

awk '{if(NR%4==1){split($1,a,"#");print a[1]"#"a[2]}else{print $0}}' RIP_qua_1.fastq > temp1.fastq

perl /data/software/ctk-master/stripBarcode.pl -format fastq -len 3 temp1.fastq debarcode_RIP_qua_1.fastq
/data/software/FastQC/fastqc debarcode_RIP_qua_1.fastq -t 4 -o ./fastqc


java -Xmx4g -jar /data/software/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33  debarcode_RIP_qua_1.fastq qf_trim_RIP.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:36
/data/software/FastQC/fastqc qf_trim_RIP.fastq -t 4 -o ./fastqc

$hisat2/hisat2 -x /data/yuzh/hisat2-build/genome_tran -p 20 -U qf_trim_RIP.fastq  -S RIP_hisat2.sam 

/data/software/samtools-1.8/samtools view RIP_hisat2.sam -q 20 -h -S > RIP_hisat2_q20.sam
/data/software/samtools-1.8/samtools view -S RIP_hisat2_q20.sam -b -o RIP_hisat2_q20.bam
/data/software/samtools-1.8/samtools sort RIP_hisat2_q20.bam -o  RIP_hisat2_q20.sorted.bam

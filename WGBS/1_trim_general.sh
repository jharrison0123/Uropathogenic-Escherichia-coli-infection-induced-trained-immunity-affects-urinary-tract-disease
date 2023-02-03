#!/usr/bin/env bash

#cutadapt                  1.18 
#fastqc                    0.11.8  

dir=/opt/home/jharrison/scripts
date

### INPUT
fq1=${1?Error: no fastq given}
echo $fq1
threads=${2?Error: no thread given}
echo $threads

fq2=${fq1/R1/R2}
sample=(${fq1%.f*})
echo $sample

cd ..
basedir=$(pwd)
echo $basedir
mkdir -p 1_trim
cd 1_trim
indir=${basedir}
outdir=${basedir}/1_trim

read1_fq=${indir}/${fq1}
read2_fq=${indir}/${fq2}

echo $read1_fq
echo $read2_fq


### OUTPUT
log_cutadapt=${sample}.cutadapt.log
trim1_fq=${outdir}/trimmed_${sample}_R1.fq.gz
trim2_fq=${outdir}/trimmed_${sample}_R2.fq.gz
trimmed_sample=trimmed_${sample}

### COMMANDS

cutadapt -a "AGATCGGAAGAGC" -A "AGATCGGAAGAGC" -q 15,10 -u 10 -U 15 --minimum-length 36  -o ${trim1_fq} -p ${trim2_fq}  \
    ${read1_fq} ${read2_fq} > ${outdir}/${log_cutadapt}



fastqc -t $threads -o ${outdir} --noextract --nogroup  ${trim1_fq} &> ${outdir}/${trimmed_sample}.fastqc.log

fastqc -t $threads -o ${outdir} --noextract --nogroup ${trim2_fq} &> ${outdir}/${trimmed_sample}.fastqc.log

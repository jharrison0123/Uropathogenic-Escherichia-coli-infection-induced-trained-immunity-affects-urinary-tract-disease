#!/usr/bin/env bash

#fastqc                    0.11.8 

### INPUT
fq1=${1?Error: no fastq given}
echo $fq1
threads=${2?Error: no thread given}
echo $threads

fq2=${fq1/R1/R2}
sample=(${fq1%.f*})
echo $sample


### OUTPUT

basedir=$(pwd)
echo $basedir
mkdir -p 0_fastqc
cd 0_fastqc
indir=${basedir}
outdir=${basedir}/0_fastqc

fastqc -t $threads -o $outdir --noextract --nogroup ${indir}/${fq1} &> ${outdir}/${i}.fastqc.log

fastqc -t $threads -o $outdir --noextract --nogroup ${indir}/${fq2} &> ${outdir}/${i}.fastqc.log

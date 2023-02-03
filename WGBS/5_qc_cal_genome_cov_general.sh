#!/bin/bash 

### INPUT
fq1=${1?Error: no fastq given}
echo $fq1
threads=${2?Error: no thread given}
echo $threads


fq2=${fq1/R1/R2}
presample=(${fq1%.f*})
sample=trimmed_${presample}
echo $sample


cd ..
basedir=$(pwd)
echo $basedir
mkdir -p 5_coverage
cd 5_coverage
indir=${basedir}/2_bismark
dir_out=${basedir}/5_coverage


### INPUT

cx_me=${indir}/${sample}_R1_R1_bismark_bt2.CXme.txt 
echo $cx_me

cnt_c=$(( 598683433+600854940 ))
echo $cnt_c # Watson strand + Crick strand	
cnt_c=$(( $cnt_c - 171823*2 )) 		# Discard chrEBV
echo $cnt_c
cnt_cg=$(( 29303965 * 2 ))
echo $cnt_cg

# OUTPUT
cov_genome=${dir_out}/${sample}.genome_cov.txt


# COMMANDS
c_cov=$( cat $cx_me | awk -F"\t" -v c=$cnt_c 'BEGIN{s=0} {s+=$2+$3} END{print s/c}' )
cg_cov=$( cat $cx_me | awk -F"\t" -v c=$cnt_cg 'BEGIN{s=0} $1=="CG" {s+=$2+$3} END{print s/c}' )

echo -e "$sample\t$c_cov\t$cg_cov" > $cov_genome

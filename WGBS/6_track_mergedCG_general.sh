#!/bin/bash

#INPUT
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
mkdir -p 6_track_mergedCG
cd 6_track_mergedCG
indir=${basedir}/2_bismark
dir_out=${basedir}/6_track_mergedCG

cx_report=${indir}/${sample}_R1_R1.CX_report.txt.gz 
# OUTPUT
methylc=${dir_out}/${sample}.CG.methylC.gz

# COMMNADS for strand-merged CG track
zcat $cx_report | awk -F"\t" 'BEGIN{OFS=FS} $6=="CG" && $4+$5>0 { if ($3=="+") {print $1,$2-1,$2+1,$4,$5} if ($3=="-") {print $1,$2-2,$2,$4,$5} }' |
	        sort -k1,1 -k2,2n |
		    groupBy -g 1,2,3 -c 4,5 -o sum,sum |
		        awk -F"\t" 'BEGIN{OFS=FS} {mcg=sprintf("%.3f", $4/($4+$5)); print $1,$2,$3,"CG",mcg,"+",$4+$5 }' |
			    bgzip >$methylc
tabix -p bed $methylc

#!/usr/bin/env bash

min_insert=0
max_insert=2000


date

### INPUT
fq1=${1?Error: no fastq given}
echo $fq1
threads=${2?Error: no thread given}
echo $threads
genome_dir=${3?Error: no genome_dir}

fq2=${fq1/R1/R2}
presample=(${fq1%.f*})
sample=trimmed_${presample}
echo $sample



### OUTPUT
cd ..
basedir=$(pwd)
echo $basedir
mkdir -p 2_bismark_lambda
cd 2_bismark_lambda
indir=${basedir}/1_trim
outdir=${basedir}/2_bismark_lambda

uncomp_read1_fq=${indir}/trimmed_${presample}_R1.fq.gz
uncomp_read2_fq=${indir}/trimmed_${presample}_R2.fq.gz
echo $uncomp_read1_fq
### OUTPUT

# 1. bismark PE files
log_bismark_pe=${outdir}/${sample}.bismark.log
bam_pe=${outdir}/${sample}_R1_bismark_bt2_pe.bam
read1_unmapped_fq=${outdir}/${sample}_R1.fq.gz_unmapped_reads_1.fq.gz
read2_unmapped_fq=${outdir}/${sample}_R2.fq.gz_unmapped_reads_2.fq.gz

# 2. deduplicate files
log_dedup_pe=${outdir}/${sample}_R1.deduplicate_bismark.log
bam_dedup_pe=${outdir}/${sample}_R1_bismark_bt2_pe.deduplicated.bam

# 3. extract files
log_methx_pe=${outdir}/${sample}_R1_pe.bismark_methylation_extractor.log

cg_pe=${outdir}/CpG_context_${sample}_R1_bismark_bt2_pe.deduplicated.txt.gz
ch_pe=${outdir}/Non_CpG_context_${sample}_R1_bismark_bt2_pe.deduplicated.txt.gz
merged=${outdir}/${sample}_R1_bismark_bt2.extracted.txt.gz

# 4. bismark2bedgraph files
log_bismark2bg=${outdir}/${sample}_R1.bismark2bedGraph.log
bedGraph=${sample}_R1.bedGraph.gz
cov=${sample}_R1.bismark.cov.gz

# 5. coverage2cytosine files
log_cov2c=${outdir}/${sample}_R1.coverage2cytosine.log
cx_report=${sample}
cx_me=${outdir}/${sample}_R1_bismark_bt2.CXme.txt


### COMMANDS
echo "-- Started on $(date) "
echo ""
echo ""



# Mapping with bismark/bowtie2
 #Note --bowtie2 and -p $nthreads are both SLOWER than single threaded bowtie1
echo "-- 1. Mapping to reference with bismark/bowtie2... started on $(date)"
bismark -q -I $min_insert -X $max_insert --parallel 2 -p 2 \
            --bowtie2 --score_min L,0,-0.6 -N 0 \
           -o $outdir --nucleotide_coverage \
            $genome_dir -1 $uncomp_read1_fq -2 $uncomp_read2_fq &>$log_bismark_pe

echo ""

# Dedpulicate reads
echo "-- 2. Deduplicating aligned reads... started on $(date)"
deduplicate_bismark -p --bam $bam_pe  &>$log_dedup_pe
echo ""

# Run methylation extractor for the sample
echo "-- 3. Analyse methylation in $bam_dedup_pe using $threads threads... started on $(date)"
bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
                   -o $outdir --gzip --parallel $threads \
                                  $bam_dedup_pe &>$log_methx_pe
echo ""

# Generate HTML Processing Report
echo "-- 4. Generate bismark HTML processing report file... started on $(date)"
bismark2report -o ${sample}_R1_bismark_bt2_PE_report.html --dir $outdir \
          --alignment_report ${outdir}/${sample}_R1_bismark_bt2_PE_report.txt \
         --dedup_report ${outdir}/${sample}_R1_bismark_bt2_pe.deduplication_report.txt \
          --splitting_report ${outdir}/${sample}_R1_bismark_bt2_pe.deduplicated_splitting_report.txt \
          --mbias_report ${outdir}/${sample}_R1_bismark_bt2_pe.deduplicated.M-bias.txt \
          --nucleotide_report ${outdir}/${sample}_R1_bismark_bt2_pe.nucleotide_stats.txt &>${outdir}/${sample}_R1_pe.bismark2report.log
echo ""

# Generate bedGraph file
echo "-- 5. Generate bedGraph file... started on $(date)"
mv $ch_pe $merged
cat $cg_pe >>$merged
bismark2bedGraph --dir $outdir --cutoff 1 --CX_context --buffer_size=75G --scaffolds \
                 -o $bedGraph $merged &>$log_bismark2bg
#rm ${outdir}/$bedGraph # $merged
echo ""

# Calculate average methylation levels per each CN context
echo "-- 6. Generate cytosine methylation file... started on $(date)"
coverage2cytosine -o $cx_report --dir $outdir --genome_folder $genome_dir --CX_context --gzip \
                  $cov &>$log_cov2c
zcat ${outdir}/${cx_report}.CX_report.txt.gz  |
    awk 'BEGIN{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}
         $7~/^CA/ {ca+=$5; mca+=$4}
         $7~/^CC/ {cc+=$5; mcc+=$4}
         $7~/^CG/ {cg+=$5; mcg+=$4}
         $7~/^CT/{ct+=$5; mct+=$4}
         END{printf("CA\t%d\t%d\t%.3f\n", ca, mca, mca/(ca+mca));
             printf("CC\t%d\t%d\t%.3f\n", cc, mcc, mcc/(cc+mcc));
             printf("CG\t%d\t%d\t%.3f\n", cg, mcg, mcg/(cg+mcg));
             printf("CT\t%d\t%d\t%.3f\n", ct, mct, mct/(ct+mct));}' >$cx_me
echo ""

# Print the files generated
echo "-- The results..."
ls -l ${outdir}/*${sample}*
echo ""
echo "-- Finished on $(date)"

##trim adapters and low quality reads
#cutadapt 1.9
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --quality-cutoff=15,10 --minimum-length=36 -u 10 -U 10 -o Trimmed_"$xbase" -p Trimmed_"${xbase/R1/R2}" "$file" "${file/R1/R2}" > "$xbase"_cutadapt.log

##MACS2 peak calling
MACS2 peaks were called using macs2 callpeak -f BEDPE --keep dup all, with treatment and control files. 
For H3K27me3, the --broad flag was added
https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4#Sec11

##H3K27Ac, H3K4me4
#macs2 2.1.1.20160309
find . -name "*_qffilter.bam" | while read file ; do if [[ $file != Trimmed_input-cutrun_R1.fastq_nochrM_nodup_qffilter.bam ]]; then echo "macs2 callpeak -t "$file" -c Trimmed_"${INPUT/.gz/_nochrM_nodup_qffilter.bam}" -f BAMPE -B -n "${file%.*}"_uniqpeak -q 0.01 --keep-dup all" >> 10_macs2Call_uniqpeak.txt; fi ; done ;

## H3K27me3
#macs2 2.1.1.20160309
find . -name "*_qffilter.bam" | while read file ; do if [[ $file != Trimmed_input-cutrun_R1.fastq_nochrM_nodup_qffilter.bam ]]; then echo "macs2 callpeak -t "$file" -c Trimmed_input-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -f BAMPE -B -n "${file%.*}" -q 0.05 --broad --keep-dup all " >> 10_macs2Call_broad_uniqpeak.txt ;fi ; done ;

#combine replicates
#samtools                  1.6 
samtools merge N_H3K27Ac_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_N1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_N3-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge N_H3K27Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_N1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_N3-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge N_H3K4Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_N1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_N3-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge R_H3K27Ac_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_R_1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_R4-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge R_H3K27Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_R_1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_R4-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge R_H3K4Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_R_1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_R4-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge S_H3K27Ac_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_S1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_S2-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge S_H3K27Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_S1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_S2-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam
samtools merge S_H3K4Me3_combined_bismark_bt2_pe.deduplicated.bam -@ 14 Trimmed_WangT_S1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam Trimmed_WangT_S2-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam

# 1/FRIP and read coverage normalization for bigwig signal (FRIP scores caluclated using consensus peak set from DiffBind-see DiffBind.R script for commands
#bamCoverage 3.4.2
#ATAC
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Nai1.filt.srt.nodup.bam -o Hultgren_Nai1.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.348 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Nai3.filt.srt.nodup.bam -o Hultgren_Nai3.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 5.263 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Res1.filt.srt.nodup.bam -o Hultgren_Res1.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 3.571 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Res4.filt.srt.nodup.bam -o Hultgren_Res4.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 3.333 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Sen1.filt.srt.nodup.bam -o Hultgren_Sen1.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 3.704 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam /scratch/jharrison/Seongmi/ATAC/processing/3_filtering/Hultgren_Sen2.filt.srt.nodup.bam -o Hultgren_Sen2.filt.srt.nodup_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.762 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads

#H3K4me3
bamCoverage --bam Trimmed_WangT_N1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_N1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 2.222 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_N3-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_N3-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 2.381 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R_1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_R1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 5.263 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R4-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_R4-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 3.448 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_S1-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 2.381 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S2-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_S2-H3K4Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw  --scaleFactor 7.143 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads

#H3K27Ac
bamCoverage --bam Trimmed_WangT_N1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_N1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.348 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_N3-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_N3-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.000 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R_1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -oTrimmed_WangT_R1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 5.882 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R4-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_R4-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.762 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_S1-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 4.348 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S2-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_S2-H3K27Ac-cutrun_R1.fastq_nochrM_nodup_qffilter_narrow_FRIP_inverse_deseq.bw --scaleFactor 5.000 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads

#H3K27Me3
bamCoverage --bam Trimmed_WangT_N1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_N1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_broad_FRIP_inverse_deseq.bw --scaleFactor 3.030 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_N3-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_N3-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_broad_FRIP_inverse_deseq.bw --scaleFactor 2.500 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R_1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam -o Trimmed_WangT_R1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilte_broad_FRIP_inverse_deseq.bw --scaleFactor 2.632 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_R4-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_R4-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_broad_FRIP_inverse_deseq.bw --scaleFactor 2.564 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_S1-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_broad_FRIP_inverse_deseq.bw --scaleFactor 2.439 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads
bamCoverage --bam Trimmed_WangT_S2-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter.bam  -o Trimmed_WangT_S2-H3K27Me3-cutrun_R1.fastq_nochrM_nodup_qffilter_broad_FRIP_inverse_deseq.bw --scaleFactor 3.030 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --binSize 10 --extendReads

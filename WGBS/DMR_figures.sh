#create bigwig files using CGs found in combined DMR comparisons

#create combined DMR list from two group comparisons using pvalue .001 

zcat *p001_twogroup.txt.gz | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3}' | sort -k1,1 -k2,2n | sed 1,3d | mergeBed -i stdin > combined_DMRs.bed

#reformat DMR two comparison group pvalue 0.001 files
zcat DMRs_Naive_1_Naive_3vsSensitized_1_Sensitized_2_p001_twogroup.txt.gz | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3}' | sort -k1,1 -k2,2n | sed 1,1d  > N_S_DMRs.bed #naive vs sensitized
zcat DMRs_Resolved_1_Resolved_4vsSensitized_1_Sensitized_2_p001_twogroup.txt.gz | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3}' | sort -k1,1 -k2,2n | sed 1,1d > R_S_DMRs.bed #resolved vs sensitized
zcat DMRs_Naive_1_Naive_3vsResolved_1_Resolved_4_p001_twogroup.txt.gz | awk 'BEGIN {OFS="\t"}; {print $1, $2, $3}' | sort -k1,1 -k2,2n | sed 1,1d  > N_R_DMRs.bed #naive vs resolved

#identify overlap between naive vs sensitized vs resolved vs sensitized DMR comparisons
bedtools intersect -a N_S_DMRs.bed -b R_S_DMRs.bed > S_specific_intersection_DMRs.bed

ml kentUCSC/334

#interesect bedgraph of CGs with combined DMR regions

find *deduplicated.bedGraph | while read file ; do base=${file%R1%deduplicated.bedGraph} ; cat $file | sort -k1,1 -k2,2n | bedtools intersect -a combined_DMRs.bed -b stdin -wb -sorted > ${base}_DMR_CG.bedGraph ; done ; \

#convert CG bedgraph within combined DMR regions to bigwig file
find *DMR_CG.bedGraph | while read file ; do base=${file%.bedGraph} ; cat $file | awk 'BEGIN{OFS="\t"} {print $4, $5, $6, $7}' | bedtools sort -i stdin > ${base}_sorted.bedGraph ; done
bedGraphToBigWig ${base}_sorted.bedGraph /scratch/jharrison/scratch_reference/mouse/mm10.chrom.sizes ${base}_sorted.bw ; done

#create PCA plot of bigwig signal of CGs in DMRs using Deeptools 3.3.0

multiBigwigSummary bins --outFileName /scratch/jharrison/Seongmi/WGBS/May_DMR/analysis/DMR_CG_multiBigwigsummary_1kb.npz  --bwfiles  \
trimmed_HHYVGDSX2_ATTCAGAAAT-ACAGGCGAAG_L003_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
trimmed_HTLTYDSX2_ATTCAGAAAT-ACATAGAGGC_L003_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
trimmed_HTLTYDSX2_ATTCAGAAAT-ACCCTATCCT_L001_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
trimmed_HHYVGDSX2_ATTCAGAAAT-ACTAATCTTA_L002_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
trimmed_HTLTYDSX2_ATTCAGAAAT-ACTATAGCCT_L004_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
trimmed_HHYVGDSX2_ATTCAGAAAT-ACGGCTCTGA_L004_R1_R1_bismark_bt2_pe.deduplicated.bedGraph_DMR_CG_sorted.bw \
 --labels N1 N3 R1 R4 S1 S2 --binSize 1000 -p 12 ; \
 
 plotPCA --transpose -in DMR_CG_multiBigwigsummary_1kb.npz  --markers o --colors blue blue orange orange purple purple -o DMR_CG_pca_1kb.pdf  --outFileNameData DMR_CG_pca_1kb \
 --plotHeight 15 --plotWidth 10 
 
 ##Intervene graph
 ###Intervene venn 

intervene venn -i \
N_S_DMRs.bed \ #naive vs sensitized
R_S_DMRs.bed \ #resolved vs sensitized
N_R_DMRs.bed \ #naive vs resolved
--output DMR_intersect -–save-overlaps


###DMR annotation of Sensitized-specific hypo-DMRs extracted from complexheatmap cluster1- hypo DMRs in sensitized samples
bedtools intersect \
-a DMRs_S_specific_cluster1_annotation.bed \
-b /scratch/jharrison/scratch_reference/mouse/Gencode.M25.mm10/gencode_M25_JH_genic_annotation.bed -f 0.2 -loj -wa -wb > DMRs_S_specific_cluster1_annotation.bed.PeakGenic_intersection ; \
groupBy -i DMRs_S_specific_cluster1_annotation.bed.PeakGenicIntersection -g 4 -c 11 -o distinct > DMRs_S_specific_cluster1_annotation.bed_distinct.PeakGenicIntersection ; \
cat DMRs_S_specific_cluster1_annotation.bed_distinct.PeakGenicIntersection | awk '$2 ~ /UpstreamBy1kb/ { print }' | awk '{print $1, "UpstreamBy1kb" }' > UpstreamBy1kb.tmp 
cat DMRs_S_specific_cluster1_annotation.bed_distinct.PeakGenicIntersection | awk '$2 !~ /UpstreamBy1kb/ { print }' > tmp1.txt
cat tmp1.txt | awk '$2 ~ /coding_exons/ { print }' | awk '{print $1, tmux a Exon" }' > coding_exons.tmp
cat tmp1.txt | awk '$2 !~ /coding_exons/ { print }' > tmp2.txt
cat tmp2.txt | awk '$2 ~ /5UTR/ { print }' | awk '{print $1, "5UTR" }'  > 5UTR.tmp
cat tmp2.txt | awk '$2 !~ /5UTR/ { print }' > tmp3.txt
cat tmp3.txt | awk '$2 ~ /3UTR/ { print }'| awk '{print $1, "3UTR" }'  > 3UTR.tmp
cat tmp3.txt | awk '$2 !~ /3UTR/ { print }' > tmp4.txt
cat tmp4.txt | awk '$2 ~ /Intron/ { print }' | awk '{print $1, "Intron" }'  > intron.tmp
cat tmp4.txt | awk '$2 !~ /Intron/ { print }' > tmp5.txt
cat tmp5.txt | awk '$2 ~ /Intergenic/ { print }' | awk '{print $1, "Intergenic" }' > intergenic.tmp
cat tmp5.txt | awk '$2 !~ /Intergenic/ { print }' > tmp6.txt ; \
cat UpstreamBy1kb.tmp 5UTR.tmp 3UTR.tmp coding_exons.tmp 5UTR.tmp 3UTR.tmp intron.tmp intergenic.tmp | awk 'BEGIN{FS=":"} {print $1, $2, $5}' \
| awk 'BEGIN{FS="-"} {print $1, $2, $3, $4}' > DMRs_S_specific_cluster1_annotation.bed.PeakGenicIntersection

### DNA methylation (CG) using over Sensitized-specific hypo-DMRs using Deeptools 3.3.0 and smoothed with ggplot2 smooth function- see R commands

file=DMRs_S_specific_cluster1_annotation.bed

computeMatrix scale-regions -S \
/scratch/jharrison/Seongmi/WGBS/2_bismark/N1_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
/scratch/jharrison/Seongmi/WGBS/2_bismark/N3_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
/scratch/jharrison/Seongmi/WGBS/2_bismark/R1_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
/scratch/jharrison/Seongmi/WGBS/2_bismark/R4_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
/scratch/jharrison/Seongmi/WGBS/2_bismark/S1_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
/scratch/jharrison/Seongmi/WGBS/2_bismark/S2_deduplicated.bismark_bt2_pe.deduplicated.bedGraph.gz_sorted.bw \
-R $file \
-p 20 --beforeRegionStartLength 5000 \
--regionBodyLength 5000 \
--afterRegionStartLength 5000 \
--skipZeros -o WGBS_DMR_${base}_clustering_matrix.mat.gz \
--outFileNameMatrix WGBS_DMR_${base}_clustering_matrix.tsv ; \

plotProfile -m WGBS_DMR_${base}_clustering_matrix.mat.gz \
--perGroup \
--colors blue mediumBlue gold khaki mediumpurple slateblue \
-out WGBS_DMR_${base}_clustering_profile.png  \
--samplesLabel N1  N3 R1 R4 S1 S2 \
--outFileNameData WGBS_DMR_${base}_clustering_profile.tsv


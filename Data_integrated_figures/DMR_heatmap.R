#R 3.6.1
# ComplexHeatmap 2.11.2
# rtracklayer 1.44.4
# dplyr 1.0.9 
# devtools)
# ggplot2)
# GenomicRanges)
# GenomicFeatures)


library(ComplexHeatmap)
library(rtracklayer)
library(dplyr)
library(devtools)
library(ggplot2)
library(GenomicRanges)
library(GenomicFeatures)

### DMR intersection with ATAC, cutandrun, DNA methylation and RNA seq

## Read in DMR regions and bigwig files to get methylation signal in DMRs
dmr.mat.wgbs <- function(bw.vec, dmr.bed, threads = 1) {
  
  dmr.gr <- rtracklayer::import.bed(dmr.bed)
  dmr.gr$dmr.locus <- paste0(seqnames(dmr.gr), ":", start(dmr.gr), "-", end(dmr.gr))
  df <- utilsFanc::safelapply(names(bw.vec), function(sample) {
    bw <- bw.vec[[sample]]
    data.in.dmr <- rtracklayer::import(con = bw)
    j <- plyranges::join_overlap_left(dmr.gr, data.in.dmr)
    
    df <- mcols(j) %>% as.data.frame() %>% 
      dplyr::group_by(dmr.locus) %>% dplyr::summarise(score= mean(score)) %>% 
      dplyr::ungroup() %>% as.data.frame()
    colnames(df) <- c("dmr.locus", sample)
    return(df)
  }, threads = threads) %>% Reduce(left_join, .)
   df <- df %>% na.omit()
  #browser()
  rownames(df) <- df$dmr.locus
  df$dmr.locus <- NULL
  mat <- as.matrix(df)
  return(mat)
}


## Read in narrowPeaks that overlap with DMR regions to calculate peakscores over DMR regions- adjust peakscore for 1/FRIP score

dmr.norm.mat.diff.locus <- function(peak.bed.vec, peak.bed.df, dmr.gr, threads =10 ) {
  

  df <- utilsFanc::safelapply(names(peak.bed.vec), function(sample) {
    browser()
    peak.bed <- peak.bed.vec[[sample]]
    #import narrowPeak file
    data.in.diffregion.gr <- rtracklayer::import(peak.bed, format="narrowPeak")
    #create diff locus = narrow peak position
    data.in.diffregion.gr$diff.locus <- paste0(seqnames(data.in.diffregion.gr), ":", start(data.in.diffregion.gr), "-", end(data.in.diffregion.gr))
    #combine sub peaks with same position into one peak score
    data.in.diffregion.df <- data.in.diffregion.gr %>% as.data.frame() %>% dplyr::group_by(diff.locus) %>% dplyr::summarise(area= sum(score)) %>% dplyr::ungroup() %>% as.data.frame() %>% dplyr::select(diff.locus, area)  %>% separate(diff.locus, c("chr","tmp"), ":") %>% separate(tmp, c("start","stop"), "-")
    #add diff locus
    data.in.diffregion.df$diff.locus <- paste0((data.in.diffregion.df$chr), ":", (data.in.diffregion.df$start), "-", (data.in.diffregion.df$stop))
    #put diff region into granges
      data.in.diffregion.gr <- makeGRangesFromDataFrame(data.in.diffregion.df, keep.extra.columns=T)
     #data.in narrow peak overlap with DMR regions (by 1bp)-keep metadata of diff locus score and region
    DiffLocus_DMR.df <- plyranges::join_overlap_left(dmr.gr,data.in.diffregion.gr) %>% as.data.frame() 
    colnames(DiffLocus_DMR.df) <- c("chr", "start","stop", "diff", "strand", "dmr.locus" , "score","diff.locus")
    DiffLocus_DMR.df[is.na(DiffLocus_DMR.df)] <- 0
    #sumarize diff locus score by DMR region (in the case of two diff locus overlapping one DMR)
    df <-  DiffLocus_DMR.df %>% dplyr::group_by(dmr.locus) %>% dplyr::summarise(area=mean(score)) %>% dplyr::ungroup() %>% as.data.frame() %>% dplyr::select(dmr.locus, area)
  #browser()
    colnames(df) <- c("dmr.locus", sample)
    return(df)
  }, threads = threads) %>% Reduce(left_join, .)
  df[is.na(df)] <- 0
  #browser()
  rownames(df) <- df$dmr.locus
  df$dmr.locus <- NULL
  mat <- as.matrix(df)
  peak.bed.df$inverse_FRIP <- (1/peak.bed.df$FRIP)
  norm_mat <- mat %*% diag(c(peak.bed.df$inverse_FRIP))
  colnames(norm_mat) <- colnames(mat)
  return(norm_mat)
}


## Read in broadPeaks that overlap with DMR regions to calculate peakscores over DMR regions- adjust peakscore for 1/FRIP score


dmr.broad.norm.mat.diff.locus <- function(peak.bed.vec, peak.bed.df, dmr.gr, threads =10 ) {
  

  df <- utilsFanc::safelapply(names(peak.bed.vec), function(sample) {
    #browser()
    peak.bed <- peak.bed.vec[[sample]]
    #import narrowPeak file
    data.in.diffregion.gr <- rtracklayer::import(peak.bed, format="broadPeak")
    #create diff locus = narrow peak position
    data.in.diffregion.gr$diff.locus <- paste0(seqnames(data.in.diffregion.gr), ":", start(data.in.diffregion.gr), "-", end(data.in.diffregion.gr))
    #combine sub peaks with same position into one peak score
    data.in.diffregion.df <- data.in.diffregion.gr %>% as.data.frame() %>% dplyr::group_by(diff.locus) %>% dplyr::summarise(area= sum(score)) %>% dplyr::ungroup() %>% as.data.frame() %>% dplyr::select(diff.locus, area)  %>% separate(diff.locus, c("chr","tmp"), ":") %>% separate(tmp, c("start","stop"), "-")
    #add diff locus
    data.in.diffregion.df$diff.locus <- paste0((data.in.diffregion.df$chr), ":", (data.in.diffregion.df$start), "-", (data.in.diffregion.df$stop))
    #put diff region into granges
      data.in.diffregion.gr <- makeGRangesFromDataFrame(data.in.diffregion.df, keep.extra.columns=T)
     #data.in narrow peak overlap with DMR regions (by 1bp)-keep metadata of diff locus score and region
    DiffLocus_DMR.df <- plyranges::join_overlap_left(dmr.gr,data.in.diffregion.gr) %>% as.data.frame() 
    colnames(DiffLocus_DMR.df) <- c("chr", "start","stop", "diff", "strand", "dmr.locus" , "score","diff.locus")
    DiffLocus_DMR.df[is.na(DiffLocus_DMR.df)] <- 0
    #sumarize diff locus score by DMR region (in the case of two diff locus overlapping one DMR)
    df <-  DiffLocus_DMR.df %>% dplyr::group_by(dmr.locus) %>% dplyr::summarise(area=mean(score)) %>% dplyr::ungroup() %>% as.data.frame() %>% dplyr::select(dmr.locus, area)
  #browser()
    colnames(df) <- c("dmr.locus", sample)
    return(df)
  }, threads = threads) %>% Reduce(left_join, .)
  df[is.na(df)] <- 0
  #browser()
  rownames(df) <- df$dmr.locus
  df$dmr.locus <- NULL
  mat <- as.matrix(df)
  peak.bed.df$inverse_FRIP <- (1/peak.bed.df$FRIP)
  norm_mat <- mat %*% diag(c(peak.bed.df$inverse_FRIP))
  colnames(norm_mat) <- colnames(mat)
  return(norm_mat)
}

#################  WGBS data ##########################
#######################################################


fastq.df <- read.table("Seongmi_sample_condition.txt", 
                    header = T, sep= "\t")[1:6, ]

fastq.df$root <- (fastq.df$File.Name) %>% sub("_R1.fastq.gz$", "", .)
bw.dir <- "/scratch/jharrison/Seongmi/WGBS/2_bismark/bedgraph"
fastq.df$bw <- sapply(fastq.df$root, function(x) {
  print(x)
  bw <- Sys.glob(paste0(bw.dir, "/*", "trimmed_", x, "*_R1_R1_bismark_bt2_pe.deduplicated._sorted..bw"))
  print(bw)
})

bw.vec <- fastq.df$bw
names(bw.vec) <- fastq.df$Sample.Name
dmr.bed <- "S_specific_intersection_DMRs.bed"

#create WGBS matrix of CG methylation signal within DMR regions
WGBS_mat <- dmr.mat.wgbs(bw.vec = bw.vec, dmr.bed = dmr.bed, threads = 16)

#create heatmap of S_specific DMRs
WGBS <- Heatmap(matrix = WGBS_mat, show_row_names = T, row_names_gp = gpar(fontsize = 5) ,clustering_distance_rows = "pearson", column_title = "DNA Methylation", column_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 9,fontface = "bold") , border=T, cluster_columns = FALSE, width = 4, height=5, cluster_rows = T, na_col = "black", heatmap_legend_param = list(title = "DNA methylation"))
WGBS

tiff("Seongmi_WGBS_S_specific_intersection.tiff", width = 8, height = 6, units = 'in', res = 300, compression = 'lzw')
WGBS  # Make plot
dev.off()

############  export DMR heatmap clusters ################
#########################################################

# set clustering to split rows by 2 - "row_split = 2"- get two clusters- cluster1 represents hypomethylated DMRs in sensitized-intersection
WGBS_cluster <- Heatmap(matrix = WGBS_mat, show_row_names = F, row_names_gp = gpar(fontsize = 5) ,clustering_distance_rows = "pearson", column_title = "DNA Methylation", column_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 9,fontface = "bold") , border=T, cluster_columns = FALSE, width = 4, height=5, cluster_rows = WGBS_r.dend, row_split = 2, na_col = "black", heatmap_legend_param = list(title = "DNA methylation"))
WGBS_cluster
WGBS_cluster_heatmap <- draw(WGBS_cluster)

WGBS_cluster_rcl.list <- row_order(WGBS_cluster_heatmap)

 lapply(WGBS_cluster_rcl.list, function(x) length(x))  #check/confirm size clusters

 
 # loop to extract genes for each cluster.
for (i in 1:length(row_order(WGBS_cluster_heatmap))){
 if (i == 1) {
  clu <- t(t(row.names(WGBS_mat[row_order(WGBS_cluster_heatmap)[[i]],])))
 out <- cbind(clu, paste("cluster", i, sep=""))
 colnames(out) <- c("DMR", "Cluster")
 } else {
 clu <- t(t(row.names(WGBS_mat[row_order(WGBS_cluster_heatmap)[[i]],])))
 clu <- cbind(clu, paste("cluster", i, sep=""))
 out <- rbind(out, clu)
 }
 }

 #export
 write.table(out, file= "DMRs_S_specific_intersection_clusters.txt", sep="\t", quote=F, row.names=FALSE)
 
 ########################## ATAC #############################
 ##########################################################
 
atac_bw.dir <- "/ATAC/processing/5_call_peaks"

atac_fastq.df <- read.table("/ATAC/1_diffbind/ATAC_narrow.txt", header = T, sep= "\t")[1:6, ]

atac_fastq.df$bw <- sapply(atac_fastq.df$File.Name, function(x) {
 atac_bw <- Sys.glob(paste0(atac_bw.dir, "/*", x))
 print(atac_bw)
})
atac_bw.vec <- atac_fastq.df$bw
names(atac_bw.vec) <- atac_fastq.df$Sample.Name

#use WGBS_DMR granges to obtain DMRs rownames in WGBS heatmap

atac_mat <- dmr.norm.mat.diff.locus(peak.bed.vec = atac_bw.vec, peak.bed.df = atac_fastq.df ,dmr.gr = WGBS_DMRs.gr, threads = 16) 

t_atac.mat <-  scale((t(as.data.frame(atac_mat))), center = F, scale = apply((t(as.data.frame(atac_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_atac.mat[is.na(t_atac.mat)] <- 0
dimnames(t_atac.mat) <- dimnames(atac_mat)
atac <- Heatmap(matrix = (t_atac.mat), show_row_names = F, column_title = "ATAC" , cluster_columns = FALSE , column_title_gp = gpar(fontsize = 9,fontface = "bold") ,column_names_gp = gpar(fontsize = 9) , heatmap_legend_param = list(title = "ATAC"),border=T, width = 4, height=5, col=( c( "white", "darkblue")))


############################ H3K4me3 #############################
##################################################################

H3K4Me3_bw.dir <- "/cutandrun/Seongmi/WangT_CutRun/aligned"
H3K4Me3_fastq.df <- read.table("/cutandrun/Seongmi/WangT_CutRun/aligned/H3K4me3_narrowPeaks.txt", header = T, sep= "\t")[1:6, ]
H3K4Me3_fastq.df$bw <- sapply(H3K4Me3_fastq.df$File.Name, function(x) {
 H3K4Me3_bw <- Sys.glob(paste0(H3K4Me3_bw.dir, "/*", x))
 print(H3K4Me3_bw)
})
H3K4Me3_bw.vec <- H3K4Me3_fastq.df$bw
names(H3K4Me3_bw.vec) <- H3K4Me3_fastq.df$Sample.Name

#use WGBS_DMR granges to obtain DMRs rownames in WGBS heatmap
H3K4Me3_mat <- dmr.norm.mat.diff.locus(peak.bed.vec = H3K4Me3_bw.vec, peak.bed.df = H3K4Me3_fastq.df,dmr.gr = WGBS_DMRs.gr, threads = 10)
t_H3K4Me3.mat <-  scale((t(as.data.frame(H3K4Me3_mat))), center = F, scale = apply((t(as.data.frame(H3K4Me3_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K4Me3.mat[is.na(t_H3K4Me3.mat)] <- 0
dimnames(t_H3K4Me3.mat) <- dimnames(H3K4Me3_mat)
f2 = colorRamp2(seq(min(t_H3K4Me3.mat), max(t_H3K4Me3.mat), length = 3), c( "#E0000D", "#EEEEEE","#003F0D"))
H3K4Me3 <- Heatmap(matrix = (t_H3K4Me3.mat), show_row_names = F ,column_title = "H3K4Me3" , column_title_gp = gpar(fontsize = 9,fontface = "bold") ,width = 4, height=5,column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE , heatmap_legend_param = list(title = "H3K4Me3"),border=T, col = ( c( "white", "darkorange")))
                     

############################ H3K27Ac ############################
#################################################################
        
H3K27Ac_bw.dir <- "/scratch/jharrison/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned"
H3K27Ac_fastq.df <- read.table("/scratch/jharrison/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned/H3K27ac_narrowPeaks.txt", header = T, sep= "\t")[1:6, ]
H3K27Ac_fastq.df$bw <- sapply(H3K27Ac_fastq.df$File.Name, function(x) {
 H3K27Ac_bw <- Sys.glob(paste0(H3K27Ac_bw.dir, "/*", x))
 print(H3K27Ac_bw)
})
H3K27Ac_bw.vec <- H3K27Ac_fastq.df$bw
names(H3K27Ac_bw.vec) <- H3K27Ac_fastq.df$Sample.Name
#use WGBS_DMR granges to obtain DMRs rownames in WGBS heatmap

H3K27Ac_mat <- dmr.norm.mat.diff.locus(peak.bed.vec = H3K27Ac_bw.vec, peak.bed.df = H3K27Ac_fastq.df,dmr.gr = WGBS_DMRs.gr, threads = 10)
col_fun = colorRamp2(c(1000, 0), c("green", "white"))

t_H3K27Ac.mat <-scale((t(as.data.frame(H3K27Ac_mat))), center = F, scale = apply((t(as.data.frame(H3K27Ac_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K27Ac.mat[is.na(t_H3K27Ac.mat)] <- 0
dimnames(t_H3K27Ac.mat) <- dimnames(H3K27Ac_mat)
f2 = colorRamp2(seq(min(t_H3K27Ac.mat), max(t_H3K27Ac.mat), length = 3), c("#FF7100","#EEEEEE","#00DADE"))
H3K27Ac <- Heatmap(matrix = (t_H3K27Ac.mat), show_row_names = F , column_title = "H3K27ac" , column_title_gp = gpar(fontsize = 9,fontface = "bold") , width = 4, height=5, column_names_gp = gpar(fontsize = 9), col = ( c( "white", "purple")), cluster_columns = FALSE ,border=T, heatmap_legend_param = list(title = "H3K27Ac"))


############################### H3K27Me3 ############################
######################################################################
        
H3K27Me3_bw.dir <- "/scratch/jharrison/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned"
H3K27Me3_fastq.df <- read.table("/scratch/jharrison/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned/H3K27me3_broadPeaks.txt", header = T, sep= "\t")[1:6, ]
H3K27Me3_fastq.df$bw <- sapply(H3K27Me3_fastq.df$File.Name, function(x) {
 H3K27Me3_bw <- Sys.glob(paste0(H3K27Me3_bw.dir, "/*", x))
 print(H3K27Me3_bw)
})
H3K27Me3_bw.vec <- H3K27Me3_fastq.df$bw
names(H3K27Me3_bw.vec) <- H3K27Me3_fastq.df$Sample.Name
H3K27Me3_mat <- dmr.broad.norm.mat.diff.locus(peak.bed.vec = H3K27Me3_bw.vec, peak.bed.df = H3K27Me3_fastq.df,dmr.gr = WGBS_DMRs.gr, threads = 10)


t_H3K27Me3.mat <-  scale((t(as.data.frame(H3K27Me3_mat))), center = F, scale = apply((t(as.data.frame(H3K27Me3_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K27Me3.mat[is.na(t_H3K27Me3.mat)] <- 0
dimnames(t_H3K27Me3.mat) <- dimnames(H3K27Me3_mat)
f2 = colorRamp2(seq(min(t_H3K27Me3.mat), max(t_H3K27Me3.mat), length = 3), c("#006DE8", "#EEEEEE", "#E87B00"))
H3K27Me3 <- Heatmap(matrix = (t_H3K27Me3.mat),width = 4, height=5,  show_row_names = F , column_title = "H3K27Me3" , column_title_gp = gpar(fontsize = 9,fontface = "bold") ,column_names_gp = gpar(fontsize = 9), col = ( c( "white", "darkred")), cluster_columns = FALSE , border=T, heatmap_legend_param = list(title = "H3K27Me3"))


        
############## DMR heatmaps combined ################
#####################################################

combined_noH3K27me3 <- WGBS + atac +  H3K27Ac + H3K4Me3 
combined_noH3K27me3

pdf("S_specific_combined_narrow_cluster_FRIP_narrow_with_norm.pdf", width = 6, height = 6)
combined_noH3K27me3
dev.off()

combined <- WGBS + atac +  H3K27Ac + H3K4Me3 + H3K27me3
combined

tiff("S_specific_combined_narrow_cluster_FRIP_narrow_with_H3K27me3_broad_norm.tiff", width = 6, height = 4)
combined
dev.off()

################################################################################################################################################################
#                                                                         
  #                                                                      Casp1 heatmap 
        
 ############################################################################################################################################################       
        
#uses normalized signal using read coverage and FRIP scores from bigwig files

 # begins after WGBS matrix/heatmap commands: use dmr.mat.wgbs to make WGBS_heatmap
        
  dmr.mat <- function(bw.vec, dmr.gr, threads = 5) {
  df <- utilsFanc::safelapply(names(bw.vec), function(sample) {
    #browser()
    bw <- bw.vec[[sample]]
    data.in.dmr <- rtracklayer::import(con = bw)
    
    j <- plyranges::join_overlap_left(dmr.gr, data.in.dmr)
    
    df <- j %>% as.data.frame() %>% 
      dplyr::group_by(dmr.locus) %>% dplyr::summarise(area= sum(score * (end - start))) %>% 
      dplyr::ungroup() %>% as.data.frame() %>% dplyr::select(dmr.locus, area)
    colnames(df) <- c("dmr.locus", sample)
    return(df)
  }, threads = threads) %>% Reduce(left_join, .)
  df[is.na(df)] <- 0
  #browser()
  rownames(df) <- df$dmr.locus
  df$dmr.locus <- NULL
  mat <- as.matrix(df)
  return(mat)
}


######## Cut&Run ##########
###########################
        
#######  H3K4me3 ######
#######################

#bigwigs normalized using read coverage and 1/FRIP score
        
H3K4Me3_fastq.df <- read.table("/cutandrun/Seongmi/WangT_CutRun/aligned/H3K4Me3_narrow_FRIP_inverse_bw.txt", 
                    header = T, sep= "\t")[1:6, ]
H3K4Me3_bw.dir <- "/cutandrun/Seongmi/WangT_CutRun/aligned"
H3K4Me3_fastq.df$bw <- sapply(H3K4Me3_fastq.df$File.Name, function(x) {
 H3K4Me3_bw <- Sys.glob(paste0(H3K4Me3_bw.dir, "/*", x))
})

H3K4Me3_bw.vec <- H3K4Me3_fastq.df$bw
names(H3K4Me3_bw.vec) <- H3K4Me3_fastq.df$Sample.Name
H3K4Me3_mat <- dmr.mat(bw.vec = H3K4Me3_bw.vec, dmr.gr = WGBS_DMRs.gr, threads = 10)
t_H3K4Me3.mat <-  scale((t(as.data.frame(H3K4Me3_mat))), center = F, scale = apply((t(as.data.frame(H3K4Me3_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K4Me3.mat[is.na(t_H3K4Me3.mat)] <- 0
dimnames(t_H3K4Me3.mat) <- dimnames(H3K4Me3_mat)
H3K4Me3 <- Heatmap(matrix = (t_H3K4Me3.mat), show_row_names = T ,width = 4, height=5,column_title = "H3K4Me3" , column_title_gp = gpar(fontsize = 10,fontface = "bold") , cluster_columns = FALSE ,column_names_gp = gpar(fontsize = 9), border=T, heatmap_legend_param = list(title = "H3K4Me3"), col = ( c( "white", "darkorange")))


############  H3K27Ac ############
#################################
        
H3K27Ac_fastq.df <- read.table("/cutandrun/Seongmi/WangT_CutRun/aligned/H3K27Ac_narrow_FRIP_inverse_bw.txt", 
                    header = T, sep= "\t")[1:6, ]
H3K27Ac_bw.dir <- "/cutandrun/Seongmi/WangT_CutRun/aligned"

H3K27Ac_fastq.df$bw <- sapply(H3K27Ac_fastq.df$File.Name, function(x) {
 H3K27Ac_bw <- Sys.glob(paste0(H3K27Ac_bw.dir, "/*", x))
})

H3K27Ac_bw.vec <- H3K27Ac_fastq.df$bw
names(H3K27Ac_bw.vec) <- H3K27Ac_fastq.df$Sample.Name
H3K27Ac_mat <- dmr.mat(bw.vec = H3K27Ac_bw.vec , dmr.gr = WGBS_DMRs.gr, threads = 16)
t_H3K27Ac.mat <-scale((t(as.data.frame(H3K27Ac_mat))), center = F, scale = apply((t(as.data.frame(H3K27Ac_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K27Ac.mat[is.na(t_H3K27Ac.mat)] <- 0
dimnames(t_H3K27Ac.mat) <- dimnames(H3K27Ac_mat)
H3K27Ac <- Heatmap(matrix = (t_H3K27Ac.mat), width = 4, height=5, show_row_names = T , column_title = "H3K27ac" , column_title_gp = gpar(fontsize = 10,fontface = "bold") , column_names_gp = gpar(fontsize = 9), col = c( "white", "#6C3483"), border=T, cluster_columns = FALSE , heatmap_legend_param = list(title = "H3K27Ac"))



############  H3K27Me3 ############
###################################
        
H3K27Me3_fastq.df <- read.table("/cutandrun/Seongmi/WangT_CutRun/aligned/H3K27Me3_broad_FRIP_inverse_bw.txt", 
                    header = T, sep= "\t")[1:6, ]
H3K27Me3_bw.dir <- "/cutandrun/Seongmi/WangT_CutRun/aligned"
H3K27Me3_fastq.df$bw <- sapply(H3K27Me3_fastq.df$File.Name, function(x) {
 H3K27Me3_bw <- Sys.glob(paste0(H3K27Me3_bw.dir, "/*", x))
 print(H3K27Me3_bw)
})

H3K27Me3_bw.vec <- H3K27Me3_fastq.df$bw
names(H3K27Me3_bw.vec) <- H3K27Me3_fastq.df$Sample.Name
H3K27Me3_mat <- dmr.mat(bw.vec = H3K27Me3_bw.vec, dmr.gr = WGBS_DMRs.gr, threads = 16)
t_H3K27Me3.mat <-  scale((t(as.data.frame(H3K27Me3_mat))), center = F, scale = apply((t(as.data.frame(H3K27Me3_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_H3K27Me3.mat[is.na(t_H3K27Me3.mat)] <- 0
dimnames(t_H3K27Me3.mat) <- dimnames(H3K27Me3_mat)
H3K27Me3 <- Heatmap(matrix = (t_H3K27Me3.mat),width = 4, height=5, show_row_names = T , column_title = "H3K27Me3" ,column_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 10,fontface = "bold") , col = c( "white", "darkred"), border=T, cluster_columns = FALSE , heatmap_legend_param = list(title = "H3K27Me3"))

############# ATAC ############
###############################
        
atac_fastq.df <- read.table("/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned/ATAC_narrow_FRIP_inverse_bw.txt", header = T, sep= "\t")[1:6, ]
atac_bw.dir <- "/Seongmi/cutandrun/Seongmi/WangT_CutRun/aligned"

atac_fastq.df$bw <- sapply(atac_fastq.df$File.Name, function(x) {
 atac_bw <- Sys.glob(paste0(atac_bw.dir, "/*", x))
 print(atac_bw)
})

atac_bw.vec <- atac_fastq.df$bw
names(atac_bw.vec) <- atac_fastq.df$Sample.Name

atac_mat <-dmr.mat(bw.vec = atac_bw.vec, dmr.gr = WGBS_DMRs.gr, threads = 16)


t_atac.mat <-  scale((t(as.data.frame(atac_mat))), center = F, scale = apply((t(as.data.frame(atac_mat))), 2, sd, na.rm = TRUE)) %>% as.data.frame %>% t %>% as.matrix 
t_atac.mat[is.na(t_atac.mat)] <- 0
dimnames(t_atac.mat) <- dimnames(atac_mat)
atac <- Heatmap(matrix = (t_atac.mat), show_row_names = T,width = 4, height=5, column_title = "ATAC" ,column_names_gp = gpar(fontsize = 9), cluster_columns = FALSE , column_title_gp = gpar(fontsize = 10,fontface = "bold") , border=T, heatmap_legend_param = list(title = "ATAC"), col=( c( "white", "#3498DB")))
 
### combine graphs ####
combined <- WGBS + atac +  H3K27Ac + H3K4Me3 
combined
        
## casp1 is row 188 

tiff("/scratch/jharrison/Seongmi/WGBS/May_DMR/analysis/casp1_DMR_S_specific_intersection_bigwig_combined_root-mean-square_FRIP_norm.pdf", width = 6, height = 4)
combined[188,]
dev.off()



   
        
        
        



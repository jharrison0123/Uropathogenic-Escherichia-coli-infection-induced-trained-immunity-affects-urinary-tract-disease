#R 3.6.1

library(ComplexHeatmap)
library(dplyr)
library(devtools)
library(ggplot2)
library(GenomicRanges)
library(tidyverse)


###### Promoter ######

### Read in differential expression data from differentiated urothelia: Sensitized vs. Resolved ###
########################

DUC <- read.table(file="/RNA_infected/NRS_DUC_deseq2_norm[17].txt", header=T, sep="\t",fill = TRUE)

# WGBS_mat is WGBS DMR matrix with methylation values for Sensitized and Resolved USCs 
## comes from WGBS DMR heamap commands (see DMR heatmap R script)

# Promoter and transcript annotation of Senstizied specific hyopDMRs- annotated using commands in DMR_figures.sh in WGBS dir under section :
## DMR annotation of Sensitized-specific hypo-DMRs extracted from complexheatmap cluster1- hypo DMRs in sensitized samples

promoter <- read.table(file="/WGBS/S_specific/DMRs_S_specific_cluster1_annotation.bed.PeakGenicIntersection_upstream_1kb.txt", quote="", header=F,fill = TRUE) %>% as.data.frame 
transcripts <- read.table(file="/WGBS/S_specific/DMRs_S_specific_cluster1_annotation.bed.PeakGenic_intersection_promoter_genes.txt") %>% as.data.frame

colnames(transcripts) <- c("dmr.locus","gene","ensemblID")
DUC_DMR_promoter <- DUC[DUC$Gene.id %in% transcripts$ensemblID,] 
DUC_DMR_promoter <- arrange(DUC_DMR_promoter, Gene.id)

rownames(DUC_DMR_promoter) <- DUC_DMR_promoter$Gene 
DUC_DMR_dmr <- transcripts[ transcripts$ensemblID %in% DUC$Gene.id,]
DUC_DMR_dmr <- arrange(DUC_DMR_dmr, dmr.locus)

DUC_DMR_WGBS_mat <- WGBS_mat[(rownames(WGBS_mat) %in% DUC_DMR_dmr$dmr.locus),]  %>% as.data.frame()
DUC_DMR_WGBS_mat$dmr.locus <- rownames(DUC_DMR_WGBS_mat)
DUC_DMR_dmr <- arrange(DUC_DMR_dmr, dmr.locus)
DUC_DMR_WGBS_mat <- arrange(DUC_DMR_WGBS_mat, dmr.locus)
DUC_DMR_WGBS_mat$gene <- DUC_DMR_dmr$gene

DUC_DMR_WGBS_mat_V1 <- aggregate(DUC_DMR_WGBS_mat, by = list(DUC_DMR_WGBS_mat$gene), FUN= mean)
DUC_DMR_WGBS_mat_V2 <- DUC_DMR_WGBS_mat_V1[,c(1:7)]

names(DUC_DMR_WGBS_mat_V2)[names(DUC_DMR_WGBS_mat_V2) == 'Group.1'] <- 'Gene'

  

#DMR vs RNA log2 fold change between sensitized and resolved

###################### RNA R vs S fold change ###################### 
#################################################################### 
#R+ PBS 9:12
#R+ infection 13:18
#S + PBS 19:24
#S + infection 25:30

R <- DUC_DMR_promoter[,c(11:14)] 
df <- rowMeans(DUC_DMR_promoter[,c(11:14)] , na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
R_list <- df$sum

R_infec <- DUC_DMR_promoter[,c(15:20)]
df <- rowMeans(DUC_DMR_promoter[,c(13:18)], na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
R_infec_list  <- df$sum

S <- DUC_DMR_promoter[,c(21:26)]
df  <- rowMeans(DUC_DMR_promoter[,c(19:24)], na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
S_list  <- df$sum

S_infec <- DUC_DMR_promoter[,c(27:32)]
df <- rowMeans(DUC_DMR_promoter[,c(25:30)], na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
S_infec_list <- df$sum



f <- foldchange(S_list, R_list)
lf <- foldchange2logratio(f, base = 2)

R_S_foldchange <- cbind(R_list, S_list, lf,1) 
row.names(R_S_foldchange) <- DUC_DMR_promoter$Gene

f <- foldchange( S_infec_list,R_infec_list) 
lf <- foldchange2logratio(f, base = 2)

R_S_infec_foldchange <- cbind(R_infec_list, S_infec_list, lf,2)
row.names(R_S_infec_foldchange) <- DUC_DMR_promoter$Gene

R_S_foldchange[is.na(R_S_foldchange)] <- 1



R_S_infec_foldchange[is.na(R_S_infec_foldchange)] <- 1

###################### WBGS R vs S fold change ###################### 
#################################################################### 
#DUC_RNA_DMR_WGBS_mat_final

df<- rowMeans(DUC_DMR_WGBS_mat_V2[,c(4:5)], na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
R_WGBS <- df$sum

df <- rowMeans(DUC_DMR_WGBS_mat_V2[,c(6:7)], na.rm=F) %>% as.data.frame() %>% mutate(sum=.+1)
S_WGBS <- df$sum

f <- data.frame("R"= R_WGBS, "S"= S_WGBS) %>%  mutate(diff = S-R)
#lf <- foldchange2logratio(f, base = 2)

R_S_WGBS_foldchange <- cbind(R_WGBS, S_WGBS, f)
row.names(R_S_WGBS_foldchange) <- DUC_DMR_WGBS_mat_V2$gene

log_fc_df <- cbind("WGBS"=R_S_WGBS_foldchange$diff,
                        "DUC_pbs"=R_S_foldchange[,c(3)],
                        "DUC_infec"=R_S_infec_foldchange[,c(3)]) %>% as.data.frame()


log_fc_df$Gene <- rownames(log_fc_df)

log_fc_df <- log_fc_df[c(1:9,12,14),]

 write.table(log_fc_df, file= "/WGBS/analysis/DMRs_S_specific_intersection_DMR_DUC_pbs_infection_log2_fc_DMR_promoter_df.txt", sep="\t", quote=F, row.names=FALSE)


log_fc_df_DMR_promoter_melt <- melt(log_fc_df,id.vars =c("Gene", "WGBS"), measure.vars = c("DUC_pbs","DUC_infec"))


pdf("/WGBS/analysis/DMR_DUC_pbs_infection_log2_fc_DMR_promoter_with_labels.pdf")


g <- ggplot(log_fc_df_DMR_promoter_melt , aes(x=WGBS, y=value)) + 
  geom_line(aes(group = WGBS,color="#E3E3E3"), size=3) +
    geom_point(aes(color=variable), size=3.8) +
  xlim(-80,80) + 
    theme_classic() +
   #xlim(-5,5) + ylim(-5,5) +
  ggtitle("Sen vs Res ")  + # for the main title
xlab("DNA methylation change") +  # for the x axis label
ylab("Log2 (RNA Fold Change)") +
    scale_color_manual(values=c("#D3C5EC","#E3E3E3", "#5C2DAF")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = 0, color="darkgrey") +
  geom_hline(yintercept = 0, color="darkgrey") +
  geom_text(data=log_fc_df,
            aes(WGBS,DUC_infec,label=Gene, hjust=-0.25, vjust=0.25),  size=2)
  
g


dev.off()




pdf("/WGBS/analysis/DMR_DUC_pbs_infection_log2_fc_DMR_promoter_no_labels.pdf", width = 8, height = 6)


g <- ggplot(log_fc_df_DMR_promoter_melt , aes(x=WGBS, y=value)) + 
  geom_line(aes(group = WGBS,color="#E3E3E3"), size=4) +
    geom_point(aes(color=variable), size=4.8) +
  xlim(-80,80) + 
    theme_classic() +
   #xlim(-5,5) + ylim(-5,5) +
  ggtitle("Sen vs Res ")  + # for the main title
xlab("DNA methylation change") +  # for the x axis label
ylab("Log2 (RNA Fold Change)") +
    scale_color_manual(values=c("#D3C5EC","#E3E3E3", "#5C2DAF")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = 0, color="darkgrey") +
  geom_hline(yintercept = 0, color="darkgrey") #+
  #geom_text(data=log_fc_df,
           # aes(WGBS,DUC_infec,label=Gene, hjust=-0.25, vjust=0.25),  size=2)
  

g


dev.off()

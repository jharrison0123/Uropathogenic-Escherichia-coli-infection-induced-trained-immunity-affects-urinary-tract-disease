#R version 3.6.1 

library(ggplot2)
#ggplot2_3.3.6

##### WGBS ###### 
d = read.delim("WGBS_DMR_DMRs_S_specific_cluster1_annotation_clustering_profile.tsv", header=F)
d2 = as.data.frame(t(as.matrix(d[,-c(1,2)]))) 
colnames(d2) = d[,1]
d2 <- d2[1:1500,]

d3 <- mutate_all(d2[,2:8], function(x) as.numeric(as.character(x)))
d3_melt <- melt(d3, id.vars="bins")

tiff("WGBS_dmr_profile_smooth_R.tiff", width = 12, height = 6, units = 'in', res = 300, compression = 'lzw')

ggplot(data=d3_melt, aes(x=bins, y=value, color=variable)) + 
  #geom_line() +
  geom_smooth(method = 'loess', formula = 'y ~ x', span = 0.17, se=F,size=1.20) +
    theme_classic() +
    theme(
plot.title = element_text(size=25, face="bold"),
axis.title.x = element_text(size=20, face="bold"),
axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
    ggtitle("DNA Methylation")  + # for the main title
#ylab("DNA Methylation") +
  theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("#2832f7", "#0812cf","#facd07","#f5da64", "#cc6cf5","#ae02f7")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = c(500,1000), alpha=0.4, color="darkgrey", size=10)  +
      scale_y_continuous("",breaks=c(0, 25, 50, 75, 100),
              labels=c("0%", "25%", "50% ","75%", "100%"), limits=c(0, 100), expand = c(0.02, 0.02)) +
    scale_x_continuous("DMRs",breaks=c(0, 500, 1000, 1500), limits=c(0, 1500), expand = c(0.02, 0.02),
              labels=c("-5kb", "Start", "End","+5kb")) 
 #+ # for the main title

     
dev.off()

 ### H3K27Ac ####                
                 
d = read.delim("H3K27Ac_DMR_DMRs_S_specific_cluster1_annotation_clustering_diff_norm_FRIP_inverse_profile.tsv", header=F)
d2 = as.data.frame(t(as.matrix(d[,-c(1,2)]))) 
colnames(d2) = d[,1]
d2 <- d2[1:1500,]

d3 <- mutate_all(d2[,2:8], function(x) as.numeric(as.character(x)))
d3_melt <- melt(d3, id.vars="bins")


tiff("H3K27Ac_dmr_profile_smooth_R.tiff", width = 12, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(data=d3_melt, aes(x=bins, y=value, color=variable)) + 
    geom_smooth(method = 'loess', formula = 'y ~ x', span = 0.05, se=F, size=1.20) +
       ylim(10,50) + 
    theme_classic() +
    theme(
plot.title = element_text(size=25, face="bold"),
axis.title.x = element_text(size=20, face="bold"),
axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
    ggtitle("H3K27ac")  + # for the main title
ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("#2832f7", "#0812cf","#facd07","#f5da64", "#cc6cf5","#ae02f7")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = c(500,1000), alpha=0.4, color="darkgrey", size=10) +
    scale_x_continuous("DMRs",breaks=c(0, 500, 1000, 1500), limits=c(0, 1500), expand = c(0.02, 0.02),
              labels=c("-5kb", "Start", "End","+5kb")) 
 + # for the main title

     
dev.off()

#### H3K4Me3  #####                
         
d = read.delim("H3K4Me3_DMR_DMRs_S_specific_cluster1_annotation_clustering_diff_norm_FRIP_inverse_profile.tsv", header=F)
d2 = as.data.frame(t(as.matrix(d[,-c(1,2)]))) 
colnames(d2) = d[,1]
d2 <- d2[1:1500,]

d3 <- mutate_all(d2[,2:8], function(x) as.numeric(as.character(x)))
d3_melt <- melt(d3, id.vars="bins")


tiff("H3K4Me3_dmr_profile_smooth_R.tiff", width = 12, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(data=d3_melt, aes(x=bins, y=value, color=variable)) + 
  #geom_line() +
    geom_smooth(method = 'loess', formula = 'y ~ x', span = 0.05, se=F, size=1.20) +
       ylim(10,50) + 
    theme_classic() +
    theme(
plot.title = element_text(size=25, face="bold"),
axis.title.x = element_text(size=20, face="bold"),
axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
    ggtitle("H3K4me3")  + # for the main title
ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("#2832f7", "#0812cf","#facd07","#f5da64", "#cc6cf5","#ae02f7")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = c(500,1000), alpha=0.4, color="darkgrey", size=10) +
    scale_x_continuous("DMRs",breaks=c(0, 500, 1000, 1500), limits=c(0, 1500), expand = c(0.02, 0.02),
              labels=c("-5kb", "Start", "End","+5kb")) 

 + # for the main title

     
dev.off()

##### H3K27Me3 #####                 
                 
d = read.delim("H3K27Me3_DMR_DMRs_S_specific_cluster1_annotation_clustering_diff_norm_FRIP_inverse_profile.tsv", header=F)
d2 = as.data.frame(t(as.matrix(d[,-c(1,2)]))) 
colnames(d2) = d[,1]
d2 <- d2[1:1500,]

d3 <- mutate_all(d2[,2:8], function(x) as.numeric(as.character(x)))
d3_melt <- melt(d3, id.vars="bins")


tiff("H3K27Me3_dmr_profile_smooth_R.tiff", width = 12, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(data=d3_melt, aes(x=bins, y=value, color=variable)) + 
  geom_smooth(method = 'loess', formula = 'y ~ x', span = 0.01, se=F,size=1.20) +
  #geom_line() +
  #geom_smooth(method = NULL, se = FALSE) +
     ylim(0,5) +
    theme_classic() +
    theme(
plot.title = element_text(size=25, face="bold"),
axis.title.x = element_text(size=20, face="bold"),
axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
    ggtitle("H3K27me3")  + # for the main title
ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("#2832f7", "#0812cf","#facd07","#f5da64", "#cc6cf5","#ae02f7")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = c(500,1000), alpha=0.4, color="darkgrey", size=10) +
    scale_x_continuous("DMRs",breaks=c(0, 500, 1000, 1500), limits=c(0, 1500), expand = c(0.02, 0.02),
              labels=c("-5kb", "Start", "End","+5kb")) 

 + # for the main title

     
dev.off()

##### ATAC ##### 
                 
d = read.delim("atac_DMR_DMRs_S_specific_cluster1_annotation_clustering_diff_norm_FRIP_inverse_profile.tsv", header=F)
d2 = as.data.frame(t(as.matrix(d[,-c(1,2)]))) 
colnames(d2) = d[,1]
d2 <- d2[1:1500,]

d3 <- mutate_all(d2[,2:8], function(x) as.numeric(as.character(x)))
d3_melt <- melt(d3, id.vars="bins")


tiff("ATAC_dmr_profile_smooth_R.tiff", width = 12, height = 6, units = 'in', res = 300, compression = 'lzw')
ggplot(data=d3_melt, aes(x=bins, y=value, color=variable)) + 
  geom_smooth(method = 'loess', formula = 'y ~ x', span = 0.01, se=F,size=1.20) +
       ylim(10,50) + 
   # ylim(0,100) + 
    theme_classic() +
    theme(
plot.title = element_text(size=25, face="bold"),
axis.title.x = element_text(size=20, face="bold"),
axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold")) +
      ggtitle("ATAC")  + # for the main title
ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values=c("#2832f7", "#0812cf","#facd07","#f5da64", "#cc6cf5","#ae02f7")) +
  theme(panel.background=element_rect(colour="black")) +
  geom_vline(xintercept = c(500,1000), alpha=0.4, color="darkgrey", size=10) +
    scale_x_continuous("DMRs",breaks=c(0, 500, 1000, 1500), limits=c(0, 1500), expand = c(0.02, 0.02),
              labels=c("-5kb", "Start", "End","+5kb")) 
  #geom_smooth(method = NULL, se = FALSE) +


 + # for the main title

     
dev.off()


                 
                 


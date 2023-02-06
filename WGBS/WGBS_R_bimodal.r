#R version 3.6.1 

library(ggplot2)

##using methylc track
mn1 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/N1._methylc_5x.txt", header=F))
mn3 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/N3._methylc_5x.txt", header=F))
mr1 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/R1._methylc_5x.txt", header=F))
mr4 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/R4._methylc_5x.txt", header=F))
ms1 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/S1._methylc_5x.txt", header=F))
ms2 <- data.frame(read.table(file="/WGBS/6_track_mergedCG/rename/S2._methylc_5x.txt", header=F))


bimodal_methylc <-ggplot(data=mn1, aes(V5)) + 
    geom_density( colour='#2832f7', fill='#2832f7',   alpha=0.05,  adjust=3) +
  geom_density(data = mn3, alpha=0.05, colour='#0812cf', fill='#0812cf',  alpha=0.05,  adjust=3) + 
  geom_density(data = mr1, alpha=0.05, colour='#facd07', fill='#facd07',  alpha=0.05,  adjust=3) +
    geom_density(data = mr4, alpha=0.05, colour='#f5da64', fill='#f5da64',  alpha=0.05,  adjust=3) +
    geom_density(data = ms1, alpha=0.05, colour='#cc6cf5', fill='#cc6cf5',  alpha=0.05,  adjust=3) +
    geom_density(data = ms2, alpha=0.05, colour='#ae02f7', fill='#ae02f7',  alpha=0.05,  adjust=3) + 
  xlab("CpG Methylation") +
  ylab("Density") + 
  scale_fill_manual(name='', breaks=c('N1', 'N3', 'R1', 'R4','S1' , 'S2'), values=c('N1' = '#2832f7', 'N3' = '#0812cf', 'R1' = '#facd07', 'R4'= '#f5da64','S1' = '#cc6cf5', 'S2' = '#ae02f7')) +
   scale_color_manual(name='', breaks=c('N1', 'N3', 'R1', 'R4','S1' , 'S2'), values=c('N1' = '#2832f7', 'N3' = '#0812cf', 'R1' = '#facd07', 'R4'= '#f5da64','S1' = '#cc6cf5', 'S2' = '#ae02f7')) +
  theme_bw() +
  scale_x_continuous("CpG Methylation",breaks=c(0, .25, .50, .75 , 1.00), limits=c(0, 1), expand = c(0.01, 0.01)) +     
  theme(text = element_text(size=30)) 

   
    #guides(fill = guide_legend(override.aes = list(color="NA", size=10, stroke=0.01) ) ) +

bimodal_methylc

pdf("bimodal_methylc_adjust3.pdf")
bimodal_methylc
dev.off()

naive1_cov.df <- as.data.frame(read.table("/scratch/jharrison/Seongmi/WGBS/2_bismark/cov_rename/N1_deduplicated.bismark.cov_coverage_5x.txt", header=F))
#R version 3.6.1 

library(ggplot2)

naive3_cov.df <- as.data.frame(read.table("N3_deduplicated.bismark.cov_coverage_5x.txt", header=F))
resolved1_cov.df <- as.data.frame(read.table("R1_deduplicated.bismark.cov_coverage_5x.txt", header=F))
resolved4_cov.df <- as.data.frame(read.table("R4_deduplicated.bismark.cov_coverage_5x.txt", header=F))
sensitized1_cov.df <- as.data.frame(read.table("S1_deduplicated.bismark.cov_coverage_5x.txt", header=F))
sensitized2_cov.df <- as.data.frame(read.table("S2_deduplicated.bismark.cov_coverage_5x.txt", header=F))



naive1_cov.plot <- ggplot(naive1_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="darkblue", fill="darkblue") 
naive3_cov.plot <- ggplot(naive3_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="darkgreen", fill="darkgreen") 
resolved1_cov.plot <- ggplot(resolved1_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="darkred", fill="darkred") 
resolved4_cov.plot <- ggplot(resolved4_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="yellow", fill="yellow") 
sensitized1_cov.plot <- ggplot(sensitized1_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="darkorange", fill="darkorange") 
sensitized2_cov.plot <- ggplot(sensitized2_cov.df, aes(V3)) + geom_histogram( binwidth = 7, alpha=0.2, colour="purple", fill="purple") 

colors <- c('N1'="darkblue", 'N3'="darkgreen", 'R1'="darkred", 'R4'="darkyellow", 'S1'="darkorange", 'S2'="darkpurple")

bimodal <-ggplot(data=naive1_cov.df, aes(V3)) + 
  geom_histogram( alpha=0.2, colour="darkblue", fill="darkblue", binwidth = 7) + 
  geom_histogram(data = naive3_cov.df, alpha=0.2, colour="darkgreen", fill="darkgreen", binwidth = 7) + 
  geom_histogram(data = resolved1_cov.df, alpha=0.2, colour="darkred", fill="darkred", binwidth = 7) +
    geom_histogram(data = resolved4_cov.df, alpha=0.2, colour="yellow", fill="yellow", binwidth = 7) +
    geom_histogram(data = sensitized1_cov.df, alpha=0.2, colour="darkorange", fill="darkorange", binwidth = 7) +
    geom_histogram(data = sensitized2_cov.df, alpha=0.2, colour="purple", fill="purple", binwidth = 7) + labs(x="CpG Methylation (%)", y="Density", colors = "Legend") + legend('topright', c('N1', 'N3', 'R1', 'R4', 'S1', 'S2'), fill=c('blue', 'green', 'red', 'yellow', 'orange', 'purple'))
   
bimodal

pdf("bimodal.pdf")
bimodal
dev.off()

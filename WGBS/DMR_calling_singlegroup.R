#R/4.1.2
#DSS v2.43.2

library(DSS)
#require(bsseq)
argv=commandArgs(trailingOnly=T)

file1 <- argv[1]
file2 <- argv[2]
dir <- argv[3]

names <- c(file1, file2)
names <- gsub("_DSS.txt", "", names)
names <- gsub("^.*/", "", names)
base  <- paste0(names[1], "_", names[2])

# Load the data
dat1 <- read.table(file1, header=TRUE)
dat2 <- read.table(file2, header=TRUE)

BSobj <- makeBSseqData( list(dat1, dat2), names )
#BSobj
print("BSobj made")

# Perform statistical test for DML
dmlTest <- DMLtest(BSobj, group1=names[1], group2=names[2], smoothing=TRUE, smoothing.span=500, ncores=8)
head(dmlTest)
wfile = paste0(dir,"/", base, ".txt.gz")
write.table(dmlTest, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

##Calling DMRs with 0.05
delta <- 0.2

# DMR detection
# minLenght>=200, minCpG>=5, p<0.05, delta=0.2

dmrs <- callDMR(dmlTest, delta=delta, p.threshold=0.05, minlen=200, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p05.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)


dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=0.001, minlen=200, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p001.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=10^-5, minlen=200, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p10_5.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=10^-16, minlen=200, minCG=5, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p10_16.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)


sessionInfo()

q()

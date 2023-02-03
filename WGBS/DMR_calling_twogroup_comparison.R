#R-4.1.2

library(DSS)
require(bsseq)
argv=commandArgs(trailingOnly=T)

file1 <- argv[1]
file2 <- argv[2]
file3 <- argv[3]
file4 <- argv[4]
dir <- argv[5]

names <- c(file1, file2, file3, file4)
names <- gsub("_R1_bismark_bt2_pe.deduplicated.bismark.cov.gz_DSS.txt", "", names)
names <- gsub("^.*/", "", names)
base  <- paste0(names[1], "_", names[2], "vs", names[3], "_", names[4])

# Load the data
dat1 <- read.table(file1, header=TRUE)
dat2 <- read.table(file2, header=TRUE)
dat3 <- read.table(file3, header=TRUE)
dat4 <- read.table(file4, header=TRUE)

BSobj <- makeBSseqData( list(dat1, dat2, dat3, dat4), names )
#BSobj
print("BSobj made")

# Perform statistical test for DML
dmlTest <- DMLtest(BSobj, group1=c(names[1], names[2]), group2=c(names[3], names[4]), smoothing=TRUE, smoothing.span=500, ncores=8)
head(dmlTest)
wfile = paste0(dir,"/", base, "twogroup.txt.gz")
write.table(dmlTest, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

##Calling DMRs with 0.05

# DMR detection
# stringent threshold: minLenght=200, minCpG=5, delta=0.2, dis.merge=50, pct.sig=0.5, p.threshold = 10^-16
#parameters based on https://www.nature.com/articles/s41586-020-2135-x

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=0.05, minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p05_twogroup.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=0.001, minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p001_twogroup.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=10^-5, minlen=50, minCG=3, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p10_5_twogroup.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)

dmrs <- callDMR(dmlTest, delta=0.2, p.threshold=10^-16, minlen=50, minCG=5, dis.merge=50, pct.sig=0.5)
dmrdir <-paste0(dir)
wfile = paste0(dmrdir, "/DMRs_", base, "_p10_16_twogroup.txt.gz")
write.table(dmrs, file = gzfile(wfile, compression = 9), quote = FALSE, sep = "\t", row.names = FALSE)


sessionInfo()

q()

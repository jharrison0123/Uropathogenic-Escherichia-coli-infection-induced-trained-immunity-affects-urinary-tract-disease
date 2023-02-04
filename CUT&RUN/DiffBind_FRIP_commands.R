#R 3.6.1
## Step 1. DiffBind
## https://bioconductor.org/packages/3.10/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/ChIPSeq/scripts/ChIP_Practical3_DiffBind.html

setwd("1_diffbind")

## 1.1. read in data
samples = read.csv("sample.csv")
samples_dba = dba(sampleSheet=samples)

## 1.2. get count matrix
samples_count = dba.count(samples_dba) 
#output normalized counts
cutandrun_counts <- dba.peakset(samples_count, bRetrieve=T, DataType=DBA_DATA_FRAME)
write.table(counts, file='/1_diffbind/normalized_count.txt', sep='\t', row.names = FALSE, col.names=TRUE, quote=FALSE)



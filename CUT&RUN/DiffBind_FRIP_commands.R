#R 3.6.1
## Step 1. DiffBind
## https://bioconductor.org/packages/3.10/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
## https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/ChIPSeq/scripts/ChIP_Practical3_DiffBind.html

setwd("1_diffbind")

##CUT&RUN 
## 1.1. read in data
samples = read.csv("CUT&RUN_sample_table.csv")
#usinng deduplicated bam and narrowpeak (H3K27Ac and H3K3me3) and broad (H3K27me3)
samples_dba = dba(sampleSheet=samples)

## 1.2. get count matrix
samples_count = dba.count(samples_dba) 

###ATAC
samples = read.csv("ATAC_sample_table.csv")
#usinng deduplicated bam and narrowpeak (H3K27Ac and H3K3me3) and broad (H3K27me3)
samples_dba = dba(sampleSheet=samples)

## 1.2. get count matrix
samples_count = dba.count(samples_dba) 


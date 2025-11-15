# load packages and data

install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"))

library(TCGAbiolinks)

query <- GDCquery(
  project = c("TCGA-LGG", "TCGA-GBM"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga_data <- GDCprepare(query)



#Loading and Processing Normal Tissue Data (GTEx)

# Read the massive GTEx TPM file (this takes a moment and requires sufficient RAM)
install.packages("data.table")
library(data.table)

gtex_tpm <- fread("C:/Users/jerry/OneDrive/Documents/bioinformatics/Final project/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct", 
                  sep = "\t", 
                  skip = 2)


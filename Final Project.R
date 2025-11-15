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

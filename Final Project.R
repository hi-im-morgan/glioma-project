install.packages("BiocManager")

BiocManager::install(c("recount3", "SummarizedExperiment", "recount"))

library(recount)
library(SummarizedExperiment)

library(recount3)

human_projects <- available_projects("human")

tcga_projects <- subset(human_projects, file_source == "tcga")
head(tcga_projects)

lgg_info <- subset(tcga_projects, project == "LGG")

gbm_info <- subset(tcga_projects, project == "GBM")

lgg_info
gbm_info

rse_lgg <- create_rse(lgg_info)
rse_gbm <- create_rse(gbm_info)

colnames(colData(rse_lgg))
colnames(colData(rse_gbm))


#GTEX Transcripts per Million

# Read the massive GTEx TPM file (this takes a moment and requires sufficient RAM)
install.packages("data.table")
library(data.table)

gtex_tpm <- fread("C:/Users/jerry/OneDrive/Documents/bioinformatics/Final project/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_tpm_non_lcm.gct", 
                  sep = "\t", 
                  skip = 2)



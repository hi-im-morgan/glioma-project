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

rse_lgg <- create_rse(lgg_info)
rse_gbm <- create_rse(gbm_info)

# Check genes line up
stopifnot(identical(rownames(rse_lgg), rownames(rse_gbm)))

rse_glioma <- cbind(rse_lgg, rse_gbm)

glioma_counts <- assay(rse_glioma)
glioma_meta <- as.data.frame(colData(rse_glioma))

library(data.table)

gtex_projects <- available_projects("human")[
  available_projects("human")$project == "BRAIN",
]

gtex_rse <- create_rse(gtex_projects)
gtex_counts <- assay(gtex_rse)
gtex_meta <- as.data.frame(colData(gtex_rse))

# Install stuff for RNA seq analysis
BiocManager::install(
  c("edgeR", "limma", "DESeq2", "biomaRt", "pheatmap", "ggplot2"),
  ask = FALSE,
  update = FALSE
)
library(edgeR)
library(limma)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(ggplot2)

dim(glioma_counts)
dim(gtex_counts)
stopifnot(ncol(glioma_counts) > 0, ncol(gtex_counts) > 0)

clean_ids <- function(genes) sub("\\..*$", "", genes) # ENSG00000.1 -> ENSG00000
glioma_ids <- clean_ids(rownames(glioma_counts))
gtex_ids <- clean_ids(rownames(gtex_counts))

# keep gene names in memory for later mapping
rownames(glioma_counts) <- glioma_ids
rownames(gtex_counts) <- gtex_ids

##### 3. Intersect genes and subset matrices #####
common_genes <- intersect(rownames(glioma_counts), rownames(gtex_counts))
length(common_genes) # how many genes are shared
if (length(common_genes) < 10000) {
  warning("Low overlap of genes — check gene IDs/annotations")
}

glioma_mat <- glioma_counts[common_genes, , drop = FALSE]
gtex_mat <- gtex_counts[common_genes, , drop = FALSE]

##### 4. Build combined matrix + metadata #####
combined <- cbind(glioma_mat, gtex_mat)
mean(combined == 0)


# create sample-level table
meta_glioma <- glioma_meta
meta_gtex <- gtex_meta

# Make sure rownames of metadata match colnames of matrices
rownames(meta_glioma) <- colnames(glioma_mat)
rownames(meta_gtex) <- colnames(gtex_mat)

# All columns that appear in either dataset
all_cols <- union(colnames(glioma_meta), colnames(gtex_meta))

# Add any missing columns with NA
glioma_meta2 <- glioma_meta
gtex_meta2 <- gtex_meta

for (col in all_cols) {
  if (!col %in% colnames(glioma_meta2)) {
    glioma_meta2[[col]] <- NA
  }
  if (!col %in% colnames(gtex_meta2)) {
    gtex_meta2[[col]] <- NA
  }
}

# Reorder columns identically
glioma_meta2 <- glioma_meta2[, all_cols]
gtex_meta2 <- gtex_meta2[, all_cols]

# Make sure rownames match sample IDs
rownames(glioma_meta2) <- colnames(glioma_mat)
rownames(gtex_meta2) <- colnames(gtex_mat)

meta_comb <- rbind(glioma_meta2, gtex_meta2)

# Create source / group variables
meta_comb$source <- ifelse(
  rownames(meta_comb) %in% colnames(glioma_mat),
  "TCGA",
  "GTEx"
)
meta_comb$group <- ifelse(meta_comb$source == "TCGA", "Glioma", "Normal")

# reorder metadata to match combined columns
meta_comb <- meta_comb[colnames(combined), , drop = FALSE]

##### 5. Basic QC #####
libsizes <- colSums(combined)
summary(libsizes)
plot(density(log10(libsizes)), main = "Library sizes (log10)")

dge_qc <- DGEList(counts = combined)
cpm_qc <- cpm(dge_qc, log = TRUE, prior.count = 1)
###############################
# SAFE PCA FOR RNA-SEQ DATA
###############################
install.packages("irlba")

library(irlba) # fast PCA

# cpm_qc = filtered CPM/TPM matrix (genes x samples)

# 1. Compute variance per gene (fast)
gene_var <- apply(cpm_qc, 1, var)

# 2. Keep top 2000 most variable genes
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
cpm_top <- cpm_qc[top_genes, ] # reduced matrix

# 3. Run fast truncated PCA (safe, won't freeze)
pca <- prcomp_irlba(
  t(cpm_top), # samples x genes
  n = 5, # top 5 PCs (enough for visualization)
  center = TRUE,
  scale. = FALSE
)

# 4. Build a data frame for plotting
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = meta_comb$group, # Normal vs Glioma
  source = meta_comb$source # GTEx vs TCGA
)

# 5. Plot with ggplot2
library(ggplot2)

ggplot(pca_df, aes(PC1, PC2, color = group, shape = source)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA of GTEx normal + TCGA glioma (top 2000 genes)")

##### 6. Filter lowly-expressed genes #####
# Keep genes with CPM > 1 in at least N samples (choose N depending on sample size)
min_samples <- 10
dge <- DGEList(counts = combined)
keep <- rowSums(cpm(dge) > 1) >= min_samples
table(keep)
dge <- dge[keep, , keep.lib.sizes = FALSE]

##### 7. Normalization & design #####
dge <- calcNormFactors(dge) # TMM normalization

# Create design (Normal = reference)
meta_comb$group <- factor(meta_comb$group, levels = c("Normal", "Glioma"))
meta_comb$source <- factor(meta_comb$source) # harmless to keep as factor

# Simple model: Glioma vs Normal only
design <- model.matrix(~group, data = meta_comb)
colnames(design)


##### 8. voom + limma #####
v <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(v, design)

# Contrast: test groupGlioma vs groupNormal while adjusting for source
# In model.matrix with "~ source + group", the coefficient for groupGlioma is the effect
fit <- eBayes(fit)

# Extract results for the group effect (name may be "groupGlioma" or "groupGlioma" depending on factor levels)
coef_name <- grep("groupGlioma", colnames(fit$coefficients), value = TRUE)
coef_name
top_all <- topTable(fit, coef = coef_name[1], number = Inf, sort.by = "P")


##### 9. Select upregulated in Glioma #####
# criteria: logFC > 1 and adj.P.Val < 0.05 (adjust as needed)
top_sig <- top_all[top_all$logFC > 1 & top_all$adj.P.Val < 0.05, ]

nrow(top_sig)
head(top_sig)

top_by_expression <- top_sig[order(-top_sig$AveExpr), ]
head(top_by_expression)

# Get the gene IDs from the rownames
candidate_genes <- rownames(top_by_expression)


library(biomaRt)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(mart)
head(attributes, 50)

proteins <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "uniprotswissprot"),
  filters = "ensembl_gene_id",
  values = candidate_genes,
  mart = mart
)
head(proteins)

library(dplyr)

# Suppose your data frame is called 'proteins'
proteins_filtered <- proteins %>%
  filter(uniprotswissprot != "") %>% # Keep only rows with non-empty UniProt IDs
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keep only the first occurrence per gene

# Check result
head(proteins_filtered)
library(httr)

top_all_df <- top_by_expression
top_all_df$ensembl_gene_id <- rownames(top_by_expression)

all_candidate_genes <- top_all_df %>%
  left_join(proteins_filtered, by = "ensembl_gene_id")

# Inspect
head(all_candidate_genes)


# Save results
write.csv(top_all, "DE_all_genes_Glioma_vs_GTex.csv", row.names = TRUE)
write.csv(top_sig, "DE_upregulated_genes_Glioma_vs_GTex.csv", row.names = TRUE)
write.csv(
  top_by_expression,
  "ORDERED_DE_upregulated_genes_Glioma_vs_GTex.csv",
  row.names = TRUE
)
write.csv(
  all_candidate_genes,
  "ORDERED_all_candidate_genes.csv",
  row.names = TRUE
)

##### 10. Quick visualizations #####
# Volcano plot
with(
  top_all,
  plot(
    logFC,
    -log10(adj.P.Val + 1e-100),
    pch = 20,
    main = "Volcano: Glioma vs Normal"
  )
)
with(
  subset(top_all, logFC > 1 & adj.P.Val < 0.05),
  points(logFC, -log10(adj.P.Val + 1e-100), col = "red", pch = 20)
)

# Heatmap of top genes (requires pheatmap package)

library(pheatmap)
topn <- rownames(head(top_all[order(top_all$adj.P.Val), ], n = 50))
z <- cpm(dge, log = TRUE)[topn, ]
pheatmap(
  z,
  annotation_col = data.frame(
    Group = meta_comb$group,
    Source = meta_comb$source
  )
)


##### 11. Save R objects for future use #####
saveRDS(dge, file = "dge_combined_filtered.rds")
saveRDS(v, file = "voom_object.rds")
saveRDS(meta_comb, file = "meta_combined.rds")

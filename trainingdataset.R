############################
# DepMap ML Training Matrix
# Using Bioconductor depmap package
############################

# 1. Install depmap package 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("depmap", ask = FALSE)

library(depmap)
library(dplyr)
library(tidyr)

# 2. List available files in the latest release
all_files <- dmfiles()
print(head(all_files$name))

# 3. Select needed files needed from their rows; latest version of these files
latest_rows <- c(
  "CRISPRGeneEffect.csv" = 10,
  "CRISPRGeneDependency.csv" = 9,
  "CCLE_expression.csv" = 239,
  "CCLE_gene_cn.csv" = 244,
  "CCLE_mutations.csv" = 245,
  "sample_info.csv" = 223
)

files_to_download <- all_files[unname(latest_rows), ]


# 4. Download files (cached locally)
local_files <- dmget(files_to_download)

# local_files is a named vector: names are file_name, values are local paths
str(local_files)

#rename the local files with their respective file name 
names(local_files) <- files_to_download$name


# 5. Load files into R
get_file <- function(target, files) {
  files[grep(target, names(files), fixed = TRUE)]
}

crispr_effect     <- readr::read_csv(get_file("CRISPRGeneEffect.csv", local_files))
crispr_dependency <- readr::read_csv(get_file("CRISPRGeneDependency.csv", local_files))
expr              <- readr::read_csv(get_file("CCLE_expression.csv", local_files))
cn                <- readr::read_csv(get_file("CCLE_gene_cn.csv", local_files))
mut               <- readr::read_csv(get_file("CCLE_mutations.csv", local_files))
meta              <- readr::read_csv(get_file("sample_info.csv", local_files))


# 6. Filter for glioma/CNS cell lines
glioma_meta <- meta %>%
  filter(
    # must be CNS
    sample_collection_site == "central_nervous_system",
    grepl("glioma",
          lineage_subtype,
          ignore.case = TRUE)
  )


glioma_ids <- glioma_meta$DepMap_ID

#Filter each using glioma ids and rename ...1 column name to DepMap_ID 
expr_glioma <- expr %>% 
  filter(...1 %in% glioma_ids) %>%
  rename(DepMap_ID = ...1)

cn_glioma <- cn %>% 
  filter(...1 %in% glioma_ids) %>%
  rename(DepMap_ID = ...1)

crispr_effect_glioma <- crispr_effect %>% 
  filter(...1 %in% glioma_ids) %>%
  rename(DepMap_ID = ...1)

mut_glioma  <- mut %>%
  filter(DepMap_ID %in% glioma_ids)


#clean gene names (remove gene id, leave Hugo symbol)
colnames(expr_glioma) <- c("DepMap_ID", sub("\\s*\\(.*", "", colnames(expr_glioma)[-1]))
colnames(cn_glioma) <- c("DepMap_ID", sub("\\s*\\(.*", "", colnames(cn_glioma)[-1]))
colnames(crispr_effect_glioma) <- c("DepMap_ID", sub("\\s*\\(.*", "", colnames(crispr_effect_glioma)[-1]))

# 7. Pivot multi-omics features to long format + Convert CERES scores to binary essentiality 
expr_long <- expr_glioma %>%
  pivot_longer(
    cols = -DepMap_ID,
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(DepMap_ID = trimws(DepMap_ID),
         gene = trimws(gene))

cn_long <- cn_glioma %>%
  pivot_longer(
    cols = -DepMap_ID,
    names_to = "gene",
    values_to = "cn"
  ) %>%
  mutate(DepMap_ID = trimws(DepMap_ID),
         gene = trimws(gene))

mut_long <- mut_glioma %>%
  mutate(mut_flag = 1,
         DepMap_ID = trimws(DepMap_ID),
         gene = trimws(Hugo_Symbol)) %>%
  select(DepMap_ID, gene, mut_flag) %>%
  distinct()

crispr_long <- crispr_effect_glioma %>%
  pivot_longer(
    cols = -DepMap_ID,
    names_to = "gene",
    values_to = "ceres"
  ) %>%
  mutate(
    DepMap_ID = trimws(DepMap_ID),
    gene = trimws(gene),
    essential_binary = if_else(ceres < -0.5, 1, 0)
  )

# 9. Merge all features + labels
training_matrix <- expr_long %>%
  left_join(cn_long, by = c("DepMap_ID", "gene")) %>%
  left_join(mut_long, by = c("DepMap_ID", "gene")) %>%
  mutate(mut_flag = replace_na(mut_flag, 0)) %>%  # Set 0 where gene not mutated
  left_join(crispr_long, by = c("DepMap_ID", "gene"))

# Remove NA targets 
training_final <- training_matrix %>%
  filter(!is.na(ceres), !is.na(essential_binary))

#most frequently mutated genes 
top_mutated_genes <- training_matrix %>%
  group_by(gene) %>%
  summarize(
    mutation_frequency = mean(mut_flag),
    times_mutated = sum(mut_flag)
  ) %>%
  arrange(desc(times_mutated)) %>%
  head(10)

print("Top 10 most frequently mutated genes:")
print(top_mutated_genes)


# See how copy number correlates with essentiality
cn_essential_corr <- training_matrix %>%
  filter(!is.na(essential_binary)) %>%
  group_by(essential_binary) %>%
  summarize(avg_cn = mean(cn, na.rm = TRUE))

print("Average copy number by essentiality:")
print(cn_essential_corr)

# See how gene expression correlates with essentiality
expr_essential_relationship <- training_matrix %>%
  filter(!is.na(essential_binary)) %>%
  group_by(essential_binary) %>%
  summarize(avg_expression = mean(expression, na.rm = TRUE))

print("Expression by essentiality:")
print(expr_essential_relationship)

# See how mutations correlates with essentiality 
mut_essential_relationship <- training_matrix %>%
  filter(!is.na(essential_binary)) %>%
  group_by(essential_binary) %>%
  summarize(mutation_rate = mean(mut_flag, na.rm = TRUE))

print("Mutation rate by essentiality:")
print(mut_essential_relationship)







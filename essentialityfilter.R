# Install if needed and load libraries
install.packages(c("tidyverse", "dplyr", "tidyr", "stringr", "ggplot2", "ggrepel"))


library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)


#Install depmap package for essentiality data 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("depmap", ask = FALSE)

library(depmap)


# Set working directory (adjust to your path)
setwd("C:/Users/jerry/OneDrive/Documents/bioinformatics/Final Project - ML")

#gene_domain_counts from the HMM is the list of proteins with druggable domains
drug_proteins <- read_csv("gene_domain_counts.csv")

#rename gene_ID column to protein_ID
drug_proteins <- drug_proteins %>% 
  rename(Protein_id = Gene_ID)






## Step 1. Creating new column for HUGO gene symbol -> used for matching to protein to gene
BiocManager::install("UniProt.ws")
library(UniProt.ws)

#set taxon to human 
up <- UniProt.ws(taxId = 9606)

#map protein ID to HUGO gene ID 
mapping <- select(
  up,
  keys = drug_proteins$Protein_id,
  columns = c("UniProtKB", "gene_primary"),
  keytype = "UniProtKB"
)

#add new column for gene ID 
drug_proteins <- drug_proteins %>%
  left_join(mapping, by = c("Protein_id" = "From"))

drug_proteins <- drug_proteins %>%
  rename(Gene_ID = Gene.Names..primary.)

## Step 2. Summarize ceres score
depmap_glioma <- read_csv("training.csv")

#summarize ceres score and arrange from most negative to least  
ceres_summary <- depmap_glioma %>%
  group_by(gene) %>%
  summarise(
    mean_ceres = mean(ceres, na.rm = TRUE),
    median_ceres = median(ceres, na.rm = TRUE),
    min_ceres = min(ceres, na.rm = TRUE),      # Most negative (most essential)
    max_ceres = max(ceres, na.rm = TRUE),
    sd_ceres = sd(ceres, na.rm = TRUE),
    n_cell_lines = n(),
    # Optional: Calculate essentiality confidence
    essential_fraction = sum(ceres < -0.5) / n()  # Fraction of cell lines where gene is essential
  ) %>%
  arrange(mean_ceres)

#merge 
merged_data <- drug_proteins %>%
  left_join(ceres_summary, by = c("Gene_ID" = "gene"))

#filter out genes without ceres score 
merged_data <- merged_data %>%
  filter(!is.na(mean_ceres))

## Step 3. Rank genes based on therapeutic potential
# Primary ranking by mean CERES (most negative first)
ranked_genes <- genes_with_ceres %>%
  arrange(mean_ceres) %>%
  mutate(
    rank_by_mean = row_number(),
    # Create therapeutic score (higher = better)
    therapeutic_score = -mean_ceres,  # Negative CERES becomes positive score
    # Confidence score combining multiple metrics
    confidence_score = (-mean_ceres) * (1 - sd_ceres) * essential_fraction
  ) %>%
  arrange(desc(therapeutic_score)) %>%
  mutate(therapeutic_rank = row_number())

# Step 4: Create a comprehensive therapeutic index
final_ranking <- genes_with_ceres %>%
  mutate(
    # Normalize metrics (optional)
    ceres_norm = scale(-mean_ceres),  # More negative = higher score
    consistency_norm = scale(1/sd_ceres),  # Lower SD = higher consistency
    coverage_norm = scale(n_cell_lines),  # More cell lines = better coverage
    
    # Create composite therapeutic index
    # Adjust weights based on priorities:
    therapeutic_index = (ceres_norm * 0.6) +  # Weight CERES most heavily
      (consistency_norm * 0.3) +  # Consistency across cell lines
      (coverage_norm * 0.1),  # Number of cell lines
    
    # Alternative simple index
    simple_therapeutic_index = -mean_ceres * essential_fraction
  ) %>%
  arrange(desc(therapeutic_index)) %>%
  mutate(
    final_rank = row_number(),
    percentile_rank = round(rank(therapeutic_index) / n() * 100, 1)
  )

# Step 6: View top therapeutic targets
top_therapeutic_targets <- final_ranking[
  order(final_ranking$final_rank),
  c("Protein_id", "Found_Domains", "Gene_ID", "mean_ceres", "therapeutic_index", "final_rank")
]
top_therapeutic_targets <- head(top_therapeutic_targets, 20)

print("Top 20 Therapeutic Targets for Glioma:")
print(top_therapeutic_targets)


#View Top Glioma Targets 
top_20_plot <- final_ranking %>%
  head(20)

ggplot(top_20_plot, aes(x = reorder(Gene_ID, therapeutic_index), 
                        y = therapeutic_index, 
                        fill = mean_ceres)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "darkred", mid = "white", high = "steelblue",
                       midpoint = 0, name = "Mean CERES") +
  coord_flip() +
  labs(
    title = "Top 20 Glioma Therapeutic Targets",
    subtitle = "Higher Therapeutic Index = Better Drug Target",
    x = "Gene Symbol",
    y = "Therapeutic Index"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )


library(tidyverse)

# 1. Load the sorted data (if not already in your environment)
candidates <- read_csv("sorted_druggable_genes.csv")

# 2. Group by Gene and Count
gene_summary <- candidates %>%
    group_by(Gene_ID) %>%
    summarise(
        # Count how many rows exist for this gene
        Domain_Count = n(),

        # Concatenate the domain names into a single string for easy reading
        # e.g., "Pkinase; Pkinase; SH2"
        Found_Domains = paste(Domain_Name, collapse = "; "),

        # Optional: Keep the best E-value found for this gene
        Best_E_Value = min(E_Value)
    ) %>%

    # 3. Sort by Count (High to Low) to see most "domain-rich" genes first
    arrange(desc(Domain_Count))

# 4. View the top results in R
print(head(gene_summary))

# 5. Save this summary to a new CSV
write_csv(gene_summary, "gene_domain_counts.csv")

# quick visualization
library(ggplot2)

p <- ggplot(gene_summary, aes(x = factor(Domain_Count))) +
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    labs(
        title = "Distribution of Druggable Domains per Gene",
        x = "Number of Domains Found",
        y = "Number of Genes"
    )

# Save the plot
ggsave("domain_distribution.png", plot = p, width = 6, height = 4, dpi = 300)

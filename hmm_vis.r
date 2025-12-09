library(ggplot2)

ggplot(gene_summary, aes(x = factor(Domain_Count))) +
    geom_bar(fill = "steelblue") +
    theme_minimal() +
    labs(
        title = "Distribution of Druggable Domains per Gene",
        x = "Number of Domains Found",
        y = "Number of Genes"
    )

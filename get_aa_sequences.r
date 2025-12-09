library(UniprotR)
library(seqinr)

# 1. Read the CSV file
df <- read.csv("ORDERED_all_candidate_genes.csv", stringsAsFactors = FALSE)

# 2. Clean the IDs
ids <- trimws(df$uniprotswissprot)

# Optional: Remove any NA or empty entries to prevent errors
ids <- ids[ids != "" & !is.na(ids)]

# 3. Retrieve the sequences
# Connect to UniProt and fetches data for the list of IDs
print(paste("Fetching sequences for", length(ids), "proteins..."))
protein_data <- GetSequences(ids)

# 4. Write to FASTA format
# The output of GetSequences is a dataframe. We need to parse it for seqinr.
# Note: protein_data row names are the Accession IDs, and the column 'Sequence' holds the AA string.

write.fasta(
    sequences = as.list(protein_data$Sequence),
    names = rownames(protein_data),
    file.out = "output_sequences.fasta",
    open = "w"
)

print("Success! File saved as output_sequences.fasta")

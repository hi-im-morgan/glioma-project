library(tidyverse)

sort_and_save_hmmer <- function(input_file, output_file) {
    # 1. Define Column Names
    hmmer_cols <- c(
        "target_name",
        "target_accession",
        "tlen",
        "query_name",
        "query_accession",
        "qlen",
        "full_evalue",
        "full_score",
        "full_bias",
        "domain_number",
        "total_domains",
        "c_evalue",
        "i_evalue",
        "score_dom",
        "bias_dom",
        "from_hmm",
        "to_hmm",
        "from_ali",
        "to_ali",
        "from_env",
        "to_env",
        "acc",
        "description"
    )

    # 2. Read and Parse (Robust method using separate)
    raw_lines <- read_lines(input_file)

    clean_data <- tibble(raw = raw_lines) %>%
        # Remove comments
        filter(!str_detect(raw, "^#")) %>%
        # Split columns, merging extra description text
        separate(raw, into = hmmer_cols, sep = "\\s+", extra = "merge") %>%

        # 3. CRITICAL: Convert E-value to numeric for correct sorting
        # If you skip this, R sorts alphabetically (e.g. "10" comes before "2")
        mutate(
            i_evalue = as.numeric(i_evalue),
            from_ali = as.numeric(from_ali),
            to_ali = as.numeric(to_ali)
        ) %>%

        # 4. Filter
        filter(i_evalue < 0.01) %>%
        filter(grepl(
            "Pkinase|7tm|Ion_trans|Nuclear_rec|Trypsin|ABC_tran",
            target_name,
            ignore.case = TRUE
        )) %>%

        # 5. Select Columns
        select(
            Gene_ID = query_name,
            Domain_Name = target_name,
            Description = description,
            E_Value = i_evalue,
            Start = from_ali,
            End = to_ali
        ) %>%

        # 6. SORT: Ascending E-value (Lowest/Best first)
        arrange(E_Value)

    # 7. Write to CSV
    write_csv(clean_data, output_file)

    message(paste(
        "Sorted and saved",
        nrow(clean_data),
        "domains to",
        output_file
    ))
    return(clean_data)
}

# --- USAGE ---
# Replace with your actual file names
data <- sort_and_save_hmmer("hmm_results.txt", "sorted_druggable_genes.csv")

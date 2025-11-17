#!/usr/bin/env Rscript

safely_source <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Cannot source required file '%s'", path))
    }
    source(path, local = FALSE)
}

safely_source("src/hmm.r")

hmm_class <- get("ProteinDomainHMM", inherits = TRUE)

sequence <- "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDD"

hmm <- hmm_class$new(
    n_states = 3,
    min_domain = 8,
    verbose = FALSE
)

result <- hmm$identify_domains(sequence)

stopifnot(
    is.list(result),
    "domains" %in% names(result),
    "loglik" %in% names(result),
    is.numeric(result$loglik)
)

if (!is.null(result$domains) && nrow(result$domains) > 0) {
    cat(sprintf("Identified %d domain candidate(s).\n", nrow(result$domains)))
} else {
    cat("No domain candidates detected for the smoke-test sequence.\n")
}

cat(sprintf("logLik: %.3f\n", result$loglik))

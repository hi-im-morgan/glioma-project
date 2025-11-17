#!/usr/bin/env Rscript

# Minimal smoke test for ProteinDomainHMM. Intended to run quickly with synthetic input.

get_script_dir <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_flag <- "--file="
    script_path <- NULL
    matches <- grep(file_flag, cmd_args)
    if (length(matches)) {
        script_path <- sub(file_flag, "", cmd_args[matches[length(matches)]])
    } else if (!is.null(sys.frames()[[1]])$ofile) {
        script_path <- sys.frames()[[1]]$ofile
    }

    if (is.null(script_path)) {
        return(normalizePath(".", winslash = "/", mustWork = FALSE))
    }
    normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
}

PROJECT_ROOT <- normalizePath(
    file.path(get_script_dir(), ".."),
    winslash = "/",
    mustWork = FALSE
)

safely_source <- function(relative_path) {
    target <- file.path(PROJECT_ROOT, relative_path)
    if (!file.exists(target)) {
        stop(sprintf("Required file '%s' is missing.", target))
    }
    source(target, local = FALSE)
}

safely_source("src/motif_hmm.r")
safely_source("src/hmm.r")

if (!exists("ProteinDomainHMM", inherits = TRUE)) {
    stop("ProteinDomainHMM is unavailable after sourcing src/hmm.r")
}
ProteinDomainHMM <- get("ProteinDomainHMM", inherits = TRUE)
if (!exists("ProteinMotifDetector", inherits = TRUE)) {
    stop("ProteinMotifDetector is unavailable after sourcing src/motif_hmm.r")
}
ProteinMotifDetector <- get("ProteinMotifDetector", inherits = TRUE)

parse_args <- function() {
    defaults <- list(
        sequence = "MEEPQSDPSVEPPLSQETF",
        n_states = 3L,
        min_domain = 8L,
        verbose = FALSE,
        output = "output/tests/hmm_smoke_test.rds"
    )

    args <- commandArgs(trailingOnly = TRUE)
    if (!length(args)) {
        return(defaults)
    }

    opts <- defaults
    for (arg in args) {
        if (!grepl("^--[A-Za-z0-9_\\-]+=", arg)) {
            warning(sprintf("Ignoring malformed argument: %s", arg))
            next
        }
        kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
        key <- kv[1]
        value <- kv[2]

        if (identical(key, "sequence")) {
            opts$sequence <- value
        } else if (identical(key, "n_states")) {
            opts$n_states <- as.integer(value)
        } else if (identical(key, "min_domain")) {
            opts$min_domain <- as.integer(value)
        } else if (identical(key, "verbose")) {
            opts$verbose <- tolower(value) %in% c("1", "true", "t", "yes", "y")
        } else if (identical(key, "output")) {
            opts$output <- value
        } else {
            warning(sprintf("Unknown argument '%s'; ignoring.", key))
        }
    }

    if (is.na(opts$n_states) || opts$n_states < 1) {
        stop("`n_states` must be an integer ≥ 1.")
    }
    if (is.na(opts$min_domain) || opts$min_domain < 1) {
        stop("`min_domain` must be an integer ≥ 1.")
    }

    opts
}

main <- function() {
    opts <- parse_args()

    hmm <- ProteinDomainHMM$new(
        n_states = opts$n_states,
        min_domain = opts$min_domain,
        verbose = opts$verbose
    )

    if (opts$verbose) {
        message(sprintf(
            "Smoke test: sequence length=%d, states=%d",
            nchar(opts$sequence),
            opts$n_states
        ))
    }

    res <- hmm$identify_domains(opts$sequence)
    motif_summary <- hmm$identify_motifs(
        sequence = opts$sequence,
        motifs = c(
            activation = "MEEPQ",
            dna_binding = "SQETF"
        ),
        probability_threshold = 0.6
    )

    detector <- ProteinMotifDetector$new(motif_threshold = 0.6)
    standalone_motifs <- detector$find_motifs(
        sequence = opts$sequence,
        motifs = c(
            activation = "MEEPQ",
            dna_binding = "SQETF"
        ),
        probability_threshold = 0.6
    )

    abs_output <- if (grepl("^/", opts$output)) {
        opts$output
    } else {
        file.path(PROJECT_ROOT, opts$output)
    }
    dir.create(dirname(abs_output), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
        list(
            domains = res,
            motifs = motif_summary,
            standalone_motifs = standalone_motifs
        ),
        abs_output
    )

    message(sprintf(
        paste0(
            "Smoke test OK: %d domain(s), %d motif hit(s) detected (logLik=%.3f). ",
            "Output saved to %s"
        ),
        if (!is.null(res$domains)) nrow(res$domains) else 0,
        if (!is.null(motif_summary$hits)) nrow(motif_summary$hits) else 0,
        res$loglik,
        abs_output
    ))
}

if (sys.nframe() == 0) {
    main()
}

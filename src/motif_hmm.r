# motif_hmm.r
if (!requireNamespace("R6", quietly = TRUE)) {
    stop("Install R6 to use ProteinMotifDetector.")
}

ProteinMotifDetector <- R6::R6Class(
    "ProteinMotifDetector",
    private = list(
        aa_alphabet = c(
            "A",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "I",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "V",
            "W",
            "Y"
        ),
        normalize_vec = function(v) v / sum(v),
        normalize_rows = function(mat) sweep(mat, 1, rowSums(mat), "/"),
        log_sum_exp = function(x) {
            m <- max(x)
            m + log(sum(exp(x - m)))
        },
        encode_sequence = function(sequence) {
            aa_index <- setNames(
                seq_along(private$aa_alphabet),
                private$aa_alphabet
            )
            obs <- aa_index[toupper(strsplit(sequence, "")[[1]])]
            if (any(is.na(obs))) {
                stop("Sequence contains unsupported symbols.")
            }
            obs
        },
        forward_backward = function(obs, model) {
            Tlen <- length(obs)
            N <- length(model$pi)
            logA <- log(model$A)
            logB <- log(model$B)
            logPi <- log(model$pi)

            alpha <- matrix(-Inf, Tlen, N)
            alpha[1, ] <- logPi + logB[, obs[1]]
            for (t in 2:Tlen) {
                for (j in 1:N) {
                    alpha[t, j] <- logB[j, obs[t]] +
                        private$log_sum_exp(alpha[t - 1, ] + logA[, j])
                }
            }

            beta <- matrix(-Inf, Tlen, N)
            beta[Tlen, ] <- 0
            for (t in seq(Tlen - 1, 1)) {
                for (i in 1:N) {
                    beta[t, i] <- private$log_sum_exp(
                        logA[i, ] + logB[, obs[t + 1]] + beta[t + 1, ]
                    )
                }
            }

            loglik <- private$log_sum_exp(alpha[Tlen, ])
            gamma <- exp(alpha + beta - loglik)

            list(loglik = loglik, gamma = gamma)
        },
        prepare_motifs = function(motifs) {
            if (is.null(motifs) || !length(motifs)) {
                stop(
                    "Provide at least one motif (character vector or named list)."
                )
            }

            if (is.list(motifs)) {
                motifs <- unlist(motifs, use.names = TRUE)
            }

            if (!is.character(motifs)) {
                stop("Motifs must be supplied as character strings.")
            }

            if (is.null(names(motifs)) || any(names(motifs) == "")) {
                names(motifs) <- sprintf("motif_%02d", seq_along(motifs))
            }

            alphabet_pattern <- paste(private$aa_alphabet, collapse = "")
            invalid <- grepl(
                sprintf("[^%s]", alphabet_pattern),
                toupper(motifs)
            )
            if (any(invalid)) {
                bad <- paste(unique(motifs[invalid]), collapse = ", ")
                stop(sprintf("Motifs contain unsupported symbols: %s", bad))
            }

            lapply(motifs, function(m) {
                strsplit(toupper(m), "", fixed = TRUE)[[1]]
            })
        },
        build_motif_hmm = function(motif_chars, entry_prob, dropout_prob) {
            entry_prob <- min(max(entry_prob, 1e-4), 0.5)
            dropout_prob <- min(max(dropout_prob, 0), 0.5)

            L <- length(motif_chars)
            if (L < 1) {
                stop("Motif length must be at least 1.")
            }

            n_symbols <- length(private$aa_alphabet)
            n_states <- L + 1

            pi <- numeric(n_states)
            pi[1] <- 1

            A <- matrix(0, n_states, n_states)
            A[1, 1] <- 1 - entry_prob
            A[1, 2] <- entry_prob
            if (L == 1) {
                A[2, 2] <- dropout_prob
                A[2, 1] <- 1 - dropout_prob
            } else {
                for (i in 2:L) {
                    A[i, i + 1] <- 1 - dropout_prob
                    A[i, 1] <- dropout_prob
                }
                A[n_states, n_states] <- dropout_prob
                A[n_states, 1] <- 1 - dropout_prob
            }

            B <- matrix(0, n_states, n_symbols)
            B[1, ] <- rep(1 / n_symbols, n_symbols)

            for (i in seq_len(L)) {
                emissions <- rep(self$motif_pseudocount, n_symbols)
                idx <- match(motif_chars[i], private$aa_alphabet)
                emissions[idx] <- emissions[idx] + 1
                B[i + 1, ] <- emissions / sum(emissions)
            }

            list(pi = pi, A = A, B = B)
        },
        extract_motif_segments = function(
            prob_series,
            threshold,
            min_len,
            max_len,
            seq_chars,
            motif_name
        ) {
            if (!length(prob_series)) {
                return(list())
            }

            if (is.null(min_len) || min_len < 1) {
                min_len <- 1
            }

            if (is.null(max_len) || !is.finite(max_len)) {
                max_len <- length(prob_series)
            }

            max_len <- max(max_len, min_len)
            hits <- list()
            cursor <- 1L
            runs <- rle(prob_series >= threshold)

            for (i in seq_along(runs$lengths)) {
                len <- runs$lengths[i]
                val <- runs$values[i]
                start_idx <- cursor
                end_idx <- cursor + len - 1L
                if (val && len >= min_len && len <= max_len) {
                    subseq <- paste(seq_chars[start_idx:end_idx], collapse = "")
                    segment_probs <- prob_series[start_idx:end_idx]
                    hits[[length(hits) + 1]] <- list(
                        motif = motif_name,
                        start = start_idx,
                        end = end_idx,
                        length = len,
                        mean_probability = mean(segment_probs),
                        peak_probability = max(segment_probs),
                        sequence = subseq
                    )
                }
                cursor <- cursor + len
            }

            hits
        }
    ),
    public = list(
        motif_entry_prob = NULL,
        motif_dropout_prob = NULL,
        motif_min_length = NULL,
        motif_max_length = NULL,
        motif_threshold = NULL,
        motif_pseudocount = NULL,
        last_hits = NULL,
        last_probabilities = NULL,

        initialize = function(
            alphabet = NULL,
            motif_entry_prob = 0.05,
            motif_dropout_prob = 0.05,
            motif_min_length = 3,
            motif_max_length = 15,
            motif_threshold = 0.65,
            motif_pseudocount = 1e-3
        ) {
            if (!is.null(alphabet)) {
                private$aa_alphabet <- toupper(alphabet)
            }
            self$motif_entry_prob <- motif_entry_prob
            self$motif_dropout_prob <- motif_dropout_prob
            self$motif_min_length <- motif_min_length
            self$motif_max_length <- motif_max_length
            self$motif_threshold <- motif_threshold
            self$motif_pseudocount <- motif_pseudocount
        },

        find_motifs = function(
            sequence,
            motifs,
            probability_threshold = NULL,
            min_length = NULL,
            max_length = NULL,
            entry_prob = NULL,
            dropout_prob = NULL
        ) {
            if (missing(sequence) || !nchar(sequence)) {
                stop("`sequence` must be a non-empty character string.")
            }

            motif_definitions <- private$prepare_motifs(motifs)
            seq_chars <- strsplit(toupper(sequence), "", fixed = TRUE)[[1]]
            obs <- private$encode_sequence(sequence)

            threshold <- if (is.null(probability_threshold)) {
                self$motif_threshold
            } else {
                probability_threshold
            }

            base_min_len <- if (is.null(min_length)) {
                self$motif_min_length
            } else {
                min_length
            }
            base_max_len <- if (is.null(max_length)) {
                self$motif_max_length
            } else {
                max_length
            }
            entry <- if (is.null(entry_prob)) {
                self$motif_entry_prob
            } else {
                entry_prob
            }
            dropout <- if (is.null(dropout_prob)) {
                self$motif_dropout_prob
            } else {
                dropout_prob
            }

            all_hits <- list()
            motif_summaries <- list()

            for (motif_name in names(motif_definitions)) {
                motif_chars <- motif_definitions[[motif_name]]
                model <- private$build_motif_hmm(motif_chars, entry, dropout)
                fb <- private$forward_backward(obs, model)
                motif_states <- 2:length(model$pi)
                motif_prob <- rowSums(fb$gamma[, motif_states, drop = FALSE])

                local_min <- if (is.null(base_min_len)) {
                    length(motif_chars)
                } else {
                    max(1, base_min_len)
                }
                local_min <- max(local_min, length(motif_chars))

                local_max <- if (is.null(base_max_len)) {
                    length(seq_chars)
                } else {
                    max(local_min, base_max_len)
                }

                hits <- private$extract_motif_segments(
                    prob_series = motif_prob,
                    threshold = threshold,
                    min_len = local_min,
                    max_len = local_max,
                    seq_chars = seq_chars,
                    motif_name = motif_name
                )

                if (length(hits)) {
                    all_hits <- c(all_hits, hits)
                }

                motif_summaries[[motif_name]] <- list(
                    loglik = fb$loglik,
                    probability = motif_prob
                )
            }

            hits_df <- if (length(all_hits)) {
                do.call(rbind, lapply(all_hits, as.data.frame))
            } else {
                data.frame()
            }

            self$last_hits <- hits_df
            self$last_probabilities <- motif_summaries

            list(
                hits = hits_df,
                per_motif = motif_summaries
            )
        }
    )
)

# Example usage:
# source("src/motif_hmm.r")
# detector <- ProteinMotifDetector$new(motif_threshold = 0.7)
# detector$find_motifs(
#     sequence = "MEEPQSDPSVEPPLSQETF",
#     motifs = c(activation = "MEEPQ", dna_binding = "SQETF")
# )

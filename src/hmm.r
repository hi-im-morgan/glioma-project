# hmm.r
if (!requireNamespace("R6", quietly = TRUE)) {
    stop("Install R6 to use ProteinDomainHMM.")
}

ProteinDomainHMM <- R6::R6Class(
    "ProteinDomainHMM",
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
        init_hmm = function(n_states, n_symbols) {
            list(
                pi = private$normalize_vec(runif(n_states)),
                A = private$normalize_rows(matrix(runif(n_states^2), n_states)),
                B = private$normalize_rows(matrix(
                    runif(n_states * n_symbols),
                    n_states
                ))
            )
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

            xi <- array(0, c(Tlen - 1, N, N))
            for (t in 1:(Tlen - 1)) {
                for (i in 1:N) {
                    for (j in 1:N) {
                        xi[t, i, j] <- alpha[t, i] +
                            logA[i, j] +
                            logB[j, obs[t + 1]] +
                            beta[t + 1, j] -
                            loglik
                    }
                }
            }
            xi <- exp(xi)

            list(loglik = loglik, gamma = gamma, xi = xi)
        },
        baum_welch = function(obs) {
            model <- private$init_hmm(
                self$n_states,
                length(private$aa_alphabet)
            )
            prev_ll <- -Inf

            for (iter in seq_len(self$max_iter)) {
                fb <- private$forward_backward(obs, model)
                gamma <- fb$gamma
                xi_sum <- apply(fb$xi, c(2, 3), sum)

                model$pi <- private$normalize_vec(gamma[1, ] + self$pseudocount)
                A_num <- xi_sum + self$pseudocount
                model$A <- private$normalize_rows(A_num)

                B_num <- matrix(
                    self$pseudocount,
                    self$n_states,
                    length(private$aa_alphabet)
                )
                for (k in seq_along(private$aa_alphabet)) {
                    mask <- obs == k
                    if (any(mask)) {
                        B_num[, k] <- B_num[, k] +
                            colSums(gamma[mask, , drop = FALSE])
                    }
                }
                model$B <- private$normalize_rows(B_num)

                if (self$verbose) {
                    message(sprintf("Iter %03d logLik=%.4f", iter, fb$loglik))
                }
                if (abs(fb$loglik - prev_ll) < self$tol) {
                    return(list(model = model, loglik = fb$loglik))
                }
                prev_ll <- fb$loglik
            }
            list(model = model, loglik = prev_ll)
        },
        viterbi = function(obs, model) {
            Tlen <- length(obs)
            N <- length(model$pi)
            logA <- log(model$A)
            logB <- log(model$B)
            logPi <- log(model$pi)

            delta <- matrix(-Inf, Tlen, N)
            psi <- matrix(0L, Tlen, N)
            delta[1, ] <- logPi + logB[, obs[1]]

            for (t in 2:Tlen) {
                for (j in 1:N) {
                    scores <- delta[t - 1, ] + logA[, j]
                    psi[t, j] <- which.max(scores)
                    delta[t, j] <- max(scores) + logB[j, obs[t]]
                }
            }

            path <- integer(Tlen)
            path[Tlen] <- which.max(delta[Tlen, ])
            for (t in seq(Tlen - 1, 1)) {
                path[t] <- psi[t + 1, path[t + 1]]
            }
            list(path = path, loglik = max(delta[Tlen, ]))
        },
        extract_domains = function(path, posterior) {
            persistence <- rowSums(posterior)
            state_scores <- persistence
            domain_states <- which(state_scores >= quantile(state_scores, 0.75))
            if (!length(domain_states)) {
                domain_states <- which.max(state_scores)
            }

            runs <- rle(path %in% domain_states)
            ends <- cumsum(runs$lengths)
            starts <- ends - runs$lengths + 1

            domains <- Map(
                function(start, end, is_domain) {
                    if (!is_domain) {
                        return(NULL)
                    }
                    if ((end - start + 1) < self$min_domain) {
                        return(NULL)
                    }
                    segment_states <- path[start:end]
                    dominant <- as.integer(names(which.max(table(
                        segment_states
                    ))))
                    score <- mean(posterior[cbind(start:end, dominant)])
                    list(
                        start = start,
                        end = end,
                        length = end - start + 1,
                        state = dominant,
                        confidence = score
                    )
                },
                starts,
                ends,
                runs$values
            )

            domains <- Filter(Negate(is.null), domains)
            if (!length(domains)) {
                data.frame()
            } else {
                do.call(rbind, lapply(domains, as.data.frame))
            }
        },
        ensure_motif_detector = function() {
            if (!exists("ProteinMotifDetector", inherits = TRUE)) {
                stop(
                    paste(
                        "ProteinMotifDetector is unavailable.",
                        "Source('src/motif_hmm.r') before calling identify_motifs."
                    )
                )
            }
            TRUE
        },
        new_motif_detector = function(
            entry_prob,
            dropout_prob,
            min_length,
            max_length,
            threshold,
            pseudocount
        ) {
            private$ensure_motif_detector()
            ProteinMotifDetector$new(
                alphabet = private$aa_alphabet,
                motif_entry_prob = entry_prob,
                motif_dropout_prob = dropout_prob,
                motif_min_length = min_length,
                motif_max_length = max_length,
                motif_threshold = threshold,
                motif_pseudocount = pseudocount
            )
        }
    ),
    public = list(
        n_states = NULL,
        min_domain = NULL,
        max_iter = NULL,
        tol = NULL,
        pseudocount = NULL,
        verbose = NULL,
        motif_entry_prob = NULL,
        motif_dropout_prob = NULL,
        motif_min_length = NULL,
        motif_max_length = NULL,
        motif_threshold = NULL,
        motif_pseudocount = NULL,
        model = NULL,
        loglik = NULL,
        viterbi_path = NULL,
        posterior = NULL,
        domains = NULL,
        motif_hits = NULL,

        initialize = function(
            n_states = 3,
            min_domain = 25,
            max_iter = 200,
            tol = 1e-4,
            pseudocount = 1e-6,
            verbose = FALSE,
            motif_entry_prob = 0.05,
            motif_dropout_prob = 0.05,
            motif_min_length = 3,
            motif_max_length = 15,
            motif_threshold = 0.65,
            motif_pseudocount = 1e-3
        ) {
            self$n_states <- n_states
            self$min_domain <- min_domain
            self$max_iter <- max_iter
            self$tol <- tol
            self$pseudocount <- pseudocount
            self$verbose <- verbose
            self$motif_entry_prob <- motif_entry_prob
            self$motif_dropout_prob <- motif_dropout_prob
            self$motif_min_length <- motif_min_length
            self$motif_max_length <- motif_max_length
            self$motif_threshold <- motif_threshold
            self$motif_pseudocount <- motif_pseudocount
        },

        identify_domains = function(sequence) {
            obs <- private$encode_sequence(sequence)
            fit <- private$baum_welch(obs)
            fb <- private$forward_backward(obs, fit$model)
            vit <- private$viterbi(obs, fit$model)

            self$model <- fit$model
            self$loglik <- fit$loglik
            self$viterbi_path <- vit$path
            self$posterior <- fb$gamma
            self$domains <- private$extract_domains(vit$path, fb$gamma)

            list(
                model = self$model,
                loglik = self$loglik,
                viterbi_path = self$viterbi_path,
                posterior = self$posterior,
                domains = self$domains
            )
        },

        identify_motifs = function(
            sequence,
            motifs,
            probability_threshold = NULL,
            min_length = NULL,
            max_length = NULL,
            entry_prob = NULL,
            dropout_prob = NULL
        ) {
            detector <- private$new_motif_detector(
                entry_prob = if (is.null(entry_prob)) {
                    self$motif_entry_prob
                } else {
                    entry_prob
                },
                dropout_prob = if (is.null(dropout_prob)) {
                    self$motif_dropout_prob
                } else {
                    dropout_prob
                },
                min_length = if (is.null(min_length)) {
                    self$motif_min_length
                } else {
                    min_length
                },
                max_length = if (is.null(max_length)) {
                    self$motif_max_length
                } else {
                    max_length
                },
                threshold = if (is.null(probability_threshold)) {
                    self$motif_threshold
                } else {
                    probability_threshold
                },
                pseudocount = self$motif_pseudocount
            )

            result <- detector$find_motifs(
                sequence = sequence,
                motifs = motifs,
                probability_threshold = probability_threshold,
                min_length = min_length,
                max_length = max_length,
                entry_prob = entry_prob,
                dropout_prob = dropout_prob
            )

            self$motif_hits <- result$hits
            result
        }
    )
)

# Example main script usage:
# source("hmm.r")
# hmm <- ProteinDomainHMM$new(n_states = 4, min_domain = 10, verbose = TRUE)
# result <- hmm$identify_domains("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDD")
# print(result$domains)

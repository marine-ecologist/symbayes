# ==============================================================================
# fit-lda.R — Latent Dirichlet Allocation (fractional membership)
# ==============================================================================

#' Fit a Latent Dirichlet Allocation topic model
#'
#' Fits LDA by Gibbs sampling. Unlike DMM, LDA gives *fractional* membership:
#' each sample is a mixture of topics, so within-sample admixture is
#' represented directly.
#'
#' Topic number selection is threshold-free-oriented: with `k = NULL` (default)
#' the function searches `k_range` and picks the **largest k with no redundant
#' topic pairs** (cosine similarity `< redundant_thresh`). This avoids the
#' "shadow topics" that perplexity-optimal K tends to produce. Pass an integer
#' `k` to fix it. The search table is stored in `sp$lda$meta$selection` and can
#' be plotted with [lda_plot_ksearch].
#'
#' Writes a membership object to `sp$lda`. `meta` contains `entropy`,
#' `max_prop`, `n_topics` (per sample), `similarity` (topic cosine matrix),
#' `redundant` (flagged pairs), and `selection` (the K-search table).
#'
#' @param sp A filtered `symbayes` object.
#' @param k Fixed number of topics, or `NULL` (default) to select adaptively.
#' @param k_range Integer vector of K to search when `k = NULL` (default 2:15).
#' @param redundant_thresh Cosine similarity above which two topics are deemed
#'   redundant (default 0.85).
#' @param burnin,iter,thin Gibbs sampler controls (defaults 1000, 2000, 5).
#' @param seed Random seed (default 42).
#' @return The `symbayes` object with `sp$lda` populated.
#' @examples
#' \dontrun{
#' sp <- fit_lda(sp)              # adaptive K
#' sp <- fit_lda(sp, k = 9)       # fixed K
#' lda_plot_ksearch(sp)
#' lda_plot_similarity(sp)
#' }
#' @export
fit_lda <- function(sp, k = NULL, k_range = 2:15,
                    redundant_thresh = 0.85,
                    burnin = 1000, iter = 2000, thin = 5, seed = 42) {
  .check_sp(sp)
  if (!isTRUE(sp$filtered))
    warning("Fitting on unfiltered data; consider filter_samples() first.",
            call. = FALSE)

  dtm <- .count_to_dtm(sp$count_mat)
  control <- list(seed = seed, burnin = burnin, iter = iter, thin = thin)

  selection <- NULL
  if (is.null(k)) {
    message(sprintf("Selecting LDA K by redundancy (cosine < %.2f) over %d:%d ...",
                    redundant_thresh, min(k_range), max(k_range)))
    selection <- .lda_search(dtm, k_range, redundant_thresh, control)
    valid <- selection[selection$n_redundant == 0, , drop = FALSE]
    k <- if (nrow(valid) > 0) max(valid$K) else
      selection$K[which.min(selection$n_redundant)]
    message(sprintf("  selected K = %d", k))
  }

  fit <- topicmodels::LDA(dtm, k = k, method = "Gibbs", control = control)
  post  <- topicmodels::posterior(fit)
  theta <- post$topics          # samples x topics
  beta  <- post$terms           # topics x sequences
  colnames(theta) <- paste0("Topic_", seq_len(k))
  rownames(beta)  <- paste0("Topic_", seq_len(k))

  dominant <- setNames(apply(theta, 1, which.max), rownames(theta))
  entropy  <- apply(theta, 1, .entropy)
  max_prop <- apply(theta, 1, max)
  n_topics <- apply(theta, 1, function(p) sum(p > 0.05))

  # Topic cosine similarity
  sim <- .cosine_matrix(beta)
  redundant <- .redundant_pairs(sim, redundant_thresh)

  md <- sp$metadata
  md$lda_dominant <- factor(dominant[md$sample_uid])
  md$lda_entropy  <- entropy[md$sample_uid]
  md$lda_max_prop <- max_prop[md$sample_uid]
  md$lda_n_topics <- n_topics[md$sample_uid]
  sp$metadata <- .map_symportal(sp, md, "lda")

  sp$lda <- .membership(
    theta = theta, beta = beta,
    dominant = factor(dominant), k = k,
    model = "lda", label_prefix = "Topic",
    meta = list(entropy = entropy, max_prop = max_prop, n_topics = n_topics,
                similarity = sim, redundant = redundant,
                selection = selection, fit = fit, dtm = dtm)
  )

  message(sprintf("  LDA: K = %d; %.0f%% samples mixed (>1 topic >5%%); mean entropy %.2f",
                  k, 100 * mean(n_topics > 1), mean(entropy)))
  sp
}


#' LDA K-search over a range, scoring redundancy and pure samples
#' @keywords internal
.lda_search <- function(dtm, k_range, redundant_thresh, control) {
  rows <- lapply(k_range, function(kk) {
    fit <- topicmodels::LDA(dtm, k = kk, method = "Gibbs", control = control)
    beta  <- topicmodels::posterior(fit)$terms
    theta <- topicmodels::posterior(fit)$topics
    sim <- .cosine_matrix(beta)
    max_sim <- max(sim[upper.tri(sim)])
    n_red   <- sum(sim[upper.tri(sim)] > redundant_thresh)
    n_pure  <- sum(vapply(seq_len(kk), function(j) any(theta[, j] > 0.9), logical(1)))
    perp    <- topicmodels::perplexity(fit, newdata = dtm)
    data.frame(K = kk, max_cosine = max_sim, n_redundant = n_red,
               n_pure_topics = n_pure, pct_pure = 100 * n_pure / kk,
               perplexity = perp)
  })
  do.call(rbind, rows)
}

#' Cosine similarity matrix among rows of a matrix
#' @keywords internal
.cosine_matrix <- function(m) {
  norms <- sqrt(rowSums(m^2))
  sim <- (m %*% t(m)) / (norms %o% norms)
  rownames(sim) <- colnames(sim) <- rownames(m)
  sim
}

#' Redundant row-pairs above a cosine threshold
#' @keywords internal
.redundant_pairs <- function(sim, thresh) {
  pr <- which(sim > thresh & upper.tri(sim), arr.ind = TRUE)
  if (nrow(pr) == 0)
    return(data.frame(a = character(), b = character(), cosine = numeric()))
  data.frame(a = rownames(sim)[pr[, 1]], b = colnames(sim)[pr[, 2]],
             cosine = sim[pr], stringsAsFactors = FALSE)[
               order(-sim[pr]), ]
}


#' Text summary of LDA topics with SymPortal / host cross-reference
#'
#' Prints, for each topic: dominant clade, number of samples, mean entropy and
#' purity, top host species (if present), top SymPortal profiles (if present),
#' and the exemplar sequence composition.
#'
#' @param sp A `symbayes` object with LDA fitted.
#' @param top_n Number of top sequences per topic (default 8).
#' @return Invisibly `NULL`; called for its printed output.
#' @examples
#' \dontrun{ lda_topic_summary(sp) }
#' @export
lda_topic_summary <- function(sp, top_n = 8) {
  m <- .get_membership(sp, "lda")
  beta <- m$beta; K <- m$k
  md <- sp$metadata; sc <- sp$seq_clades

  cat(sprintf("LDA topic exemplars (K = %d)\n\n", K))
  for (k in seq_len(K)) {
    top <- utils::head(sort(beta[k, ], decreasing = TRUE), top_n)
    dom_clade <- sc[names(top)[1]]; if (is.na(dom_clade)) dom_clade <- "?"
    samples <- md[md$lda_dominant == k, , drop = FALSE]
    n <- nrow(samples)
    cat(sprintf("Topic %d [Clade %s] - %d samples (entropy %.2f, purity %.2f)\n",
                k, dom_clade, n,
                mean(samples$lda_entropy, na.rm = TRUE),
                mean(samples$lda_max_prop, na.rm = TRUE)))
    if ("host_species" %in% colnames(md) && n > 0) {
      ht <- table(samples$host_species)
      cat(sprintf("  hosts: %s\n",
                  paste(names(ht), ht, sep = ":", collapse = ", ")))
    }
    if ("sp_profile" %in% colnames(md) && n > 0) {
      pt <- sort(table(samples$sp_profile), decreasing = TRUE)
      cat(sprintf("  top profiles: %s\n",
                  paste(utils::head(names(pt), 3), utils::head(pt, 3),
                        sep = ":", collapse = ", ")))
    }
    for (i in seq_along(top))
      cat(sprintf("    %-15s [%s] %.1f%%\n", names(top)[i],
                  ifelse(is.na(sc[names(top)[i]]), "?", sc[names(top)[i]]),
                  100 * top[i]))
    cat("\n")
  }
  invisible(NULL)
}

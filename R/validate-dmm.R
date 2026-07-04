# ==============================================================================
# validate-dmm.R — DMM prediction, cross-validation, simulation
# ==============================================================================
# These operate on the DMM closed-form predictive and so are DMM-specific.
# ==============================================================================

#' Dirichlet-Multinomial log-likelihood of a count vector under one component
#' @keywords internal
.dmm_loglik <- function(counts, alpha) {
  N <- sum(counts); A <- sum(alpha)
  lgamma(A) - lgamma(N + A) + sum(lgamma(counts + alpha) - lgamma(alpha))
}

#' Predict DMM component membership for new samples
#'
#' Scores new count vectors against a fitted DMM's components using the
#' closed-form Dirichlet-Multinomial predictive, returning posterior membership.
#' Columns are aligned to the fitted sequences (missing zero-filled, extras
#' dropped). No refitting is required.
#'
#' @param sp A `symbayes` object with DMM fitted.
#' @param new_counts Matrix (samples x sequences) or a single named vector.
#' @return A data frame with `assigned`, `certainty`, `total_reads`, and per
#'   component posterior probabilities.
#' @examples
#' \dontrun{ predict_dmm(sp, new_counts) }
#' @export
predict_dmm <- function(sp, new_counts) {
  m <- .get_membership(sp, "dmm")
  fit <- m$meta$best_fit; K <- m$k
  alpha_mat <- DirichletMultinomial::fitted(fit)
  pi_k <- DirichletMultinomial::mixturewt(fit)$pi
  model_seqs <- rownames(alpha_mat)

  if (is.null(dim(new_counts)))
    new_counts <- matrix(new_counts, nrow = 1,
                         dimnames = list("new_1", names(new_counts)))

  aligned <- matrix(0L, nrow(new_counts), length(model_seqs),
                    dimnames = list(rownames(new_counts), model_seqs))
  shared <- intersect(colnames(new_counts), model_seqs)
  aligned[, shared] <- new_counts[, shared]

  out <- lapply(seq_len(nrow(aligned)), function(i) {
    x <- aligned[i, ]
    ll <- vapply(seq_len(K), function(k)
      .dmm_loglik(x, alpha_mat[, k]) + log(pi_k[k]), numeric(1))
    mx <- max(ll); post <- exp(ll - (mx + log(sum(exp(ll - mx)))))
    data.frame(sample = rownames(aligned)[i], assigned = which.max(post),
               certainty = max(post), total_reads = sum(x), t(post))
  })
  res <- do.call(rbind, out)
  colnames(res)[5:(4 + K)] <- paste0("P_comp_", seq_len(K))
  rownames(res) <- NULL
  res
}


#' Leave-one-out cross-validation for DMM
#'
#' For each sample: refit the DMM on the other N-1 samples and predict the
#' held-out sample. Reports concordance with the full-model assignment and
#' flags samples that switch component (candidate boundary/mixed samples).
#'
#' @param sp A `symbayes` object with DMM fitted.
#' @param seed Random seed (default 42).
#' @return A data frame of per-sample LOO results, returned invisibly.
#' @examples
#' \dontrun{ loo <- dmm_loo(sp) }
#' @export
dmm_loo <- function(sp, seed = 42) {
  m <- .get_membership(sp, "dmm"); K <- m$k
  set.seed(seed)
  mat <- sp$count_mat; N <- nrow(mat)

  rows <- lapply(seq_len(N), function(i) {
    fit <- tryCatch(DirichletMultinomial::dmn(mat[-i, , drop = FALSE], k = K,
                                              verbose = FALSE),
                    error = function(e) NULL)
    if (is.null(fit)) return(data.frame(sample_uid = rownames(mat)[i],
                                        loo_assigned = NA, loo_certainty = NA))
    am <- DirichletMultinomial::fitted(fit)
    pk <- DirichletMultinomial::mixturewt(fit)$pi
    x <- mat[i, ]
    ll <- vapply(seq_len(K), function(k)
      .dmm_loglik(x, am[, k]) + log(pk[k]), numeric(1))
    mx <- max(ll); post <- exp(ll - (mx + log(sum(exp(ll - mx)))))
    data.frame(sample_uid = rownames(mat)[i], loo_assigned = which.max(post),
               loo_certainty = max(post))
  })
  loo <- do.call(rbind, rows)
  full <- as.integer(sp$metadata$dmm_dominant[
    match(loo$sample_uid, sp$metadata$sample_uid)])
  loo$full_assigned <- full
  loo$concordant <- loo$loo_assigned == loo$full_assigned
  message(sprintf("LOO concordance: %.1f%% (%d/%d); %d samples switch",
                  100 * mean(loo$concordant, na.rm = TRUE),
                  sum(loo$concordant, na.rm = TRUE),
                  sum(!is.na(loo$concordant)),
                  sum(!loo$concordant, na.rm = TRUE)))
  invisible(loo)
}


#' Repeated subsampling cross-validation for DMM
#'
#' Holds out a fraction of samples, refits, predicts held-out; repeats. Reports
#' overall concordance and per-sample stability.
#'
#' @param sp A `symbayes` object with DMM fitted.
#' @param frac_test Fraction held out per replicate (default 0.2).
#' @param n_reps Replicates (default 50).
#' @param seed Random seed (default 42).
#' @return A data frame of per-sample stability, returned invisibly.
#' @examples
#' \dontrun{ dmm_subsample_cv(sp) }
#' @export
dmm_subsample_cv <- function(sp, frac_test = 0.2, n_reps = 50, seed = 42) {
  m <- .get_membership(sp, "dmm"); K <- m$k
  set.seed(seed)
  mat <- sp$count_mat; N <- nrow(mat)
  n_test <- max(1, round(N * frac_test))

  reps <- lapply(seq_len(n_reps), function(r) {
    idx <- sample(N, n_test)
    fit <- tryCatch(DirichletMultinomial::dmn(mat[-idx, , drop = FALSE], k = K,
                                              verbose = FALSE),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    am <- DirichletMultinomial::fitted(fit)
    pk <- DirichletMultinomial::mixturewt(fit)$pi
    do.call(rbind, lapply(idx, function(i) {
      x <- mat[i, ]
      ll <- vapply(seq_len(K), function(k)
        .dmm_loglik(x, am[, k]) + log(pk[k]), numeric(1))
      mx <- max(ll); post <- exp(ll - (mx + log(sum(exp(ll - mx)))))
      data.frame(sample_uid = rownames(mat)[i], predicted = which.max(post))
    }))
  })
  cv <- do.call(rbind, reps)
  full <- as.integer(sp$metadata$dmm_dominant[
    match(cv$sample_uid, sp$metadata$sample_uid)])
  cv$concordant <- cv$predicted == full
  stab <- dplyr::summarise(dplyr::group_by(cv, .data$sample_uid),
                           n = dplyr::n(),
                           pct_concordant = 100 * mean(.data$concordant),
                           .groups = "drop")
  message(sprintf("Subsample CV concordance: %.1f%%",
                  100 * mean(cv$concordant, na.rm = TRUE)))
  invisible(stab[order(stab$pct_concordant), ])
}


#' Simulate samples and test DMM assignment (blend / perturb / from_alpha)
#'
#' Three modes for probing DMM behaviour:
#' `from_alpha` draws from each component's Dirichlet then multinomial;
#' `blend` mixes two components across a gradient to trace decision boundaries;
#' `perturb` reassigns a fraction of reads in real samples to test robustness.
#'
#' @param sp A `symbayes` object with DMM fitted.
#' @param mode One of `"from_alpha"`, `"blend"`, `"perturb"`.
#' @param n_sim Number of simulated samples (default 100).
#' @param lib_size Reads per simulated sample (default 5000).
#' @param perturb_frac Fraction of reads to reassign for `"perturb"`
#'   (default 0.1).
#' @param blend_comps Two component indices to mix for `"blend"`
#'   (default `c(1, 2)`).
#' @param blend_range Mixing proportions for `"blend"`
#'   (default `seq(0, 1, 0.1)`).
#' @param seed Random seed (default 42).
#' @return A data frame of simulated samples with predicted assignments.
#' @examples
#' \dontrun{
#' dmm_simulate(sp, mode = "blend", blend_comps = c(1, 4))
#' dmm_simulate(sp, mode = "perturb", perturb_frac = 0.1)
#' }
#' @export
dmm_simulate <- function(sp, mode = c("from_alpha", "blend", "perturb"),
                         n_sim = 100, lib_size = 5000, perturb_frac = 0.1,
                         blend_comps = c(1, 2), blend_range = seq(0, 1, 0.1),
                         seed = 42) {
  mode <- match.arg(mode)
  m <- .get_membership(sp, "dmm"); K <- m$k
  set.seed(seed)
  alpha_mat <- DirichletMultinomial::fitted(m$meta$best_fit)
  seq_names <- rownames(alpha_mat)

  draw <- function(alpha) {
    props <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
    props <- props / sum(props)
    setNames(stats::rmultinom(1, lib_size, props)[, 1], seq_names)
  }

  if (mode == "from_alpha") {
    sims <- lapply(seq_len(n_sim), function(i) {
      k <- sample(K, 1)
      list(counts = draw(alpha_mat[, k]), true_comp = k,
           label = sprintf("sim_%03d_c%d", i, k))
    })
  } else if (mode == "blend") {
    c1 <- blend_comps[1]; c2 <- blend_comps[2]
    sims <- unlist(lapply(blend_range, function(w)
      lapply(seq_len(ceiling(n_sim / length(blend_range))), function(i)
        list(counts = draw(w * alpha_mat[, c1] + (1 - w) * alpha_mat[, c2]),
             true_comp = NA, blend_w = w,
             label = sprintf("blend_%.2f_%03d", w, i)))),
      recursive = FALSE)
  } else {
    mat <- sp$count_mat
    sims <- lapply(seq_len(n_sim), function(i) {
      idx <- sample(nrow(mat), 1); x <- mat[idx, ]
      n_move <- round(sum(x) * perturb_frac)
      rem <- which(x > 0)
      for (j in seq_len(n_move)) {
        pos <- sample(rem, 1, prob = x[rem]); x[pos] <- x[pos] - 1
        if (x[pos] == 0) rem <- rem[rem != pos]
      }
      for (pos in sample(length(x), n_move, TRUE)) x[pos] <- x[pos] + 1
      aligned <- setNames(rep(0L, length(seq_names)), seq_names)
      shared <- intersect(names(x), seq_names); aligned[shared] <- x[shared]
      list(counts = aligned,
           true_comp = as.integer(sp$metadata$dmm_dominant[
             match(rownames(mat)[idx], sp$metadata$sample_uid)]),
           label = sprintf("perturb_%s", rownames(mat)[idx]))
    })
  }

  sim_mat <- do.call(rbind, lapply(sims, `[[`, "counts"))
  rownames(sim_mat) <- vapply(sims, `[[`, character(1), "label")
  pred <- predict_dmm(sp, sim_mat)
  pred$true_comp <- vapply(sims, function(s)
    if (is.null(s$true_comp)) NA_integer_ else s$true_comp, integer(1))
  if (mode == "blend") pred$blend_w <- vapply(sims, `[[`, numeric(1), "blend_w")

  if (mode %in% c("from_alpha", "perturb") && any(!is.na(pred$true_comp)))
    message(sprintf("%s accuracy: %.1f%%", mode,
                    100 * mean(pred$assigned == pred$true_comp, na.rm = TRUE)))
  pred
}

# ==============================================================================
# fit-hdp.R — Hierarchical Dirichlet Process (threshold-free, mixture-aware)
# ==============================================================================

#' Fit a Hierarchical Dirichlet Process
#'
#' The threshold-free, mixture-aware core of `symbayes`. Unlike DMM and LDA,
#' the number of community types is **inferred from the data** rather than
#' specified: the DP concentration parameters (with Gamma hyperpriors) govern
#' component number, and samples may be admixed. Uses the \pkg{hdp} package
#' (Roberts). See <https://github.com/nicolaroberts/hdp>.
#'
#' With `group_var = NULL` the structure is flat (one parent DP, one child DP
#' per sample) and assignment is covariate-free. Supplying a metadata column
#' (e.g. `"host_species"`) inserts a middle DP layer so samples within a group
#' share strength preferentially. **Caution:** grouping by a covariate you
#' later want to *test* (e.g. host effect on symbiont community) introduces
#' circularity; use flat structure for host-blind claims, grouped only as an
#' explicitly-labelled predictive enhancement.
#'
#' Writes a membership object to `sp$hdp`. `theta` excludes the residual
#' (unassigned) mass, which is kept in `meta$residual`. `meta` also holds
#' `entropy`, `max_prop`, `n_comp` per sample, `raw_ncl` (sampler cluster
#' count), and the `multi` chain object.
#'
#' @param sp A filtered `symbayes` object.
#' @param group_var Optional metadata column defining a middle DP layer, or
#'   `NULL` (default) for a flat, covariate-free structure.
#' @param n_chains Independent Gibbs chains (default 4).
#' @param burnin Burn-in iterations per chain (default 5000).
#' @param n Posterior samples collected per chain (default 100).
#' @param space Iterations between collected samples (default 50).
#' @param initcc Initial cluster count (default 10; a starting point only).
#' @param alphaa,alphab Gamma hyperprior (shape, rate) on DP concentration
#'   (default 1, 1).
#' @param cos_merge Cosine similarity above which raw clusters are merged during
#'   extraction (default 0.90; the \pkg{hdp} default). Higher retains more
#'   components. Do not set very low or everything collapses.
#' @param min_sample Minimum posterior samples a component must appear in to be
#'   retained (default 1).
#' @param seed Base random seed (chains use `seed + chain`).
#' @return The `symbayes` object with `sp$hdp` populated.
#' @examples
#' \dontrun{
#' sp <- fit_hdp(sp)                       # flat, threshold-free
#' hdp_plot_numcomp(sp)
#' hdp_vs_symportal(sp)
#' }
#' @export
fit_hdp <- function(sp, group_var = NULL,
                    n_chains = 4, burnin = 5000, n = 100, space = 50,
                    initcc = 10, alphaa = 1, alphab = 1,
                    cos_merge = 0.90, min_sample = 1, seed = 42) {
  .check_sp(sp)
  if (!requireNamespace("hdp", quietly = TRUE))
    stop("Install hdp: remotes::install_github('nicolaroberts/hdp')",
         call. = FALSE)
  if (!isTRUE(sp$filtered))
    warning("Fitting on unfiltered data; consider filter_samples() first.",
            call. = FALSE)

  mat <- sp$count_mat
  storage.mode(mat) <- "integer"
  n_samp  <- nrow(mat)
  n_categ <- ncol(mat)
  message(sprintf("Fitting HDP (%d samples x %d sequences, %d chains) ...",
                  n_samp, n_categ, n_chains))

  built <- .hdp_build(sp, mat, group_var, alphaa, alphab, n_categ, n_samp)
  hdp_obj <- hdp::dp_activate(built$obj, built$activate, initcc, seed = seed)

  chains <- lapply(seq_len(n_chains), function(ch) {
    hdp::hdp_posterior(hdp_obj, burnin = burnin, n = n, space = space,
                       cpiter = 3, seed = seed + ch)
  })

  multi <- hdp::hdp_multi_chain(chains)
  multi <- hdp::hdp_extract_components(multi, cos.merge = cos_merge,
                                       min.sample = min_sample)

  n_comp <- length(hdp::comp_categ_counts(multi)) - 1   # exclude residual
  raw_ncl <- tryCatch(median(unlist(lapply(hdp::chains(multi), hdp::numcluster))),
                      error = function(e) NA_real_)
  message(sprintf("  HDP: %d components inferred (raw sampler median %.0f)",
                  n_comp, raw_ncl))
  if (!is.na(raw_ncl) && raw_ncl > 2 * n_comp && n_comp < raw_ncl - 2)
    message("  note: large raw-vs-extracted gap; if under-resolved try cos_merge ~0.95")

  # --- Extract per-sample and per-component distributions ---------------------
  # comp_dp_distn()$mean can be a single matrix OR a list of matrices depending
  # on the hdp version / object structure. Coerce to a matrix.
  dp_mean <- hdp::comp_dp_distn(multi)$mean
  if (is.list(dp_mean) && !is.data.frame(dp_mean)) {
    # list of per-DP row vectors -> bind into a matrix
    dp_distn <- do.call(rbind, dp_mean)
  } else {
    dp_distn <- as.matrix(dp_mean)
  }
  sample_rows <- (nrow(dp_distn) - n_samp + 1):nrow(dp_distn)
  theta_full <- dp_distn[sample_rows, , drop = FALSE]
  storage.mode(theta_full) <- "double"
  rownames(theta_full) <- rownames(mat)

  residual <- as.numeric(theta_full[, 1])
  names(residual) <- rownames(theta_full)
  theta <- theta_full[, -1, drop = FALSE]
  storage.mode(theta) <- "double"
  colnames(theta) <- paste0("HDP_", seq_len(ncol(theta)))

  cd_mean <- hdp::comp_categ_distn(multi)$mean
  cd_mean <- if (is.list(cd_mean) && !is.data.frame(cd_mean))
    do.call(rbind, cd_mean) else as.matrix(cd_mean)
  beta <- cd_mean[-1, , drop = FALSE]
  storage.mode(beta) <- "double"
  rownames(beta) <- paste0("HDP_", seq_len(nrow(beta)))
  colnames(beta) <- colnames(mat)

  dominant <- setNames(as.integer(apply(theta, 1, which.max)), rownames(theta))
  entropy  <- as.numeric(apply(theta, 1, .entropy))
  max_prop <- as.numeric(apply(theta, 1, max))
  n_comp_s <- as.integer(apply(theta, 1, function(p) sum(p > 0.05)))
  names(entropy) <- names(max_prop) <- names(n_comp_s) <- rownames(theta)

  md <- sp$metadata
  md$hdp_dominant <- factor(dominant[md$sample_uid])
  md$hdp_entropy  <- entropy[md$sample_uid]
  md$hdp_max_prop <- max_prop[md$sample_uid]
  md$hdp_n_comp   <- n_comp_s[md$sample_uid]
  md$hdp_residual <- residual[md$sample_uid]
  sp$metadata <- .map_symportal(sp, md, "hdp")

  sp$hdp <- .membership(
    theta = theta, beta = beta,
    dominant = factor(dominant), k = n_comp,
    model = "hdp", label_prefix = "HDP",
    meta = list(entropy = entropy, max_prop = max_prop, n_comp = n_comp_s,
                residual = residual, raw_ncl = raw_ncl,
                group_var = group_var, multi = multi)
  )
  sp
}


#' Build the HDP tree structure (flat or grouped)
#' @keywords internal
.hdp_build <- function(sp, mat, group_var, alphaa, alphab, n_categ, n_samp) {
  if (is.null(group_var)) {
    ppindex <- c(0, rep(1, n_samp))
    cpindex <- c(1, rep(2, n_samp))
    obj <- hdp::hdp_init(ppindex = ppindex, cpindex = cpindex,
                         hh = rep(1, n_categ),
                         alphaa = rep(alphaa, 2), alphab = rep(alphab, 2))
    obj <- hdp::hdp_setdata(obj, 2:(n_samp + 1), mat)
    list(obj = obj, activate = 1:(n_samp + 1))
  } else {
    if (!group_var %in% colnames(sp$metadata))
      stop(sprintf("group_var '%s' not in metadata", group_var), call. = FALSE)
    grp <- droplevels(factor(sp$metadata[rownames(mat), group_var]))
    n_grp <- nlevels(grp)
    sample_parent <- as.integer(grp) + 1
    ppindex <- c(0, rep(1, n_grp), sample_parent)
    cpindex <- c(1, rep(2, n_grp), rep(3, n_samp))
    obj <- hdp::hdp_init(ppindex = ppindex, cpindex = cpindex,
                         hh = rep(1, n_categ),
                         alphaa = rep(alphaa, 3), alphab = rep(alphab, 3))
    sample_dp <- (n_grp + 2):(n_grp + 1 + n_samp)
    obj <- hdp::hdp_setdata(obj, sample_dp, mat)
    list(obj = obj, activate = 1:(n_grp + 1 + n_samp))
  }
}

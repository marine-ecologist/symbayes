# ==============================================================================
# fit-dmm.R — Dirichlet-Multinomial Mixture
# ==============================================================================

#' Fit a Dirichlet-Multinomial Mixture model
#'
#' Fits DMMs for `k = 1:max_k`, selects the number of components by the Laplace
#' approximation to the model evidence (with BIC reported alongside), extracts
#' posterior membership, and, if SymPortal profiles are present, maps each
#' sample's dominant profile and clade.
#'
#' DMM is a *hard-assignment* model: each sample is assigned to one component.
#' Its posterior certainty is often near 1.0 even for mixed samples, because it
#' models between-sample variation, not within-sample mixing. Use [fit_lda] or
#' [fit_hdp] to represent admixture. See [dmm_loo] and [dmm_simulate] to probe
#' assignment stability.
#'
#' Writes a membership object to `sp$dmm` with the standard schema
#' (`theta`, `beta`, `dominant`, `k`, `meta`). The `meta` list contains
#' `certainty` (per-sample max posterior), `model_fit` (the K-selection table),
#' and `best_fit` (the fitted `DMN` object, used by [dmm_predict] etc.).
#'
#' @param sp A filtered `symbayes` object.
#' @param max_k Maximum number of components to evaluate (default 10).
#' @param seed Random seed (default 42).
#' @return The `symbayes` object with `sp$dmm` populated.
#' @examples
#' \dontrun{
#' sp <- fit_dmm(sp, max_k = 10)
#' sp$dmm$k
#' dmm_plot_selection(sp)
#' }
#' @export
fit_dmm <- function(sp, max_k = 10, seed = 42) {
  .check_sp(sp)
  if (!isTRUE(sp$filtered))
    warning("Fitting on unfiltered data; consider filter_samples() first.",
            call. = FALSE)
  set.seed(seed)

  mat <- sp$count_mat
  message(sprintf("Fitting DMM for k = 1:%d ...", max_k))

  fits <- lapply(seq_len(max_k), function(k) {
    DirichletMultinomial::dmn(mat, k = k, verbose = FALSE)
  })

  laplace <- vapply(fits, DirichletMultinomial::laplace, numeric(1))
  aic_val <- vapply(fits, DirichletMultinomial::AIC, numeric(1))
  bic_val <- vapply(fits, DirichletMultinomial::BIC, numeric(1))
  model_fit <- data.frame(K = seq_len(max_k), Laplace = laplace,
                          AIC = aic_val, BIC = bic_val)

  best_K  <- which.min(laplace)
  best    <- fits[[best_K]]

  # Posterior membership (samples x components)
  theta <- DirichletMultinomial::mixture(best, assign = FALSE)
  dominant  <- setNames(apply(theta, 1, which.max), rownames(theta))
  certainty <- setNames(apply(theta, 1, max), rownames(theta))

  # Exemplars: fitted Dirichlet parameters, normalised to compositions
  fp <- DirichletMultinomial::fitted(best)     # sequences x components
  beta <- t(sweep(fp, 2, colSums(fp), "/"))    # components x sequences

  # Attach dominant + certainty to metadata for convenience
  md <- sp$metadata
  md$dmm_dominant  <- factor(dominant[md$sample_uid])
  md$dmm_certainty <- certainty[md$sample_uid]
  sp$metadata <- .map_symportal(sp, md, "dmm")

  sp$dmm <- .membership(
    theta = theta, beta = beta,
    dominant = factor(dominant), k = best_K,
    model = "dmm", label_prefix = "Comp",
    meta = list(certainty = certainty, model_fit = model_fit,
                best_fit = best, fits = fits)
  )

  message(sprintf("  DMM: best k = %d (Laplace); mean certainty = %.3f",
                  best_K, mean(certainty)))
  sp
}


#' Map SymPortal dominant profile and clade onto metadata
#' @keywords internal
.map_symportal <- function(sp, md, model) {
  if (is.null(sp$prof_mat) || is.null(sp$prof_meta)) return(md)

  prof_aligned <- sp$prof_mat[rownames(sp$count_mat), , drop = FALSE]
  dom_idx <- apply(prof_aligned, 1, which.max)
  dom_uid <- setNames(colnames(prof_aligned)[dom_idx], rownames(prof_aligned))

  pm <- sp$prof_meta
  name_map  <- setNames(pm$`ITS2 type profile`,
                        as.character(pm$`ITS2 type profile UID`))
  clade_map <- setNames(pm$Clade,
                        as.character(pm$`ITS2 type profile UID`))

  if (is.null(md$sp_profile)) md$sp_profile <- name_map[dom_uid[md$sample_uid]]
  if (is.null(md$sp_clade))   md$sp_clade   <- clade_map[dom_uid[md$sample_uid]]
  md
}

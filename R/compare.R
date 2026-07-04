# ==============================================================================
# compare.R — cross-model comparison
# ==============================================================================

#' Compare group counts and selection strategy across fitted models
#'
#' Tabulates, for each fitted model (and SymPortal if present), the number of
#' groups, the selection strategy, and whether the model represents
#' within-sample mixing.
#'
#' @param sp A `symbayes` object with one or more models fitted.
#' @return A data frame, printed and returned invisibly.
#' @examples
#' \dontrun{ compare_models(sp) }
#' @export
compare_models <- function(sp) {
  .check_sp(sp)
  rows <- list()
  if (!is.null(sp$prof_mat)) {
    pa <- sp$prof_mat[rownames(sp$count_mat), , drop = FALSE]
    rows[["SymPortal"]] <- data.frame(
      method = "SymPortal", n_groups = sum(colSums(pa) > 0),
      selection = "deterministic co-occurrence",
      mixture_aware = "partial (multi-profile)")
  }
  if (!is.null(sp$dmm)) rows[["dmm"]] <- data.frame(
    method = "DMM", n_groups = sp$dmm$k,
    selection = "Laplace/BIC (K specified)", mixture_aware = "no")
  if (!is.null(sp$lda)) rows[["lda"]] <- data.frame(
    method = "LDA", n_groups = sp$lda$k,
    selection = "cosine redundancy (K specified)", mixture_aware = "yes")
  if (!is.null(sp$hdp)) rows[["hdp"]] <- data.frame(
    method = "HDP", n_groups = sp$hdp$k,
    selection = "non-parametric (K inferred)", mixture_aware = "yes")

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  print(out)
  invisible(out)
}


#' Compare within-sample mixing across models and SymPortal
#'
#' Computes per-sample "mixedness" (1 - dominant proportion) and the effective
#' number of groups (Hill number, q = 1) for LDA, HDP, DMM, and SymPortal
#' (multi-profile proportions), and reports correlations between them. DMM
#' mixedness is ~0 by construction (hard assignment).
#'
#' @param sp A `symbayes` object.
#' @return A data frame of per-sample mixing metrics, returned invisibly.
#' @examples
#' \dontrun{ mix <- compare_mixing(sp) }
#' @export
compare_mixing <- function(sp) {
  .check_sp(sp)
  mix <- data.frame(sample_uid = sp$metadata$sample_uid,
                    stringsAsFactors = FALSE)
  hill <- function(p) { p <- p[p > 0]; if (length(p) == 0) return(NA_real_)
  exp(-sum(p * log(p))) }

  # Safe row lookup: returns model rows aligned to mix$sample_uid, with NA for
  # any sample absent from that model's theta (metadata and a model's theta can
  # diverge if the object was subset/filtered after fitting).
  align_rows <- function(theta, uids) {
    idx <- match(uids, rownames(theta))
    theta[idx, , drop = FALSE]   # NA rows where idx is NA
  }

  for (model in c("lda", "hdp", "dmm")) {
    if (is.null(sp[[model]])) next
    th <- align_rows(sp[[model]]$theta, mix$sample_uid)
    mix[[paste0(model, "_mix")]]  <- 1 - apply(th, 1, function(p)
      if (all(is.na(p))) NA_real_ else max(p, na.rm = TRUE))
    mix[[paste0(model, "_hill")]] <- apply(th, 1, function(p)
      if (all(is.na(p))) NA_real_ else hill(p))
  }

  if (!is.null(sp$prof_mat)) {
    pr_full <- sp$prof_mat
    pr_full <- pr_full / pmax(rowSums(pr_full), 1)
    pr <- align_rows(pr_full, mix$sample_uid)
    mix$sp_mix  <- apply(pr, 1, function(p)
      if (all(is.na(p))) NA_real_ else 1 - max(p, na.rm = TRUE))
    mix$sp_hill <- apply(pr, 1, function(p)
      if (all(is.na(p))) NA_real_ else hill(p))
  }

  cat("Mixing concordance\n")
  if (all(c("lda_mix", "sp_mix") %in% names(mix)))
    cat(sprintf("  LDA vs SymPortal mixedness: r = %.3f\n",
                stats::cor(mix$lda_mix, mix$sp_mix, use = "complete.obs")))
  if (all(c("hdp_mix", "sp_mix") %in% names(mix)))
    cat(sprintf("  HDP vs SymPortal mixedness: r = %.3f\n",
                stats::cor(mix$hdp_mix, mix$sp_mix, use = "complete.obs")))
  if ("dmm_mix" %in% names(mix))
    cat(sprintf("  DMM mean mixedness: %.4f (~0 = hard assignment)\n",
                mean(mix$dmm_mix, na.rm = TRUE)))
  invisible(mix)
}


#' Compare a soft (fractional) model against SymPortal's profiles via the
#' fuzzy Rand index
#'
#' Quantifies how much a fractional model (LDA or HDP) agrees with SymPortal's
#' ITS2 type-profile assignments at the level of *sample pairs*, using the
#' fuzzy / probabilistic Rand index (Campello 2007; Huellermeier et al. 2011).
#' This is a principled upgrade over a scalar mixedness correlation: it compares
#' the full relational structure (which samples cluster with which, fractionally)
#' rather than a per-sample summary.
#'
#' For each pair of samples the probability they are co-assigned is the dot
#' product of their membership vectors (`theta %*% t(theta)`). SymPortal's
#' multi-profile proportions are treated the same way, giving two soft
#' co-membership matrices. The fuzzy Rand index is `1 - mean|EQ_model - EQ_sp|`
#' over all pairs; the adjusted version corrects for chance agreement (the raw
#' index is inflated by the many pairs that are clearly separate under both).
#'
#' Also returns per-sample co-membership divergence (where the two methods
#' disagree most — typically the admixed samples) and a disagreement breakdown
#' separating "SymPortal missed mixing the model found" from "SymPortal invented
#' a compound profile the model sees as pure".
#'
#' @param sp A `symbayes` object with the chosen model fitted and SymPortal
#'   profiles present.
#' @param model One of `"lda"` or `"hdp"` (fractional models). DMM is hard-
#'   assignment and not meaningful here.
#' @param mixed_thresh A sample is "mixed" if its dominant proportion is below
#'   this (default 0.7). Applied to BOTH the model's `theta` and SymPortal's
#'   multi-profile proportions, so "SymPortal mixed" means SymPortal assigned
#'   the sample multiple profiles with non-trivial mass -- not whether the
#'   profile name contains a separator (SymPortal names list intragenomic DIVs
#'   and almost always contain '-' or '/').
#' @return Invisibly, a list with `fuzzy_rand`, `adjusted_fuzzy_rand`,
#'   `sample_divergence` (data frame, most-divergent first, with `model_mixed`
#'   and `sp_mixed` flags), and `disagreement` (counts of the two disagreement
#'   types: SymPortal missed mixing vs SymPortal over-split).
#' @examples
#' \dontrun{
#' res <- compare_soft_hard(sp, model = "lda")
#' res$fuzzy_rand
#' head(res$sample_divergence)
#' }
#' @export
compare_soft_hard <- function(sp, model = c("lda", "hdp"),
                              mixed_thresh = 0.7) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  if (is.null(sp$prof_mat))
    stop("No SymPortal profiles in this object.", call. = FALSE)

  # --- Align samples present in both model theta and SymPortal profiles -------
  theta <- m$theta
  prof  <- sp$prof_mat

  # Guard against duplicate rownames, which would inflate the co-membership
  # matrices (a common symptom of metadata/matrix desync).
  if (anyDuplicated(rownames(theta)))
    theta <- theta[!duplicated(rownames(theta)), , drop = FALSE]
  if (anyDuplicated(rownames(prof)))
    prof  <- prof[!duplicated(rownames(prof)), , drop = FALSE]

  common <- intersect(rownames(theta), rownames(prof))
  common <- unique(common)
  if (length(common) < 3)
    stop("Fewer than 3 samples shared between model and SymPortal.",
         call. = FALSE)
  theta <- theta[common, , drop = FALSE]
  prof  <- prof[common,  , drop = FALSE]

  # Drop SymPortal profiles absent from these samples
  prof <- prof[, colSums(prof) > 0, drop = FALSE]

  # Row-normalise both to proper membership vectors (soft co-membership)
  theta <- theta / pmax(rowSums(theta), .Machine$double.eps)
  prof  <- prof  / pmax(rowSums(prof),  .Machine$double.eps)

  # --- Soft co-membership matrices: EQ_ij = probability i,j co-assigned ------
  eq_model <- theta %*% t(theta)
  eq_sp    <- prof  %*% t(prof)

  # --- Fuzzy Rand index = 1 - mean pairwise disagreement (upper triangle) ----
  ut <- upper.tri(eq_model)
  disagree <- abs(eq_model[ut] - eq_sp[ut])
  fuzzy_rand <- 1 - mean(disagree)

  # --- Chance-adjusted fuzzy Rand --------------------------------------------
  # Expected disagreement if the two co-membership patterns were independent,
  # using the observed marginal co-membership levels. ARI-style correction:
  # (index - expected) / (max - expected), max = 1 (perfect agreement).
  em <- eq_model[ut]; es <- eq_sp[ut]
  expected_agree <- mean(em) * mean(es) + (1 - mean(em)) * (1 - mean(es))
  adjusted_fuzzy_rand <- (fuzzy_rand - expected_agree) / (1 - expected_agree)

  # --- Per-sample co-membership divergence -----------------------------------
  # Mean absolute difference in how each sample co-clusters with all others.
  n <- length(common)
  row_div <- (rowSums(abs(eq_model - eq_sp)) - 0) / (n - 1)  # exclude self ~0

  md <- sp$metadata
  dom_prop <- apply(theta, 1, max)
  sp_dom_prop <- apply(prof, 1, max)   # prof is row-normalised profile props

  # "SymPortal mixed" = SymPortal assigned this sample MULTIPLE profiles with
  # non-trivial proportions (its dominant profile holds < mixed_thresh of the
  # sample's profile mass). This is the correct analogue of model mixedness --
  # NOT whether the dominant profile's NAME contains a separator (SymPortal
  # profile names list intragenomic DIVs and almost always contain '-' or '/',
  # so a name-based test flags nearly everything).
  sp_mixed_vec <- sp_dom_prop < mixed_thresh

  sample_div <- data.frame(
    sample_uid    = common,
    sample_name   = if ("sample_name" %in% colnames(md))
      md$sample_name[match(common, md$sample_uid)] else common,
    divergence    = round(row_div, 4),
    model_domprop = round(dom_prop, 3),
    model_mixed   = dom_prop < mixed_thresh,
    sp_domprop    = round(sp_dom_prop, 3),
    sp_mixed      = sp_mixed_vec,
    sp_profile    = if ("sp_profile" %in% colnames(md))
      md$sp_profile[match(common, md$sample_uid)] else NA_character_,
    stringsAsFactors = FALSE)
  sample_div <- sample_div[order(-sample_div$divergence), ]

  # --- Disagreement breakdown -------------------------------------------------
  model_mixed <- sample_div$model_mixed
  sp_mixed    <- sample_div$sp_mixed

  # Two asymmetric disagreement types, both defined by MULTI-PROFILE structure
  # (not profile-name shape):
  #  (a) model mixed, SymPortal single-profile -> SymPortal MISSED mixing
  #  (b) model pure,  SymPortal multi-profile   -> SymPortal SPLIT a pure sample
  missed_mixing  <- sum(model_mixed & !sp_mixed, na.rm = TRUE)
  false_split    <- sum(!model_mixed & sp_mixed, na.rm = TRUE)
  both_mixed     <- sum(model_mixed & sp_mixed, na.rm = TRUE)
  both_pure      <- sum(!model_mixed & !sp_mixed, na.rm = TRUE)

  disagreement <- data.frame(
    category = c("both mixed (agree)",
                 "both single (agree)",
                 "model mixed / SymPortal single (SymPortal missed mixing)",
                 "model single / SymPortal multi (SymPortal over-split)"),
    n = c(both_mixed, both_pure, missed_mixing, false_split))

  cat(sprintf("Soft-vs-hard comparison: %s vs SymPortal (%d samples)\n",
              toupper(model), length(common)))
  cat(sprintf("  Fuzzy Rand index:          %.3f\n", fuzzy_rand))
  cat(sprintf("  Adjusted fuzzy Rand:       %.3f  (chance-corrected)\n",
              adjusted_fuzzy_rand))
  cat(sprintf("  both mixed / both single (agree):  %d / %d\n",
              both_mixed, both_pure))
  cat(sprintf("  %s mixed / SymPortal single (SymPortal missed mixing): %d\n",
              toupper(model), missed_mixing))
  cat(sprintf("  %s single / SymPortal multi (SymPortal over-split):    %d\n",
              toupper(model), false_split))

  invisible(list(
    fuzzy_rand          = fuzzy_rand,
    adjusted_fuzzy_rand = adjusted_fuzzy_rand,
    sample_divergence   = sample_div,
    disagreement        = disagreement,
    eq_model            = eq_model,
    eq_sp               = eq_sp))
}


#' Visualise soft-vs-hard co-membership agreement as heatmaps
#'
#' Renders the pairwise co-membership structure behind [compare_soft_hard]:
#' the model's soft co-membership matrix, SymPortal's, and their difference.
#' Samples are ordered by the model's dominant group so block structure is
#' visible. In the difference panel, warm cells are pairs the model
#' co-clusters more than SymPortal (mixing SymPortal splits apart); cool cells
#' are the reverse (SymPortal groups pairs the model separates).
#'
#' @param sp A `symbayes` object with the chosen model fitted and SymPortal
#'   profiles present.
#' @param model One of `"lda"` or `"hdp"`.
#' @param panels Which panels to show: any of `"model"`, `"sp"`, `"diff"`
#'   (default all three).
#' @param mixed_thresh Passed to [compare_soft_hard] (default 0.7).
#' @return A \pkg{patchwork} object (or a single ggplot if one panel).
#' @examples
#' \dontrun{ plot_soft_hard(sp, model = "lda") }
#' @export
plot_soft_hard <- function(sp, model = c("lda", "hdp"),
                           panels = c("model", "sp", "diff"),
                           mixed_thresh = 0.7) {
  model <- match.arg(model)
  panels <- match.arg(panels, c("model", "sp", "diff"), several.ok = TRUE)

  cmp <- compare_soft_hard(sp, model = model, mixed_thresh = mixed_thresh)
  eq_model <- cmp$eq_model
  eq_sp    <- cmp$eq_sp
  common   <- rownames(eq_model)

  # Order samples by the model's dominant group, then dominant proportion
  m <- .get_membership(sp, model)
  th <- m$theta[common, , drop = FALSE]
  dom  <- apply(th, 1, which.max)
  domp <- apply(th, 1, max)
  ord  <- common[order(dom, -domp)]

  eq_model <- eq_model[ord, ord]
  eq_sp    <- eq_sp[ord, ord]
  eq_diff  <- eq_model - eq_sp

  # Long-form helper for a co-membership matrix
  melt_eq <- function(mat) {
    df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
    names(df) <- c("i", "j", "value")
    df$i <- factor(df$i, levels = ord)
    df$j <- factor(df$j, levels = rev(ord))
    df
  }

  heat <- function(mat, title, diverging = FALSE) {
    df <- melt_eq(mat)
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$i, .data$j, fill = .data$value)) +
      ggplot2::geom_raster() +
      ggplot2::labs(title = title, x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     legend.key.size = ggplot2::unit(0.35, "cm")) +
      ggplot2::coord_equal()
    if (diverging) {
      lim <- max(abs(df$value), na.rm = TRUE)
      p + ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
        limits = c(-lim, lim), name = "model - SP")
    } else {
      p + ggplot2::scale_fill_viridis_c(option = "magma", limits = c(0, 1),
                                        name = "P(co-assign)")
    }
  }

  plist <- list()
  if ("model" %in% panels)
    plist[["model"]] <- heat(eq_model, sprintf("%s co-membership", toupper(model)))
  if ("sp" %in% panels)
    plist[["sp"]] <- heat(eq_sp, "SymPortal co-membership")
  if ("diff" %in% panels)
    plist[["diff"]] <- heat(eq_diff, sprintf("Difference (%s - SymPortal)",
                                             toupper(model)), diverging = TRUE)

  if (length(plist) == 1) return(plist[[1]])
  patchwork::wrap_plots(plist, nrow = 1) +
    patchwork::plot_annotation(
      title = sprintf("Soft-vs-hard co-membership: %s vs SymPortal",
                      toupper(model)),
      subtitle = sprintf(
        "Fuzzy Rand = %.3f (adjusted %.3f); samples ordered by %s dominant group",
        cmp$fuzzy_rand, cmp$adjusted_fuzzy_rand, toupper(model)),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")))
}

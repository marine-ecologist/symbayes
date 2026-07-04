# ==============================================================================
# plot-hdp.R — HDP-specific diagnostics and over-splitting quantification
# ==============================================================================

#' Plot the HDP posterior distribution of component number
#'
#' Shows two quantities honestly: the *raw* cluster count the Gibbs sampler
#' used at each posterior sample (a distribution), and the number of *extracted*
#' stable components after merging cosine-similar and dropping transient
#' clusters (a single line). The raw-to-extracted reduction is expected, not a
#' loss.
#'
#' @param sp A `symbayes` object with HDP fitted.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ hdp_plot_numcomp(sp) }
#' @export
hdp_plot_numcomp <- function(sp) {
  m <- .get_membership(sp, "hdp")
  multi <- m$meta$multi
  chs <- hdp::chains(multi)
  ncl <- do.call(rbind, lapply(seq_along(chs), function(i)
    data.frame(chain = i, ncluster = hdp::numcluster(chs[[i]]))))

  extracted <- m$k
  raw_med   <- median(ncl$ncluster)
  x_min <- min(extracted, min(ncl$ncluster)) - 1
  x_max <- max(ncl$ncluster) + 1

  ggplot2::ggplot(ncl, ggplot2::aes(.data$ncluster)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "#2c7fb8", colour = "white",
                            alpha = 0.85) +
    ggplot2::geom_vline(xintercept = extracted, linetype = "dashed",
                        colour = "#d62728", linewidth = 0.8) +
    ggplot2::annotate("text", x = extracted, y = Inf,
                      label = sprintf("extracted = %d", extracted),
                      colour = "#d62728", hjust = 1.1, vjust = 1.5, size = 3.5) +
    ggplot2::annotate("text", x = raw_med, y = Inf,
                      label = sprintf("raw median = %.0f", raw_med),
                      colour = "#2c7fb8", hjust = -0.1, vjust = 1.5, size = 3.5) +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max)) +
    ggplot2::labs(
      title = "HDP component number: inferred, not fixed",
      subtitle = "Blue: raw clusters per posterior sample. Red: stable extracted components.",
      x = "Number of clusters / components", y = "Posterior samples") +
    ggplot2::theme_minimal(base_size = 12)
}


#' Quantify SymPortal over-classification relative to HDP
#'
#' The core method metric. Compares SymPortal's profile count to HDP's inferred
#' component count and, crucially, partitions SymPortal profiles into three
#' mutually exclusive classes so they do not contaminate one another:
#'
#' \describe{
#'   \item{Over-splitting}{Several profiles that *confidently* map to one HDP
#'     component (`dominant_frac >= over_split_min_frac`). SymPortal split one
#'     real community type into multiple named profiles.}
#'   \item{Mixtures}{Profiles whose samples spread across components
#'     (`dominant_frac < mixture_thresh`). SymPortal called an admixture a
#'     unique type; these are not force-assigned.}
#'   \item{Orphans}{Not confidently in any component and not clear mixtures
#'     (usually rare, n = 1). Reported separately and excluded from the
#'     over-splitting count.}
#' }
#'
#' The partition matters: a rare profile whose single sample has its largest
#' component at only a few percent must not be counted as "over-split into"
#' that component.
#'
#' @param sp A `symbayes` object with HDP fitted and SymPortal profiles present.
#' @param over_split_min_frac Minimum dominant-component share for a profile to
#'   count as over-split (default 0.7).
#' @param mixture_thresh Dominant-component share below which a profile is a
#'   mixture (default 0.7).
#' @return Invisibly, a list with `n_sp_profiles`, `n_hdp_comp`,
#'   `oversplit_ratio`, `mapping`, `redundant_splits`, `mixture_profiles`
#'   (profile-level), `mixture_samples` (one row per sample whose dominant
#'   profile is a mixture, with ranked `comp_i`/`frac_i` columns showing which
#'   HDP components co-occur in that individual sample), and `orphan_profiles`.
#' @examples
#' \dontrun{ res <- hdp_vs_symportal(sp) }
#' @export
hdp_vs_symportal <- function(sp, over_split_min_frac = 0.7,
                             mixture_thresh = 0.7) {
  m <- .get_membership(sp, "hdp")
  if (is.null(sp$prof_mat)) stop("No SymPortal profiles in this object.",
                                 call. = FALSE)
  theta <- m$theta
  md <- sp$metadata

  prof_aligned <- sp$prof_mat[rownames(sp$count_mat), , drop = FALSE]
  prof_present <- colnames(prof_aligned)[colSums(prof_aligned) > 0]
  n_sp <- length(prof_present)
  n_hdp <- m$k

  cat(sprintf("SymPortal over-classification vs HDP\n"))
  cat(sprintf("  SymPortal profiles: %d\n  HDP components: %d\n  Ratio: %.1fx\n",
              n_sp, n_hdp, n_sp / n_hdp))

  dom_idx <- apply(prof_aligned, 1, which.max)
  dom_prof <- setNames(colnames(prof_aligned)[dom_idx], rownames(prof_aligned))
  name_map <- setNames(sp$prof_meta$`ITS2 type profile`,
                       as.character(sp$prof_meta$`ITS2 type profile UID`))

  rows <- lapply(prof_present, function(puid) {
    samps <- names(dom_prof)[dom_prof == puid]
    if (length(samps) == 0) return(NULL)
    mt <- colMeans(theta[samps, , drop = FALSE])
    data.frame(sp_profile = name_map[puid], n_samples = length(samps),
               dominant_hdp = which.max(mt), dominant_frac = max(mt),
               hdp_entropy = .entropy(mt),
               is_mixture = max(mt) < mixture_thresh,
               stringsAsFactors = FALSE)
  })
  mapping <- do.call(rbind, rows)
  mapping$confident_assign <- mapping$dominant_frac >= over_split_min_frac

  confident <- mapping[mapping$confident_assign, ]
  h2p <- dplyr::summarise(dplyr::group_by(confident, .data$dominant_hdp),
                          n_sp_profiles = dplyr::n(),
                          profiles = paste(.data$sp_profile, collapse = ", "),
                          .groups = "drop")
  redundant <- h2p[h2p$n_sp_profiles > 1, ]
  mixtures  <- mapping[mapping$is_mixture, ]
  orphans   <- mapping[!mapping$confident_assign & !mapping$is_mixture, ]

  cat(sprintf("\n  Over-splitting (confident): %d HDP components absorb >1 profile\n",
              nrow(redundant)))
  if (nrow(redundant) > 0)
    for (i in seq_len(nrow(redundant)))
      cat(sprintf("    HDP_%d <- %d profiles: %s\n",
                  redundant$dominant_hdp[i], redundant$n_sp_profiles[i],
                  redundant$profiles[i]))
  cat(sprintf("\n  Mixture-as-taxa profiles (<%.0f%% one comp): %d\n",
              100 * mixture_thresh, nrow(mixtures)))
  cat(sprintf("  Orphan profiles (excluded from over-split count): %d\n",
              nrow(orphans)))

  # --- Per-sample mixture breakdown -------------------------------------------
  # For every sample whose dominant SymPortal profile is a flagged mixture,
  # report which HDP components co-occur in that individual sample, as ranked
  # columns (comp_1 / frac_1 = largest, comp_2 / frac_2 = second, ...).
  mixture_samples <- .mixture_sample_table(
    sp, theta, dom_prof, name_map,
    mixture_uids = names(dom_prof)[
      name_map[dom_prof] %in% mixtures$sp_profile],
    top_n = 3, min_frac = 0.05)

  invisible(list(n_sp_profiles = n_sp, n_hdp_comp = n_hdp,
                 oversplit_ratio = n_sp / n_hdp, mapping = mapping,
                 redundant_splits = redundant, mixture_profiles = mixtures,
                 mixture_samples = mixture_samples,
                 orphan_profiles = orphans))
}


#' Build a per-sample HDP mixture breakdown table
#'
#' One row per sample, listing the top contributing HDP components and their
#' proportions as ranked columns.
#'
#' @param sp A `symbayes` object.
#' @param theta HDP sample x component matrix.
#' @param dom_prof Named vector: sample_uid -> dominant profile UID.
#' @param name_map Named vector: profile UID -> profile name.
#' @param mixture_uids Sample UIDs whose dominant profile is a flagged mixture.
#' @param top_n Number of ranked component columns to emit (default 3).
#' @param min_frac Minimum proportion to include a component (default 0.05).
#' @return A data frame: sample_uid, sample_name, sp_profile, n_comp, then
#'   comp_1/frac_1 ... comp_top_n/frac_top_n.
#' @keywords internal
.mixture_sample_table <- function(sp, theta, dom_prof, name_map,
                                  mixture_uids, top_n = 3, min_frac = 0.05) {
  if (length(mixture_uids) == 0) return(NULL)
  md <- sp$metadata
  have_name <- "sample_name" %in% colnames(md)

  rows <- lapply(mixture_uids, function(uid) {
    if (!uid %in% rownames(theta)) return(NULL)
    v <- sort(theta[uid, ], decreasing = TRUE)
    v <- v[v >= min_frac]
    if (length(v) == 0) v <- sort(theta[uid, ], decreasing = TRUE)[1]

    base <- data.frame(
      sample_uid  = uid,
      sample_name = if (have_name)
        md$sample_name[match(uid, md$sample_uid)] else NA_character_,
      sp_profile  = unname(name_map[dom_prof[uid]]),
      n_comp      = length(v),
      stringsAsFactors = FALSE)

    # Ranked component / fraction columns
    for (i in seq_len(top_n)) {
      base[[paste0("comp_", i)]] <- if (i <= length(v)) names(v)[i] else NA_character_
      base[[paste0("frac_", i)]] <- if (i <= length(v))
        round(unname(v[i]), 4) else NA_real_
    }
    base
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(-out$n_comp, out$sp_profile), ]
}

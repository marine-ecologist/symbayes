# ==============================================================================
# plot-generalist.R — visualisations that apply to any fitted model
# ==============================================================================
# These read the standardised membership object and dispatch on `model=`.
# ==============================================================================

#' PCoA ordination coloured by a fitted model's grouping (or any variable)
#'
#' Principal coordinates analysis of Bray-Curtis distances among samples,
#' coloured by the dominant group of a fitted model, or by any metadata column.
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`. Colours points by that
#'   model's dominant group. Ignored if `colour_by` is set.
#' @param colour_by Optional metadata column name to colour by instead
#'   (e.g. `"host_species"`, `"sp_clade"`).
#' @param colours Optional named colour vector.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{
#' plot_pcoa(sp, model = "hdp")
#' plot_pcoa(sp, colour_by = "sp_clade")
#' }
#' @export
plot_pcoa <- function(sp, model = "dmm", colour_by = NULL, colours = NULL) {
  .check_sp(sp)
  bc   <- vegan::vegdist(sp$count_mat, method = "bray")
  pc   <- stats::cmdscale(bc, k = 2, eig = TRUE)
  ve   <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)

  df <- data.frame(sample_uid = rownames(pc$points),
                   PC1 = pc$points[, 1], PC2 = pc$points[, 2])
  df <- dplyr::left_join(df, sp$metadata, by = "sample_uid")

  if (is.null(colour_by)) {
    .get_membership(sp, model)               # validate model is fitted
    colour_by <- paste0(model, "_dominant")
    title <- sprintf("PCoA coloured by %s dominant group", toupper(model))
  } else {
    title <- sprintf("PCoA coloured by %s", colour_by)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$PC1, .data$PC2,
                                        colour = .data[[colour_by]])) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::labs(x = sprintf("PCoA 1 (%.1f%%)", ve[1]),
                  y = sprintf("PCoA 2 (%.1f%%)", ve[2]),
                  title = title, colour = colour_by) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::coord_equal()
  if (!is.null(colours))
    p <- p + ggplot2::scale_colour_manual(values = colours, na.value = "grey60")
  p
}


#' Stacked sequence-composition barplot faceted by a model's groups
#'
#' ITS2 sequence composition per sample, faceted by the dominant group of a
#' fitted model, using SymPortal colours where available. Optionally annotates
#' each bar with the model's assignment confidence (DMM certainty / LDA-HDP
#' dominant proportion).
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param min_pct Minimum percent in any sample to show a sequence (default 1).
#' @param sample_label Metadata column for x labels (default `"sample_name"`).
#' @param show_conf Annotate bars with confidence (default `TRUE`).
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_barplot(sp, model = "hdp") }
#' @export
plot_barplot <- function(sp, model = "dmm", min_pct = 1,
                         sample_label = "sample_name", show_conf = TRUE) {
  .check_sp(sp); m <- .get_membership(sp, model)
  md <- sp$metadata
  if (!sample_label %in% colnames(md)) sample_label <- "sample_uid"

  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")

  rel <- sp$count_mat / rowSums(sp$count_mat)
  show_seqs <- names(which(apply(rel, 2, max) * 100 >= min_pct))
  rel_show  <- rel[, show_seqs, drop = FALSE]

  bar_df <- as.data.frame(rel_show, check.names = FALSE)
  bar_df$sample_uid <- rownames(rel_show)
  bar_df$Other <- 1 - rowSums(rel_show)
  bar_df <- dplyr::left_join(
    bar_df, md[, c("sample_uid", sample_label, dom_col, conf_col)],
    by = "sample_uid")
  bar_df <- tidyr::pivot_longer(bar_df,
                                cols = c(dplyr::all_of(show_seqs), "Other"),
                                names_to = "sequence", values_to = "rel_abund")

  ord <- md[order(md[[dom_col]], -md[[conf_col]]), ]
  bar_df[[sample_label]] <- factor(bar_df[[sample_label]],
                                   levels = ord[[sample_label]])
  seq_order <- c(names(sort(colSums(rel_show), decreasing = TRUE)), "Other")
  bar_df$sequence <- factor(bar_df$sequence, levels = seq_order)

  seq_cols <- .build_colours(show_seqs, sp$col_dict)
  n_per <- table(md[[dom_col]])

  p <- ggplot2::ggplot(bar_df, ggplot2::aes(.data[[sample_label]],
                                            .data$rel_abund,
                                            fill = .data$sequence)) +
    ggplot2::geom_bar(stat = "identity", width = 1, linewidth = 0.1,
                      colour = "white") +
    ggplot2::scale_fill_manual(values = seq_cols,
                               guide = ggplot2::guide_legend(ncol = 4)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, if (show_conf) 0.12 else 0.02))) +
    ggplot2::facet_wrap(ggplot2::vars(.data[[dom_col]]),
                        nrow = 1, scales = "free_x") +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::labs(title = sprintf("ITS2 composition by %s group (k = %d)",
                                  toupper(model), m$k),
                  x = NULL, y = "Relative abundance", fill = "Sequence") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                          size = 5),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = ggplot2::unit(0.35, "cm"))

  if (show_conf) {
    cf <- md[, c(sample_label, dom_col, conf_col)]
    cf[[sample_label]] <- factor(cf[[sample_label]],
                                 levels = levels(bar_df[[sample_label]]))
    cf$label <- sprintf("%.2f", cf[[conf_col]])
    p <- p + ggplot2::geom_text(
      data = cf, ggplot2::aes(x = .data[[sample_label]], y = 1.02,
                              label = .data$label),
      inherit.aes = FALSE, size = 1.8, angle = 90, hjust = 0, colour = "grey30")
  }
  p
}


#' Exemplar composition per group (what a "pure" sample looks like)
#'
#' Stacked bar of each group's exemplar sequence distribution (`beta`), with a
#' clade annotation and dominant-sample count per group.
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param min_pct Minimum percent in any group to show a sequence (default 1).
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_exemplars(sp, model = "lda") }
#' @export
plot_exemplars <- function(sp, model = "dmm", min_pct = 1) {
  .check_sp(sp); m <- .get_membership(sp, model)
  beta <- m$beta; K <- m$k; sc <- sp$seq_clades
  dom_col <- paste0(model, "_dominant")

  show_seqs <- names(which(apply(beta, 2, max) * 100 >= min_pct))
  beta_show <- beta[, show_seqs, drop = FALSE]
  other <- 1 - rowSums(beta_show)

  n_per <- table(factor(sp$metadata[[dom_col]], levels = seq_len(K)))

  ex <- as.data.frame(beta_show, check.names = FALSE)
  ex$group <- rownames(beta_show); ex$Other <- other
  ex <- tidyr::pivot_longer(ex, -"group", names_to = "sequence",
                            values_to = "proportion")
  ex <- ex[ex$proportion > 0, ]

  seq_order <- c(names(sort(colSums(beta_show), decreasing = TRUE)), "Other")
  ex$sequence <- factor(ex$sequence, levels = rev(seq_order))
  ex$sequence <- droplevels(ex$sequence)

  labs_v <- vapply(seq_len(K), function(k)
    sprintf("%s\n(n=%d)", rownames(beta)[k],
            ifelse(is.na(n_per[as.character(k)]), 0, n_per[as.character(k)])),
    character(1))
  names(labs_v) <- rownames(beta)
  ex$group_label <- factor(labs_v[ex$group], levels = labs_v)

  topclade <- sc[apply(beta, 1, function(r) names(which.max(r)))]
  seq_cols <- .build_colours(show_seqs, sp$col_dict)

  p <- ggplot2::ggplot(ex, ggplot2::aes(.data$group_label, .data$proportion,
                                        fill = .data$sequence)) +
    ggplot2::geom_bar(stat = "identity", width = 0.8, linewidth = 0.2,
                      colour = "white") +
    ggplot2::scale_fill_manual(values = seq_cols, drop = TRUE,
                               guide = ggplot2::guide_legend(ncol = 4)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05)) +
    ggplot2::labs(title = sprintf("%s exemplars (k = %d)", toupper(model), K),
                  subtitle = "Each bar = group's archetypal composition; n = dominant samples",
                  x = NULL, y = "P(sequence | group)", fill = "ITS2 sequence") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = ggplot2::unit(0.3, "cm"))

  clade_annot <- data.frame(group_label = factor(labs_v, levels = labs_v),
                            clade = topclade, y = 1.02)
  p + ggplot2::geom_text(data = clade_annot,
                         ggplot2::aes(.data$group_label, .data$y,
                                      label = .data$clade),
                         inherit.aes = FALSE, size = 3.5, fontface = "bold",
                         colour = "grey30")
}


#' Membership heatmap (per-sample group proportions)
#'
#' Heatmap of `theta`: near-binary for DMM, gradient for LDA/HDP.
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_membership(sp, model = "lda") }
#' @export
plot_membership <- function(sp, model = "dmm") {
  .check_sp(sp); m <- .get_membership(sp, model)
  theta <- m$theta
  dom_col <- paste0(model, "_dominant")
  ord <- sp$metadata$sample_uid[order(sp$metadata[[dom_col]])]

  df <- as.data.frame(theta); df$sample_uid <- rownames(theta)
  df <- tidyr::pivot_longer(df, -"sample_uid", names_to = "group",
                            values_to = "prop")
  df$sample_uid <- factor(df$sample_uid, levels = ord)
  df$group <- factor(df$group, levels = rev(colnames(theta)))

  ggplot2::ggplot(df, ggplot2::aes(.data$sample_uid, .data$group,
                                   fill = .data$prop)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    ggplot2::labs(title = sprintf("%s membership (k = %d)", toupper(model), m$k),
                  x = "Sample", y = NULL, fill = "Proportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}


#' Contingency heatmap: a model's groups vs SymPortal profiles
#'
#' Tile grid of SymPortal ITS2 type profiles (rows, grouped by clade) against a
#' model's groups (columns); tile colour = mean assignment confidence, number =
#' sample count.
#'
#' @param sp A `symbayes` object with SymPortal profiles present.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_contingency(sp, model = "hdp") }
#' @export
plot_contingency <- function(sp, model = "dmm") {
  .check_sp(sp); .get_membership(sp, model)
  md <- sp$metadata
  if (is.null(md$sp_profile))
    stop("No SymPortal profiles in this object.", call. = FALSE)

  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")

  cont <- md[!is.na(md$sp_profile), ]
  cont <- dplyr::summarise(
    dplyr::group_by(cont, .data$sp_profile, .data$sp_clade,
                    grp = .data[[dom_col]]),
    n = dplyr::n(), conf = mean(.data[[conf_col]], na.rm = TRUE),
    .groups = "drop")

  ord <- dplyr::arrange(
    dplyr::count(md[!is.na(md$sp_profile), ], .data$sp_profile, .data$sp_clade),
    .data$sp_clade, dplyr::desc(.data$n))
  cont$sp_profile <- factor(cont$sp_profile, levels = ord$sp_profile)

  ggplot2::ggplot(cont, ggplot2::aes(.data$grp, .data$sp_profile,
                                     fill = .data$conf)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n), size = 3.5,
                       fontface = "bold") +
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = c(0, 1),
                                  name = "Mean\nconfidence") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$sp_clade),
                        scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::labs(title = sprintf("SymPortal profile vs %s group", toupper(model)),
                  x = sprintf("%s group", toupper(model)),
                  y = "ITS2 type profile") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(size = 7),
      panel.grid = ggplot2::element_blank())
}


#' Top characterising sequences per group
#'
#' Faceted bar chart of the highest-weight sequences in each group's exemplar,
#' coloured by clade.
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param top_n Sequences per group (default 8).
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_top_seqs(sp, model = "dmm") }
#' @export
plot_top_seqs <- function(sp, model = "dmm", top_n = 8) {
  .check_sp(sp); m <- .get_membership(sp, model)
  beta <- m$beta; K <- m$k; sc <- sp$seq_clades

  rows <- lapply(seq_len(K), function(k) {
    top <- utils::head(sort(beta[k, ], decreasing = TRUE), top_n)
    data.frame(group = rownames(beta)[k], sequence = names(top),
               weight = as.numeric(top), clade = sc[names(top)])
  })
  df <- do.call(rbind, rows)
  df$sequence <- stats::reorder(df$sequence, df$weight)

  ggplot2::ggplot(df, ggplot2::aes(.data$weight, .data$sequence,
                                   fill = .data$clade)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$group), scales = "free_y") +
    ggplot2::scale_fill_manual(values = .clade_palette(), na.value = "grey60") +
    ggplot2::labs(title = sprintf("Top sequences per %s group", toupper(model)),
                  x = "P(sequence | group)", y = NULL, fill = "Clade") +
    ggplot2::theme_minimal(base_size = 11)
}

#' Admixture plot: mixed samples, sequence mix over boxed component structure
#'
#' Filters to samples that are genuinely *mixed* under a fitted model (at least
#' `min_groups` groups above `min_frac`) and returns a two-panel figure sharing
#' one x-axis:
#'
#' \describe{
#'   \item{Top}{ITS2 **sequence** composition of the admixed samples (a single
#'     stacked barplot, no facet) -- the raw mixture SymPortal sees.}
#'   \item{Bottom}{**Component** proportions (`theta`) for the same samples,
#'     with every contributing segment outlined so the within-sample admixture
#'     is demarcated.}
#' }
#'
#' Together they show what SymPortal would collapse to a single compound profile
#' name (top) versus the fractional community structure the mixture model
#' recovers (bottom). DMM is near-degenerate here (hard assignment), so this is
#' most informative for `"lda"` and `"hdp"`.
#'
#' @param sp A `symbayes` object.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param min_frac Minimum group proportion to count as present (default 0.05).
#' @param min_groups Minimum number of present groups for a sample to be shown
#'   (default 2; i.e. only mixed samples).
#' @param seq_min_pct Minimum percent in any admixed sample to show a sequence
#'   individually in the top panel (default 1).
#' @param sample_label Metadata column for x labels (default `"sample_name"`).
#' @param box_colour Outline colour for each within-sample component segment
#'   (default `"grey15"`).
#' @param label_profile If `TRUE` (default) and SymPortal profiles are present,
#'   annotate the bottom panel with each sample's dominant SymPortal profile.
#' @param combine If `TRUE` (default) return a \pkg{patchwork} of both panels;
#'   if `FALSE` return a list with `seq` and `admix` ggplots.
#' @return A \pkg{patchwork} object (or a list if `combine = FALSE`).
#' @examples
#' \dontrun{ plot_admixtures(sp, model = "lda") }
#' @export
plot_admixtures <- function(sp, model = "hdp", min_frac = 0.05,
                            min_groups = 2, seq_min_pct = 1,
                            sample_label = "sample_name",
                            box_colour = "grey15", label_profile = TRUE,
                            combine = TRUE) {
  .check_sp(sp); m <- .get_membership(sp, model)
  md <- sp$metadata
  if (!sample_label %in% colnames(md)) sample_label <- "sample_uid"

  theta <- m$theta
  n_present <- apply(theta, 1, function(p) sum(p >= min_frac))
  mixed <- names(n_present)[n_present >= min_groups]
  if (length(mixed) == 0)
    stop(sprintf("No samples with >= %d groups above %.2f under model '%s'.",
                 min_groups, min_frac, model), call. = FALSE)

  # --- Shared sample ordering: dominant group, then dominant proportion -------
  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")
  sub <- md[match(mixed, md$sample_uid), , drop = FALSE]
  ord_uid <- sub$sample_uid[order(sub[[dom_col]], -sub[[conf_col]])]
  ord_lab <- md[[sample_label]][match(ord_uid, md$sample_uid)]

  labeller_fac <- function(uids) {
    factor(md[[sample_label]][match(uids, md$sample_uid)], levels = ord_lab)
  }

  # === TOP PANEL: ITS2 sequence composition (no facet) ========================
  rel <- sp$count_mat[mixed, , drop = FALSE]
  rel <- rel / rowSums(rel)
  show_seqs <- names(which(apply(rel, 2, max) * 100 >= seq_min_pct))
  rel_show <- rel[, show_seqs, drop = FALSE]

  sdf <- as.data.frame(rel_show, check.names = FALSE)
  sdf$sample_uid <- rownames(rel_show)
  sdf$Other <- 1 - rowSums(rel_show)
  sdf <- tidyr::pivot_longer(sdf, cols = c(dplyr::all_of(show_seqs), "Other"),
                             names_to = "sequence", values_to = "rel_abund")
  sdf$.x <- labeller_fac(sdf$sample_uid)
  seq_order <- c(names(sort(colSums(rel_show), decreasing = TRUE)), "Other")
  sdf$sequence <- factor(sdf$sequence, levels = seq_order)
  seq_cols <- .build_colours(show_seqs, sp$col_dict)

  p_seq <- ggplot2::ggplot(sdf, ggplot2::aes(.data$.x, .data$rel_abund,
                                             fill = .data$sequence)) +
    ggplot2::geom_col(width = 1, linewidth = 0.1, colour = "white") +
    ggplot2::scale_fill_manual(values = seq_cols, name = "ITS2 sequence",
                               guide = ggplot2::guide_legend(ncol = 3)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(
      title = sprintf("Admixed samples under %s (%d samples)",
                      toupper(model), length(mixed)),
      subtitle = "Top: ITS2 sequence mix (what SymPortal profiles). Bottom: recovered component structure.",
      x = NULL, y = "Seq.\nabundance") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = ggplot2::unit(0.3, "cm"))

  # === BOTTOM PANEL: boxed component proportions ==============================
  th <- theta[mixed, , drop = FALSE]
  long <- as.data.frame(th, check.names = FALSE)
  long$sample_uid <- rownames(th)
  long <- tidyr::pivot_longer(long, -"sample_uid", names_to = "group",
                              values_to = "prop")
  long <- long[long$prop >= min_frac, ]
  long$.x <- labeller_fac(long$sample_uid)
  grp_levels <- colnames(theta)
  long$group <- factor(long$group, levels = grp_levels)
  grp_cols <- setNames(scales::hue_pal()(length(grp_levels)), grp_levels)

  p_admix <- ggplot2::ggplot(long, ggplot2::aes(.data$.x, .data$prop,
                                                fill = .data$group)) +
    ggplot2::geom_col(width = 0.9, colour = box_colour, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = grp_cols, name = toupper(model)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(subtitle = "Each boxed segment = one component's share within that sample",
                  x = NULL, y = "Component\nproportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                          size = 6),
      panel.grid.major.x = ggplot2::element_blank())

  if (label_profile && "sp_profile" %in% colnames(md)) {
    ann <- md[md$sample_uid %in% ord_uid, c("sample_uid", "sp_profile")]
    ann$.x <- labeller_fac(ann$sample_uid)
    p_admix <- p_admix + ggplot2::geom_text(
      data = ann, ggplot2::aes(x = .data$.x, y = 1.02, label = .data$sp_profile),
      inherit.aes = FALSE, angle = 90, hjust = 0, size = 1.9, colour = "grey30")
  }

  if (!combine) return(list(seq = p_seq, admix = p_admix))
  patchwork::wrap_plots(list(p_seq, p_admix), ncol = 1, heights = c(1, 1.4))
}

#' Four-row SymPortal-vs-model comparison over all samples
#'
#' Stacks four aligned panels sharing one sample order (all samples, ordered by
#' the chosen model's dominant group then confidence):
#'
#' \enumerate{
#'   \item **Model assignment** -- a single colour band per sample showing the
#'     model's dominant group.
#'   \item **Model admixture** -- boxed component proportions (`theta`); each
#'     outlined segment is one component's share within the sample.
#'   \item **SymPortal assignment** -- a single colour band per sample showing
#'     its dominant SymPortal profile.
#'   \item **SymPortal mix** -- stacked SymPortal multi-profile proportions;
#'     samples with more than one profile show the admixture SymPortal records.
#' }
#'
#' Rows 1-2 are the model's view (hard then soft); rows 3-4 are SymPortal's
#' (hard then soft). Reading a column top-to-bottom shows, for one sample, how
#' the model and SymPortal each assign and each mix it.
#'
#' @param sp A `symbayes` object with the chosen model fitted and SymPortal
#'   profiles present.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param sample_label Metadata column for x labels (default `"sample_name"`).
#' @param box_colour Outline colour for admixture segments (default `"grey15"`).
#' @param show_x Draw sample labels on the bottom row (default `FALSE`; with
#'   many samples they are unreadable).
#' @return A \pkg{patchwork} object (four rows).
#' @examples
#' \dontrun{ plot_symportal_comparison(sp, model = "hdp") }
#' @export
plot_symportal_comparison <- function(sp, model = "hdp",
                                      sample_label = "sample_name",
                                      box_colour = "grey15", show_x = FALSE) {
  .check_sp(sp); m <- .get_membership(sp, model)
  if (is.null(sp$prof_mat))
    stop("No SymPortal profiles in this object.", call. = FALSE)
  md <- sp$metadata
  if (!sample_label %in% colnames(md)) sample_label <- "sample_uid"

  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")

  # --- Shared sample order: all samples, by model dominant group then conf ----
  theta <- m$theta
  common <- rownames(theta)
  sub <- md[match(common, md$sample_uid), , drop = FALSE]
  ord_uid <- sub$sample_uid[order(sub[[dom_col]], -sub[[conf_col]])]
  ord_lab <- md[[sample_label]][match(ord_uid, md$sample_uid)]
  xf <- function(uids) factor(md[[sample_label]][match(uids, md$sample_uid)],
                              levels = ord_lab)

  # --- Group boundaries (for boxes spanning all rows) -------------------------
  # Runs of consecutive samples sharing the model's dominant group, in x-order.
  ord_dom <- as.integer(sub[[dom_col]])[match(ord_uid, sub$sample_uid)]
  rle_dom <- rle(ord_dom)
  ends   <- cumsum(rle_dom$lengths)
  starts <- c(1, utils::head(ends, -1) + 1)
  # Rect x-bounds in discrete-axis units (bars centred on integers 1..n)
  box_df <- data.frame(xmin = starts - 0.5, xmax = ends + 0.5)

  # A rectangle layer drawn identically in every panel so the boxes line up
  # into continuous columns across the stacked figure.
  box_layer <- function() {
    ggplot2::geom_rect(
      data = box_df,
      mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                             ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = NA, colour = "grey20", linewidth = 0.4)
  }

  # Common x-axis theme; only the bottom row shows labels
  x_theme <- function(bottom = FALSE) {
    if (bottom && show_x)
      ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5, size = 5))
    else
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
  }

  # === ROW 1: model dominant-group assignment band ============================
  r1 <- data.frame(sample_uid = common,
                   group = factor(paste0(m$label_prefix, "_",
                                         as.integer(sub[[dom_col]])),
                                  levels = colnames(theta)))
  r1$.x <- xf(r1$sample_uid)
  grp_cols <- setNames(scales::hue_pal()(ncol(theta)), colnames(theta))

  p1 <- ggplot2::ggplot(r1, ggplot2::aes(.data$.x, 1, fill = .data$group)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = grp_cols, name = toupper(model)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = sprintf("%s assignment", toupper(model)),
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) + x_theme()

  # === ROW 2: model admixture (boxed component proportions) ===================
  a <- as.data.frame(theta, check.names = FALSE)
  a$sample_uid <- rownames(theta)
  a <- tidyr::pivot_longer(a, -"sample_uid", names_to = "group",
                           values_to = "prop")
  a <- a[a$prop > 0, ]
  a$.x <- xf(a$sample_uid)
  a$group <- factor(a$group, levels = colnames(theta))

  p2 <- ggplot2::ggplot(a, ggplot2::aes(.data$.x, .data$prop,
                                        fill = .data$group)) +
    ggplot2::geom_col(width = 0.9, colour = box_colour, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = grp_cols, name = toupper(model)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = sprintf("%s admixture", toupper(model)),
                  x = NULL, y = "Proportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) + x_theme()

  # === SymPortal proportions (shared by rows 3 and 4) =========================
  prof <- sp$prof_mat[common, , drop = FALSE]
  prof <- prof[, colSums(prof) > 0, drop = FALSE]
  prof <- prof / pmax(rowSums(prof), .Machine$double.eps)

  # Map profile UID -> name for legible labels
  pm <- sp$prof_meta
  pname <- if (!is.null(pm))
    setNames(pm$`ITS2 type profile`, as.character(pm$`ITS2 type profile UID`))
  else setNames(colnames(prof), colnames(prof))
  prof_names <- ifelse(colnames(prof) %in% names(pname),
                       pname[colnames(prof)], colnames(prof))

  # === ROW 3: SymPortal dominant-profile assignment band ======================
  sp_dom <- colnames(prof)[apply(prof, 1, which.max)]
  r3 <- data.frame(sample_uid = common,
                   profile = factor(pname[sp_dom], levels = sort(unique(prof_names))))
  r3$.x <- xf(r3$sample_uid)
  prof_cols <- setNames(scales::hue_pal()(length(unique(prof_names))),
                        sort(unique(prof_names)))

  p3 <- ggplot2::ggplot(r3, ggplot2::aes(.data$.x, 1, fill = .data$profile)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = prof_cols, name = "SymPortal profile",
                               guide = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "SymPortal assignment", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) + x_theme()

  # === ROW 4: SymPortal multi-profile mix =====================================
  sdf <- as.data.frame(prof, check.names = FALSE)
  sdf$sample_uid <- rownames(prof)
  sdf <- tidyr::pivot_longer(sdf, -"sample_uid", names_to = "profile_uid",
                             values_to = "prop")
  sdf <- sdf[sdf$prop > 0, ]
  sdf$profile <- factor(pname[sdf$profile_uid],
                        levels = sort(unique(prof_names)))
  sdf$.x <- xf(sdf$sample_uid)

  p4 <- ggplot2::ggplot(sdf, ggplot2::aes(.data$.x, .data$prop,
                                          fill = .data$profile)) +
    ggplot2::geom_col(width = 0.9, colour = box_colour, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = prof_cols, name = "SymPortal profile",
                               guide = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "SymPortal mix", x = NULL, y = "Proportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) +
    x_theme(bottom = TRUE)

  # Overlay the group boxes on every panel so they line up into columns
  p1 <- p1 + box_layer()
  p2 <- p2 + box_layer()
  p3 <- p3 + box_layer()
  p4 <- p4 + box_layer()

  patchwork::wrap_plots(list(p1, p2, p3, p4), ncol = 1,
                        heights = c(0.4, 1, 0.4, 1)) +
    patchwork::plot_annotation(
      title = sprintf("%s vs SymPortal: assignment and admixture (%d samples)",
                      toupper(model), length(common)),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")))
}

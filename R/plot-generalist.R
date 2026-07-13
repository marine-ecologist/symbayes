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
#' @param scale `"relative"` (default; each sample sums to 1) or `"absolute"`
#'   (raw read counts; bars reach different heights, revealing depth). In
#'   absolute mode the confidence annotation is disabled.
#' @param panel_width `"proportional"` (default; panel width scales with the
#'   number of samples in each group, via `facet_grid` + `space = "free_x"`) or
#'   `"equal"` (all group panels the same width, via `facet_wrap`).
#' @examples
#' \dontrun{ plot_barplot(sp, model = "hdp") }
#' @export
plot_barplot <- function(sp, model = "dmm", min_pct = 1,
                         sample_label = "sample_name", show_conf = TRUE,
                         scale = c("relative", "absolute"),
                         panel_width = c("proportional", "equal")) {
  scale <- match.arg(scale)
  panel_width <- match.arg(panel_width)
  if (scale == "absolute") show_conf <- FALSE
  .check_sp(sp); m <- .get_membership(sp, model)
  md <- sp$metadata
  if (!sample_label %in% colnames(md)) sample_label <- "sample_uid"

  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")

  rel <- sp$count_mat / rowSums(sp$count_mat)
  # Sequence selection is always on relative abundance (which seqs to show).
  show_seqs <- names(which(apply(rel, 2, max) * 100 >= min_pct))

  # Plotted values: relative (sum to 1) or absolute (raw counts).
  if (scale == "relative") {
    mat_show <- rel[, show_seqs, drop = FALSE]
    other    <- 1 - rowSums(mat_show)
    y_lab    <- "Relative abundance"
  } else {
    mat_show <- sp$count_mat[, show_seqs, drop = FALSE]
    other    <- rowSums(sp$count_mat) - rowSums(mat_show)
    y_lab    <- "Reads"
  }
  rel_show <- mat_show

  bar_df <- as.data.frame(rel_show, check.names = FALSE)
  bar_df$sample_uid <- rownames(rel_show)
  bar_df$Other <- other
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
    (if (panel_width == "proportional")
      ggplot2::facet_grid(cols = ggplot2::vars(.data[[dom_col]]),
                          scales = "free_x", space = "free_x")
     else
       ggplot2::facet_wrap(ggplot2::vars(.data[[dom_col]]),
                           nrow = 1, scales = "free_x")) +
    ggplot2::scale_x_discrete(drop = TRUE) +
    ggplot2::labs(title = sprintf("ITS2 composition by %s group (k = %d)",
                                  toupper(model), m$k),
                  x = NULL, y = y_lab, fill = "Sequence") +
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
#' @param scale Sequence panel (top) scaling: `"relative"` (default, sums to 1)
#'   or `"absolute"` (raw reads, revealing depth). The component panel (bottom)
#'   is always proportional (`theta`).
#' @return A \pkg{patchwork} object (or a list if `combine = FALSE`).
#' @examples
#' \dontrun{ plot_admixtures(sp, model = "lda") }
#' @export
plot_admixtures <- function(sp, model = "hdp", min_frac = 0.05,
                            min_groups = 2, seq_min_pct = 1,
                            sample_label = "sample_name",
                            box_colour = "grey15", label_profile = TRUE,
                            combine = TRUE,
                            scale = c("relative", "absolute")) {
  scale <- match.arg(scale)
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
  cm <- sp$count_mat[mixed, , drop = FALSE]
  rel <- cm / rowSums(cm)
  show_seqs <- names(which(apply(rel, 2, max) * 100 >= seq_min_pct))

  if (scale == "relative") {
    seq_mat <- rel[, show_seqs, drop = FALSE]
    seq_other <- 1 - rowSums(seq_mat)
    seq_ylab <- "Seq.\nabundance"
  } else {
    seq_mat <- cm[, show_seqs, drop = FALSE]
    seq_other <- rowSums(cm) - rowSums(seq_mat)
    seq_ylab <- "Reads"
  }
  rel_show <- seq_mat

  sdf <- as.data.frame(rel_show, check.names = FALSE)
  sdf$sample_uid <- rownames(rel_show)
  sdf$Other <- seq_other
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
      x = NULL, y = seq_ylab) +
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
  # Clade-grouped, deterministic palette (stable across plots)
  tmap <- tryCatch({
    nmn <- name_topics(sp, model = model)
    setNames(nmn$name, nmn$group)
  }, error = function(e) setNames(grp_levels, grp_levels))
  grp_cols <- .topic_palette(grp_levels, names_map = tmap)

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
#' @param screen If `TRUE`, row 1 uses *screen-aware* assignment
#'   ([screen_intragenomic]): intragenomic topic pairs are collapsed to one unit
#'   before labelling, and samples that genuinely mix distinct units (a
#'   "mixing"-verdict pair) are shown as a distinct "true mix" band rather than
#'   forced to a single plurality topic. Default `FALSE` (raw `which.max`).
#' @param majority Minimum merged-unit proportion for a sample to be called that
#'   unit under `screen = TRUE`; below it the sample is "true mix" (default 0.6).
#' @param max_legend Maximum number of SymPortal profiles to list in the legend
#'   (the most frequent dominant profiles). All profiles are still coloured;
#'   only the legend is capped, to stay on-panel for large datasets (default 15).
#' @param residual How to show model mass assigned to no component in the
#'   admixture row: `"renormalise"` (default; rescale each bar to sum to 1) or
#'   `"show"` (grey Residual segment so bars reach 1.0).
#' @param use_names If `TRUE`, label topics/components in the model legend with
#'   their SymPortal-style composition names (from [name_topics]) instead of
#'   `Topic_N` / `HDP_N` (default `FALSE`).
#' @return A \pkg{patchwork} object (four rows).
#' @examples
#' \dontrun{ plot_symportal_comparison(sp, model = "hdp") }
#' @export
plot_symportal_comparison <- function(sp, model = "hdp",
                                      sample_label = "sample_name",
                                      box_colour = "grey15", show_x = FALSE,
                                      screen = FALSE, majority = 0.6,
                                      max_legend = 15, n_seqs = 20,
                                      residual = c("renormalise", "show"),
                                      use_names = FALSE) {
  residual <- match.arg(residual)
  .check_sp(sp); m <- .get_membership(sp, model)
  if (is.null(sp$prof_mat))
    stop("No SymPortal profiles in this object.", call. = FALSE)
  md <- sp$metadata
  if (!sample_label %in% colnames(md)) sample_label <- "sample_uid"

  dom_col  <- paste0(model, "_dominant")
  conf_col <- switch(model, dmm = "dmm_certainty",
                     lda = "lda_max_prop", hdp = "hdp_max_prop")

  # Optional SymPortal-style names for topics/components (from name_topics()).
  # tname maps each component label (Topic_N / HDP_N) -> composition name.
  tname <- setNames(colnames(m$theta), colnames(m$theta))
  if (use_names) {
    nm <- name_topics(sp, model = model)
    tname <- setNames(nm$name, nm$group)
    # Guard: any component missing a name keeps its original label
    miss <- setdiff(colnames(m$theta), names(tname))
    if (length(miss) > 0) tname[miss] <- miss
    tname <- tname[colnames(m$theta)]
  }

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

  # === ROW 1: model assignment band (raw dominant, or screen-aware) ===========
  mix_uids <- character(0)   # samples flagged as genuine mixes (for row-2 boxes)
  if (!screen) {
    r1 <- data.frame(sample_uid = common,
                     group = factor(paste0(m$label_prefix, "_",
                                           as.integer(sub[[dom_col]])),
                                    levels = colnames(theta)))
    r1$.x <- xf(r1$sample_uid)
    grp_cols <- .topic_palette(colnames(theta), names_map = tname)
    row1_title <- sprintf("%s assignment", toupper(model))
    row1_name  <- toupper(model)
    r1_stacked <- FALSE
  } else {
    # Screen-aware: collapse intragenomic pairs into units; show each sample's
    # UNIT composition as a stack (pure = one colour, mix = multi-colour).
    scr <- suppressMessages(screen_intragenomic(sp, model = model))
    gmap <- setNames(scr$groups$group, scr$groups$component)
    units <- sort(unique(gmap))
    umat <- vapply(units, function(g)
      rowSums(theta[, names(gmap)[gmap == g], drop = FALSE]),
      numeric(nrow(theta)))
    unit_label <- vapply(units, function(g)
      paste(names(gmap)[gmap == g], collapse = "+"), character(1))
    colnames(umat) <- unit_label
    rownames(umat) <- rownames(theta)

    # Genuine-mix flag: a sample is boxed only if (a) no unit reaches
    # `majority`, AND (b) its top two units are connected by a "mixing"-verdict
    # topic pair. Falling below majority alone is not enough -- the co-occurring
    # units must be a confirmed genuine-mixing pair, not merely distinct or
    # ambiguous topics that happen to co-occur in one sample.
    unit_of <- gmap                                   # topic -> unit id
    mixing_pairs <- scr$pairs[scr$pairs$verdict == "mixing", , drop = FALSE]
    # Unordered unit-id pairs that are genuine mixing (drop within-unit pairs,
    # which would be intragenomic-merged topics, not between-unit mixing).
    mup <- lapply(seq_len(nrow(mixing_pairs)), function(i) {
      u <- sort(c(unit_of[[mixing_pairs$comp_a[i]]],
                  unit_of[[mixing_pairs$comp_b[i]]]))
      if (u[1] == u[2]) NULL else u
    })
    mup <- mup[!vapply(mup, is.null, logical(1))]
    mixing_unit_pairs <- if (length(mup)) unique(do.call(rbind, mup)) else
      matrix(numeric(0), ncol = 2)
    is_mixing_unit_pair <- function(u1, u2) {
      if (nrow(mixing_unit_pairs) == 0) return(FALSE)
      key <- sort(c(u1, u2))
      any(mixing_unit_pairs[, 1] == key[1] & mixing_unit_pairs[, 2] == key[2])
    }

    top_p <- apply(umat, 1, max)
    sub_majority <- rownames(umat)[top_p < majority]
    mix_uids <- Filter(function(uid) {
      ord2 <- order(umat[uid, ], decreasing = TRUE)
      u1 <- units[ord2[1]]; u2 <- units[ord2[2]]
      umat[uid, ord2[2]] >= 0.05 && is_mixing_unit_pair(u1, u2)
    }, sub_majority)
    mix_uids <- unlist(mix_uids)
    if (is.null(mix_uids)) mix_uids <- character(0)

    # Row 1 = single colour per sample (its dominant unit), EXCEPT confirmed
    # genuine mixes, which show their stacked unit composition. Pure/intragenomic
    # samples are one solid colour; only true admixtures reveal their mix here.
    dom_unit <- unit_label[apply(umat, 1, which.max)]
    names(dom_unit) <- rownames(umat)

    non_mix <- setdiff(rownames(umat), mix_uids)
    # Non-mix samples: a single full-height (prop = 1) row of the dominant unit.
    r1_solid <- data.frame(
      sample_uid = non_mix,
      group = dom_unit[non_mix],
      prop = 1, stringsAsFactors = FALSE)
    # Mix samples: stacked unit composition (renormalised for a clean 0-1 bar).
    if (length(mix_uids) > 0) {
      um <- umat[mix_uids, , drop = FALSE]
      um <- um / pmax(rowSums(um), .Machine$double.eps)
      r1_mix <- as.data.frame(um, check.names = FALSE)
      r1_mix$sample_uid <- rownames(um)
      r1_mix <- tidyr::pivot_longer(r1_mix, -"sample_uid", names_to = "group",
                                    values_to = "prop")
      r1_mix <- r1_mix[r1_mix$prop > 0, ]
      r1 <- rbind(r1_solid, as.data.frame(r1_mix))
    } else {
      r1 <- r1_solid
    }
    r1$.x <- xf(r1$sample_uid)
    r1$group <- factor(r1$group, levels = unit_label)
    # Unit -> name of its first constituent topic (for clade-based colouring)
    unit_name_map <- setNames(
      vapply(units, function(g) {
        first_topic <- names(gmap)[gmap == g][1]
        tname[[first_topic]]
      }, character(1)),
      unit_label)
    grp_cols <- .topic_palette(unit_label, names_map = unit_name_map)
    row1_title <- sprintf(
      "%s assignment (screen-aware: %d units; %d genuine-mix samples, boxed by combination)",
      toupper(model), length(unit_label), length(mix_uids))
    row1_name  <- "Unit"
    r1_stacked <- TRUE   # r1 now carries an explicit prop column in both cases
  }

  p1 <- ggplot2::ggplot(r1, ggplot2::aes(.data$.x, if (r1_stacked) .data$prop else 1,
                                         fill = .data$group)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = grp_cols, name = row1_name,
                               guide = "none") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = row1_title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank()) + x_theme()

  # === ROW 2: model admixture (boxed component proportions) ===================
  # `residual` controls how sub-1.0 theta rows (mass HDP assigns to no extracted
  # component) are shown: "show" caps bars with a grey Residual segment (honest
  # but visually ragged); "renormalise" rescales each row to sum to 1 (clean
  # bars; relative component proportions preserved, absolute residual dropped).
  theta_plot <- theta
  if (residual == "renormalise")
    theta_plot <- theta / pmax(rowSums(theta), .Machine$double.eps)
  a <- as.data.frame(theta_plot, check.names = FALSE)
  a$sample_uid <- rownames(theta_plot)
  if (residual == "show")
    a$Residual <- pmax(0, 1 - rowSums(theta))
  a <- tidyr::pivot_longer(a, -"sample_uid", names_to = "group",
                           values_to = "prop")
  a <- a[a$prop > 0, ]
  a$.x <- xf(a$sample_uid)
  lvls <- if (residual == "show") c(colnames(theta), "Residual") else colnames(theta)
  a$group <- factor(a$group, levels = lvls)
  # Topic-keyed palette: clade-grouped, deterministic (stable across plots)
  topic_cols <- .topic_palette(colnames(theta), names_map = tname)
  if (residual == "show") topic_cols <- c(topic_cols, Residual = "grey90")

  p2_labels <- c(tname, Residual = "Residual")
  p2 <- ggplot2::ggplot(a, ggplot2::aes(.data$.x, .data$prop,
                                        fill = .data$group)) +
    ggplot2::geom_col(width = 1, linewidth = 0) +
    ggplot2::scale_fill_manual(values = topic_cols, name = toupper(model),
                               labels = p2_labels,
                               guide = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = sprintf("%s admixture", toupper(model)),
                  x = NULL, y = "Proportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) + x_theme()

  # Box consistent mix combinations on ROW 1 (screen mode): group contiguous
  # runs of samples sharing the same unit combination and draw one box per run.
  if (screen && length(mix_uids) > 0) {
    # Combination signature per mix sample: sorted units above min_frac (0.05)
    combo_of <- vapply(mix_uids, function(uid) {
      present <- unit_label[umat[uid, ] >= 0.05]
      paste(sort(present), collapse = " + ")
    }, character(1))
    # Position each mix sample on the shared x-axis
    mix_lab <- md[[sample_label]][match(mix_uids, md$sample_uid)]
    mix_pos <- match(mix_lab, ord_lab)
    keep <- !is.na(mix_pos)
    mix_pos <- mix_pos[keep]; combo_of <- combo_of[keep]
    ord2 <- order(mix_pos)
    mix_pos <- mix_pos[ord2]; combo_of <- combo_of[ord2]

    if (length(mix_pos) > 0) {
      # Contiguous run = consecutive x-positions sharing the same combination
      brk <- c(TRUE, diff(mix_pos) != 1 | combo_of[-1] != combo_of[-length(combo_of)])
      run_id <- cumsum(brk)
      runs <- lapply(split(seq_along(mix_pos), run_id), function(ix) {
        data.frame(xmin = min(mix_pos[ix]) - 0.5,
                   xmax = max(mix_pos[ix]) + 0.5)
      })
      pbox <- do.call(rbind, runs)
      p1 <- p1 + ggplot2::geom_rect(
        data = pbox, inherit.aes = FALSE,
        mapping = ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                               ymin = 0, ymax = 1),
        fill = NA, colour = "grey50", linewidth = 0.5, linetype = "dashed")
    }
  }

  # === MIDDLE PANEL: raw ITS2 sequence composition (SymPortal colours) ========
  # The underlying data both the model and SymPortal are built from. Coloured
  # with sp$col_dict (SymPortal's own post-MED sequence palette). ALL sequences
  # are coloured; the legend is limited to the top `n_seqs` by abundance.
  cm_seq <- sp$count_mat[common, , drop = FALSE]
  rel_seq <- cm_seq / pmax(rowSums(cm_seq), .Machine$double.eps)
  all_seqs <- colnames(rel_seq)
  seq_rank <- sort(colSums(rel_seq), decreasing = TRUE)
  legend_seqs <- names(utils::head(seq_rank, n_seqs))   # legend breaks only

  seqdf <- as.data.frame(rel_seq, check.names = FALSE)
  seqdf$sample_uid <- rownames(rel_seq)
  seqdf <- tidyr::pivot_longer(seqdf, cols = dplyr::all_of(all_seqs),
                               names_to = "sequence", values_to = "rel")
  seqdf <- seqdf[seqdf$rel > 0, ]
  seqdf$.x <- xf(seqdf$sample_uid)
  seq_order <- names(seq_rank)                            # abundance order
  seqdf$sequence <- factor(seqdf$sequence, levels = seq_order)
  seqcols <- .build_colours(all_seqs, sp$col_dict)

  pseq <- ggplot2::ggplot(seqdf, ggplot2::aes(.data$.x, .data$rel,
                                              fill = .data$sequence)) +
    ggplot2::geom_col(width = 1, linewidth = 0) +
    ggplot2::scale_fill_manual(values = seqcols, name = "ITS2 sequence",
                               breaks = legend_seqs,
                               guide = ggplot2::guide_legend(ncol = 3)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "Raw ITS2 sequences", x = NULL, y = "Proportion") +
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

  # Colour EVERY profile from SymPortal's own palette (prof_col_dict, keyed by
  # UID) so the plot matches SymPortal visually; the legend is then restricted
  # to the profiles that form a dominant assignment (below) to stay legible.
  pcd <- sp$prof_col_dict
  all_names <- sort(unique(prof_names))
  prof_cols <- setNames(rep(NA_character_, length(all_names)), all_names)
  if (!is.null(pcd)) {
    for (uid in colnames(prof)) {
      nm <- pname[as.character(uid)]
      if (!is.na(nm) && nm %in% names(prof_cols) &&
          as.character(uid) %in% names(pcd))
        prof_cols[nm] <- pcd[[as.character(uid)]]
    }
  }
  miss <- names(prof_cols)[is.na(prof_cols)]
  if (length(miss) > 0)
    prof_cols[miss] <- scales::hue_pal()(length(miss))

  # Legend breaks: profiles that are some sample's dominant assignment, capped
  # at the `max_legend` most frequent so the legend stays on-panel for large
  # datasets. All profiles remain coloured; only the legend is capped.
  # Samples with no SymPortal profile (all-zero prof row) are flagged "None"
  # rather than given a spurious which.max dominant.
  has_prof <- rowSums(prof) > 0
  sp_dom_uid <- rep(NA_character_, nrow(prof))
  sp_dom_uid[has_prof] <-
    colnames(prof)[apply(prof[has_prof, , drop = FALSE], 1, which.max)]
  dom_names <- ifelse(is.na(sp_dom_uid), "None", pname[as.character(sp_dom_uid)])
  all_names_n <- c(all_names, "None")
  prof_cols_n <- c(prof_cols, None = "grey90")
  dom_freq <- sort(table(dom_names[dom_names != "None"]), decreasing = TRUE)
  legend_breaks <- names(utils::head(dom_freq, max_legend))

  # === ROW 3: SymPortal dominant-profile assignment band ======================
  r3 <- data.frame(sample_uid = common,
                   profile = factor(dom_names, levels = all_names_n))
  r3$.x <- xf(r3$sample_uid)

  p3 <- ggplot2::ggplot(r3, ggplot2::aes(.data$.x, 1, fill = .data$profile)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::scale_fill_manual(values = prof_cols_n, name = "SymPortal profile",
                               breaks = legend_breaks,
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
  sdf$profile <- factor(pname[as.character(sdf$profile_uid)], levels = all_names_n)
  sdf$.x <- xf(sdf$sample_uid)

  p4 <- ggplot2::ggplot(sdf, ggplot2::aes(.data$.x, .data$prop,
                                          fill = .data$profile)) +
    ggplot2::geom_col(width = 1, linewidth = 0) +
    ggplot2::scale_fill_manual(values = prof_cols_n, guide = "none") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(title = "SymPortal mix", x = NULL, y = "Proportion") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) +
    x_theme(bottom = TRUE)

  # Overlay the group boxes on every panel so they line up into columns
  p1 <- p1 + box_layer()
  p2 <- p2 + box_layer()
  pseq <- pseq + box_layer()
  p3 <- p3 + box_layer()
  p4 <- p4 + box_layer()

  patchwork::wrap_plots(list(p1, p2, pseq, p3, p4), ncol = 1,
                        heights = c(0.4, 1, 1, 0.4, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("%s vs SymPortal: assignment and admixture (%d samples)",
                      toupper(model), length(common)),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold")))
}

# ==============================================================================
# plot-lda.R — LDA-specific diagnostics
# ==============================================================================

#' Plot the LDA K-selection search
#'
#' Four panels across the searched K: maximum pairwise topic cosine similarity
#' (with the redundancy threshold marked), number of redundant pairs, percent
#' of topics with pure samples, and perplexity. The selected K (largest with no
#' redundant pairs) is marked. Requires `fit_lda(sp)` to have been run with
#' adaptive selection (i.e. `k = NULL`).
#'
#' @param sp A `symbayes` object with LDA fitted adaptively.
#' @return A \pkg{patchwork} composition of four \pkg{ggplot2} panels.
#' @examples
#' \dontrun{ lda_plot_ksearch(sp) }
#' @export
lda_plot_ksearch <- function(sp) {
  m <- .get_membership(sp, "lda")
  df <- m$meta$selection
  if (is.null(df))
    stop("No K-search table; fit_lda() was called with a fixed k.",
         call. = FALSE)
  optk <- m$k

  base <- function(y, title, hline = NULL) {
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$K, .data[[y]])) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_vline(xintercept = optk, linetype = "dashed",
                          colour = "blue") +
      ggplot2::labs(title = title, y = title) +
      ggplot2::theme_minimal(base_size = 11)
    if (!is.null(hline))
      p <- p + ggplot2::geom_hline(yintercept = hline, linetype = "dashed",
                                   colour = "red", linewidth = 0.4)
    p
  }

  p1 <- base("max_cosine", "Max pairwise cosine", hline = 0.85)
  p2 <- base("n_redundant", "Redundant topic pairs")
  p3 <- base("pct_pure", "% topics with pure samples")
  p4 <- base("perplexity", "Perplexity")

  patchwork::wrap_plots(list(p1, p2, p3, p4), ncol = 2) +
    patchwork::plot_annotation(
      title = sprintf("LDA K selection: optimal K = %d (no redundant topics)",
                      optk),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13)))
}


#' Plot the LDA topic cosine-similarity matrix
#'
#' Heatmap of pairwise topic cosine similarity, annotated with values. Pairs
#' above ~0.85 indicate redundant (over-split) topics.
#'
#' @param sp A `symbayes` object with LDA fitted.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ lda_plot_similarity(sp) }
#' @export
lda_plot_similarity <- function(sp) {
  m <- .get_membership(sp, "lda")
  sim <- m$meta$similarity
  K <- nrow(sim)

  long <- as.data.frame(as.table(sim))
  names(long) <- c("A", "B", "Cosine")
  long$A <- factor(long$A, levels = rev(rownames(sim)))
  long$B <- factor(long$B, levels = colnames(sim))

  ggplot2::ggplot(long, ggplot2::aes(.data$B, .data$A, fill = .data$Cosine)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(
      label = ifelse(.data$Cosine > 0.01, sprintf("%.2f", .data$Cosine), "")),
      size = 2.8) +
    ggplot2::scale_fill_gradient2(low = "white", mid = "#FDE725",
                                  high = "#440154", midpoint = 0.5,
                                  limits = c(0, 1), name = "Cosine") +
    ggplot2::labs(
      title = sprintf("LDA topic cosine similarity (K = %d)", K),
      subtitle = sprintf("Max off-diagonal = %.2f", max(sim[upper.tri(sim)])),
      x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_equal()
}

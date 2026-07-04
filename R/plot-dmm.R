# ==============================================================================
# plot-dmm.R — DMM-specific diagnostics
# ==============================================================================

#' Plot DMM model selection (Laplace, AIC, BIC vs K)
#'
#' Shows the three information criteria across the candidate component numbers,
#' with the Laplace-selected K marked. Laplace (the model evidence
#' approximation) is the primary criterion.
#'
#' @param sp A `symbayes` object with DMM fitted.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ dmm_plot_selection(sp) }
#' @export
dmm_plot_selection <- function(sp) {
  m <- .get_membership(sp, "dmm")
  mf <- m$meta$model_fit
  long <- tidyr::pivot_longer(mf, -"K", names_to = "Criterion",
                              values_to = "Score")
  ggplot2::ggplot(long, ggplot2::aes(.data$K, .data$Score)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(ggplot2::vars(.data$Criterion), scales = "free_y") +
    ggplot2::geom_vline(xintercept = m$k, linetype = "dashed",
                        colour = "red") +
    ggplot2::labs(title = "DMM model selection",
                  subtitle = sprintf("Selected K = %d (Laplace)", m$k),
                  x = "K (components)", y = "Score (lower = better)") +
    ggplot2::theme_minimal(base_size = 12)
}

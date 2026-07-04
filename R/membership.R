# ==============================================================================
# membership.R — standardised membership schema and accessors
# ==============================================================================
# Every fit_*() writes a membership object with a common schema into sp$<model>.
# Generalist functions read models through .get_membership() so they never need
# to know which model produced the object.
# ==============================================================================

#' Construct a standardised membership object
#'
#' Internal constructor used by the `fit_*()` functions. All three models emit
#' this same structure so generalist plots can consume any of them.
#'
#' @param theta samples x groups membership matrix (rows sum to ~1)
#' @param beta groups x sequences exemplar matrix
#' @param dominant named integer/factor of per-sample dominant group
#' @param k number of groups
#' @param model one of "dmm", "lda", "hdp"
#' @param label_prefix prefix for group labels (e.g. "Comp", "Topic", "HDP")
#' @param meta named list of model-specific extras
#' @return A list of class `symbayes_membership`
#' @keywords internal
.membership <- function(theta, beta, dominant, k, model,
                        label_prefix = "Group", meta = list()) {

  stopifnot(is.matrix(theta), is.matrix(beta))

  # Standardise group labels on theta columns and beta rows
  grp_labels <- paste0(label_prefix, "_", seq_len(k))
  if (ncol(theta) == k) colnames(theta) <- grp_labels
  if (nrow(beta)  == k) rownames(beta)  <- grp_labels

  structure(
    list(
      theta        = theta,
      beta         = beta,
      dominant     = dominant,
      k            = k,
      model        = model,
      label_prefix = label_prefix,
      meta         = meta
    ),
    class = "symbayes_membership"
  )
}

#' Retrieve a fitted model's membership object from an `sp` object
#'
#' @param sp A `symbayes` object
#' @param model One of "dmm", "lda", "hdp"
#' @return The membership object
#' @keywords internal
.get_membership <- function(sp, model) {
  model <- match.arg(model, c("dmm", "lda", "hdp"))
  m <- sp[[model]]
  if (is.null(m)) {
    stop(sprintf("No %s fit in this object. Run fit_%s() first.",
                 toupper(model), model), call. = FALSE)
  }
  m
}

#' Check an object is a symbayes object
#' @keywords internal
.check_sp <- function(sp) {
  if (!inherits(sp, "symbayes"))
    stop("Expected a 'symbayes' object from import().", call. = FALSE)
  invisible(TRUE)
}

#' Print method for symbayes objects
#' @param x A `symbayes` object
#' @param ... Ignored
#' @export
print.symbayes <- function(x, ...) {
  cat("<symbayes>\n")
  cat(sprintf("  samples:   %d\n", nrow(x$count_mat)))
  cat(sprintf("  sequences: %d%s\n", ncol(x$count_mat),
              if (isTRUE(x$filtered)) " (filtered)" else " (unfiltered)"))
  if (!is.null(x$prof_mat))
    cat(sprintf("  SymPortal profiles: %d\n", ncol(x$prof_mat)))
  fitted <- intersect(c("dmm", "lda", "hdp"), names(x))
  if (length(fitted) > 0) {
    cat("  fitted models:\n")
    for (m in fitted) {
      cat(sprintf("    %-4s k = %d\n", toupper(m), x[[m]]$k))
    }
  } else {
    cat("  fitted models: none (run fit_dmm/fit_lda/fit_hdp)\n")
  }
  invisible(x)
}

#' Print method for membership objects
#' @param x A `symbayes_membership` object
#' @param ... Ignored
#' @export
print.symbayes_membership <- function(x, ...) {
  cat(sprintf("<symbayes_membership: %s, k = %d>\n", toupper(x$model), x$k))
  cat(sprintf("  theta: %d samples x %d groups\n",
              nrow(x$theta), ncol(x$theta)))
  cat(sprintf("  beta:  %d groups x %d sequences\n",
              nrow(x$beta), ncol(x$beta)))
  if (length(x$meta) > 0)
    cat(sprintf("  meta:  %s\n", paste(names(x$meta), collapse = ", ")))
  invisible(x)
}

# ==============================================================================
# utils.R — internal helpers
# ==============================================================================

#' Infer Symbiodiniaceae clade from sequence names
#'
#' SymPortal names sequences by clade (leading letter, e.g. `C3`, `D1`, `A6b`)
#' or by clade suffix for named intragenomic variants (e.g. `15443_A`).
#'
#' @param seq_names Character vector of sequence names
#' @return Named character vector of clade letters (A-D, F, G) or NA
#' @keywords internal
.infer_clade <- function(seq_names) {
  clade <- rep(NA_character_, length(seq_names))
  for (cl in c("A", "B", "C", "D", "F", "G")) {
    clade[grepl(paste0("^", cl), seq_names)] <- cl
    clade[grepl(paste0("_", cl, "$"), seq_names)] <- cl
  }
  setNames(clade, seq_names)
}

#' Build a colour vector from a SymPortal colour dictionary plus fallback
#'
#' @param seq_names Character vector of sequence names to colour
#' @param col_dict Named list/vector mapping sequence names to hex colours
#'   (from SymPortal's `color_dict_post_med.json`); may be `NULL`
#' @return Named character vector of hex colours, including a grey `Other`
#' @keywords internal
.build_colours <- function(seq_names, col_dict = NULL) {
  if (!is.null(col_dict)) {
    cols <- unlist(col_dict[seq_names])
    missing <- seq_names[!seq_names %in% names(col_dict)]
    if (length(missing) > 0) {
      fb <- setNames(rep("#AAAAAA", length(missing)), missing)
      cols <- c(cols, fb)
    }
  } else {
    cols <- setNames(scales::hue_pal()(length(seq_names)), seq_names)
  }
  c(cols, Other = "#D4D4D4")
}

#' Convert a count matrix to a topicmodels DocumentTermMatrix
#'
#' Builds the DTM directly as a `slam::simple_triplet_matrix` with the
#' `DocumentTermMatrix` class and weighting attribute, avoiding a hard
#' dependency on `tm::as.DocumentTermMatrix`.
#'
#' @param count_mat Integer sample x sequence matrix
#' @return A `DocumentTermMatrix`
#' @keywords internal
.count_to_dtm <- function(count_mat) {
  rs <- rowSums(count_mat)
  cs <- colSums(count_mat)
  mat <- count_mat[rs > 0, cs > 0, drop = FALSE]

  ij <- which(mat > 0, arr.ind = TRUE)
  stm <- slam::simple_triplet_matrix(
    i = ij[, 1], j = ij[, 2], v = as.numeric(mat[ij]),
    nrow = nrow(mat), ncol = ncol(mat),
    dimnames = list(Docs = rownames(mat), Terms = colnames(mat))
  )
  class(stm) <- c("DocumentTermMatrix", "simple_triplet_matrix")
  attr(stm, "weighting") <- c("term frequency", "tf")
  stm
}

#' Shannon entropy of a proportion vector (log2)
#' @param p Numeric vector of proportions
#' @return Scalar entropy in bits
#' @keywords internal
.entropy <- function(p) {
  p <- p[p > 0]
  if (length(p) == 0) return(0)
  -sum(p * log2(p))
}

#' Standard clade colour palette
#' @keywords internal
.clade_palette <- function() {
  c(A = "#FF8C00", B = "#8B4513", C = "#20B2AA",
    D = "#9370DB", F = "#CD853F", G = "#808000")
}

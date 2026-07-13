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

#' Stable, distinguishable colour palette for topics / components
#'
#' Assigns each topic/component an evenly-spaced hue across the full colour
#' wheel so all topics are maximally distinguishable, while keeping the
#' assignment DETERMINISTIC by label (same label -> same colour across plots,
#' models, and runs). Hue order is set by a deterministic hash of the label, so
#' colours are stable but spread out rather than clustered.
#'
#' Clade is not used to anchor hue: in datasets dominated by one clade,
#' clade-anchoring collapses most topics onto near-identical colours. This
#' palette prioritises telling topics apart. (Clade grouping remains available
#' via [.clade_palette] for clade-level fills such as [plot_top_seqs].)
#'
#' @param labels Character vector of topic/component labels.
#' @param names_map Accepted for API compatibility; unused here.
#' @return A named character vector of hex colours, one per label.
#' @keywords internal
.topic_palette <- function(labels, names_map = NULL) {
  n <- length(labels)
  if (n == 0) return(setNames(character(0), character(0)))

  # Deterministic hash of a string -> value in [0, 1) (overflow-safe).
  hash01 <- function(s) {
    if (is.na(s) || nchar(s) == 0) return(0.5)
    bytes <- as.integer(charToRaw(s))
    h <- 0
    for (byte in bytes) h <- (h * 131 + byte) %% 1000003
    h / 1000003
  }

  # Order labels by their hash, then assign evenly-spaced hues in that order.
  # Deterministic (hash-driven) yet maximally separated (even spacing).
  hv <- vapply(labels, hash01, numeric(1))
  ord <- order(hv)
  hues <- seq(0, 1, length.out = n + 1)[seq_len(n)]
  # Slight deterministic per-label jitter in saturation/value for extra contrast
  cols <- character(n)
  for (k in seq_len(n)) {
    i <- ord[k]
    s <- 0.65 + 0.2 * ((hv[i] * 7) %% 1)     # 0.65-0.85
    v <- 0.80 + 0.15 * ((hv[i] * 13) %% 1)   # 0.80-0.95
    cols[i] <- grDevices::hsv(hues[k], min(1, s), min(1, v))
  }
  setNames(cols, labels)
}

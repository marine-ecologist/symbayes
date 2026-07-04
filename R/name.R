# ==============================================================================
# name.R — SymPortal-style names for model topics / components
# ==============================================================================
# Convenience labelling: render each group's exemplar (beta) in SymPortal's
# naming convention so model groups are directly comparable to SymPortal
# profiles. The name is for DISPLAY and comparison only -- a group's identity is
# its full distribution and its samples' membership, not this string. Two groups
# sharing a backbone but differing in minor variants may receive similar names;
# that is a property of the naming convention (and part of why SymPortal
# over-splits), not of the groups themselves.
# ==============================================================================

#' Name a single exemplar vector in SymPortal style
#'
#' `/` separates co-dominant sequences (within `codominant_ratio` of the top);
#' `-` appends less-abundant associates in descending order, down to
#' `min_weight`.
#'
#' @param beta_row Named numeric vector (sequence -> weight), sums to ~1.
#' @param codominant_ratio Sequences whose weight is >= this fraction of the
#'   single largest are treated as co-dominant and joined with `/` (default 0.5).
#' @param min_weight Minimum weight to include a sequence at all (default 0.05).
#' @param max_seqs Maximum sequences in the name (default 6).
#' @return A SymPortal-style name string.
#' @keywords internal
.name_exemplar <- function(beta_row, codominant_ratio = 0.5,
                           min_weight = 0.05, max_seqs = 6) {
  b <- sort(beta_row[beta_row >= min_weight], decreasing = TRUE)
  if (length(b) == 0) {                 # nothing above threshold: use the top 1
    b <- sort(beta_row, decreasing = TRUE)[1]
  }
  b <- utils::head(b, max_seqs)
  top <- b[1]

  # Co-dominant set: sequences within codominant_ratio of the largest
  is_codom <- b >= top * codominant_ratio
  codom <- names(b)[is_codom]
  assoc <- names(b)[!is_codom]

  name <- paste(codom, collapse = "/")
  if (length(assoc) > 0)
    name <- paste0(name, "-", paste(assoc, collapse = "-"))
  name
}

#' Name a fitted model's topics/components in SymPortal style
#'
#' Renders each group's exemplar sequence distribution (`beta`) as a
#' SymPortal-style name, so model groups can be compared to SymPortal profiles
#' in the same convention. Returns a data frame mapping the group label to its
#' name, dominant clade, and the exemplar's Shannon entropy (a diversity /
#' potential-over-naming flag: high-entropy exemplars produce long names).
#'
#' The names are a labelling convenience for display and comparison. A group's
#' identity is its full distribution and its samples' membership, not this
#' string. Groups sharing a backbone but differing in minor variants can receive
#' similar names -- the same mechanism behind SymPortal over-splitting -- so do
#' not treat similar names as evidence groups are redundant (use the cosine
#' similarity from [lda_plot_similarity] or the HDP component distances for
#' that).
#'
#' @param sp A `symbayes` object with the chosen model fitted.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param codominant_ratio Co-dominance cutoff for `/` (default 0.5).
#' @param min_weight Minimum sequence weight to include (default 0.05).
#' @param max_seqs Maximum sequences per name (default 6).
#' @param assign If `TRUE`, also writes the names onto `sp` as
#'   `sp$<model>$meta$names` and returns the updated `sp`; if `FALSE` (default)
#'   returns the naming data frame.
#' @return A data frame (group, name, clade, entropy, n_samples), or the updated
#'   `sp` object if `assign = TRUE`.
#' @examples
#' \dontrun{
#' name_topics(sp, model = "lda")
#' sp <- name_topics(sp, model = "hdp", assign = TRUE)
#' sp$hdp$meta$names
#' }
#' @export
name_topics <- function(sp, model = c("dmm", "lda", "hdp"),
                        codominant_ratio = 0.5, min_weight = 0.05,
                        max_seqs = 6, assign = FALSE) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  beta <- m$beta
  sc <- sp$seq_clades
  dom_col <- paste0(model, "_dominant")

  n_per <- table(factor(sp$metadata[[dom_col]], levels = seq_len(m$k)))

  rows <- lapply(seq_len(m$k), function(k) {
    nm <- .name_exemplar(beta[k, ], codominant_ratio, min_weight, max_seqs)
    top_seq <- names(which.max(beta[k, ]))
    data.frame(
      group    = rownames(beta)[k],
      name     = nm,
      clade    = if (!is.na(sc[top_seq])) sc[top_seq] else "?",
      entropy  = round(.entropy(beta[k, ]), 3),
      n_samples = as.integer(n_per[k]),
      stringsAsFactors = FALSE)
  })
  naming <- do.call(rbind, rows)

  if (assign) {
    sp[[model]]$meta$names <- setNames(naming$name, naming$group)
    return(sp)
  }
  naming
}


#' Match model topics/components to SymPortal profiles by composition
#'
#' String-matching topic names to SymPortal profile names fails: the same
#' community can be written `C22-C3-C22ad` (topic) and `C22-C22ad-C3-C22ah`
#' (SymPortal) -- overlapping sequences, different order and set. This function
#' instead matches by *composition*: it represents each SymPortal profile as a
#' sequence vector (from the profile's defining sequences, via `prof_meta` if
#' available, else from the mean sequence composition of samples dominant for
#' that profile) and computes cosine similarity to each topic's exemplar
#' (`beta`). Each topic is matched to its most similar profile.
#'
#' @param sp A `symbayes` object with the model fitted and SymPortal profiles
#'   present.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param top_n Number of best-matching profiles to report per topic (default 3).
#' @return A data frame: one row per topic-profile match, with `group`,
#'   `topic_name`, `sp_profile`, `cosine`, and match rank.
#' @examples
#' \dontrun{ match_profiles(sp, model = "lda") }
#' @export
match_profiles <- function(sp, model = c("dmm", "lda", "hdp"), top_n = 3) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  if (is.null(sp$prof_mat))
    stop("No SymPortal profiles in this object.", call. = FALSE)

  beta <- m$beta                       # topics x sequences
  seqs <- colnames(beta)
  topic_names <- name_topics(sp, model)$name

  # --- Build a sequence-space vector for each SymPortal profile ---------------
  # Preferred: mean sequence composition of samples whose dominant profile is p.
  prof <- sp$prof_mat[rownames(sp$count_mat), , drop = FALSE]
  prof <- prof[, colSums(prof) > 0, drop = FALSE]
  rel  <- sp$count_mat / rowSums(sp$count_mat)
  dom_prof <- colnames(prof)[apply(prof, 1, which.max)]

  pm <- sp$prof_meta
  pname <- if (!is.null(pm))
    setNames(pm$`ITS2 type profile`, as.character(pm$`ITS2 type profile UID`))
  else setNames(colnames(prof), colnames(prof))

  prof_vecs <- lapply(colnames(prof), function(puid) {
    samps <- rownames(rel)[dom_prof == puid]
    if (length(samps) == 0) return(NULL)
    v <- colMeans(rel[samps, , drop = FALSE])
    v[seqs]                            # align to beta's sequence order
  })
  names(prof_vecs) <- colnames(prof)
  prof_vecs <- prof_vecs[!vapply(prof_vecs, is.null, logical(1))]
  prof_mat <- do.call(rbind, prof_vecs)
  prof_mat[is.na(prof_mat)] <- 0

  # --- Cosine similarity: each topic (beta row) vs each profile vector --------
  cos_sim <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (d == 0) 0 else sum(a * b) / d
  }

  rows <- lapply(seq_len(m$k), function(k) {
    sims <- vapply(seq_len(nrow(prof_mat)),
                   function(i) cos_sim(beta[k, ], prof_mat[i, ]), numeric(1))
    names(sims) <- rownames(prof_mat)
    ord <- order(sims, decreasing = TRUE)[seq_len(min(top_n, length(sims)))]
    data.frame(
      group      = rownames(beta)[k],
      topic_name = topic_names[k],
      rank       = seq_along(ord),
      sp_profile = pname[rownames(prof_mat)[ord]],
      cosine     = round(sims[ord], 3),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Relabel dominant assignment with a minimum-dominance threshold
#'
#' The hard `dominant` label is `which.max(theta)` -- a plurality that can be as
#' low as ~1/k. For near-boundary samples this label is fragile and overstates
#' confidence (a 0.35/0.34 split is labelled definitively by its 0.35 winner).
#' This function recomputes the dominant label so that samples whose top
#' component is below `min_dominance` are labelled `"mixed"` instead of forced
#' to a plurality winner. It changes only the *label* (the `<model>_dominant`
#' metadata column and the membership `dominant` factor); `theta` -- the model's
#' actual fractional output -- is untouched.
#'
#' There is no "correct" threshold: this is a reporting-honesty choice, not a
#' model parameter. `0.5` (call a sample mixed unless one component holds an
#' outright majority) is a natural, explainable default. Note SymPortal has no
#' equivalent control, because it never produces a proportion to threshold --
#' it assigns each sample one profile (possibly a compound name that already
#' absorbs the mixing). This function is symbayes acknowledging admixture
#' explicitly where SymPortal cannot.
#'
#' @param sp A `symbayes` object with the model fitted.
#' @param model One of `"dmm"`, `"lda"`, `"hdp"`.
#' @param min_dominance Minimum top-component proportion to keep a hard label;
#'   below this the sample is labelled `"mixed"` (default 0.5).
#' @return The updated `sp` object. Adds a `<model>_dominant_thr` factor to
#'   metadata (levels = group labels plus `"mixed"`) and stores the threshold
#'   in `sp$<model>$meta$min_dominance`. The original `<model>_dominant` is left
#'   intact for backwards compatibility.
#' @examples
#' \dontrun{
#' sp <- set_dominance(sp, model = "hdp", min_dominance = 0.5)
#' table(sp$metadata$hdp_dominant_thr)   # how many samples are "mixed"
#' }
#' @export
set_dominance <- function(sp, model = c("dmm", "lda", "hdp"),
                          min_dominance = 0.5) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta

  top_grp  <- apply(theta, 1, which.max)
  top_prop <- apply(theta, 1, max)
  grp_labels <- colnames(theta)

  lab <- ifelse(top_prop >= min_dominance, grp_labels[top_grp], "mixed")
  lab <- factor(lab, levels = c(grp_labels, "mixed"))
  names(lab) <- rownames(theta)

  # Write a thresholded label column alongside the original
  col <- paste0(model, "_dominant_thr")
  sp$metadata[[col]] <- lab[sp$metadata$sample_uid]

  n_mixed <- sum(lab == "mixed", na.rm = TRUE)
  message(sprintf("%s: %d/%d samples below dominance %.2f labelled 'mixed'",
                  toupper(model), n_mixed, nrow(theta), min_dominance))

  sp[[model]]$meta$min_dominance <- min_dominance
  sp
}

# ==============================================================================
# validate-mixing.R — depth-aware intragenomic-vs-mixing discrimination
# ==============================================================================
# The central confound in claiming a probabilistic model detects "more mixing"
# than SymPortal: the ITS2 marker is multi-copy, so a single symbiont genotype
# produces several sequences in FIXED intragenomic ratios. A sample that is
# biologically one symbiont can therefore look admixed across topics.
#
# Intragenomic variation and genuine community mixing are distinguishable by
# RATIO CONSISTENCY across samples -- but only once read depth is accounted for,
# because at low depth sampling noise inflates apparent ratio variability.
#
# This tests, for each co-occurring pair of components, whether the observed
# cross-sample variance in their ratio EXCEEDS the multinomial sampling
# variance expected at each sample's actual read depth (overdispersion). Ratio
# consistent with sampling noise => fixed ratio => intragenomic. Overdispersed
# => genuinely variable ratio => community mixing.
# ==============================================================================

#' Assign defining sequences to each component from its exemplar
#'
#' Each component's defining sequences are those where it holds the largest
#' share of that sequence's weight across components (argmax over the beta
#' columns), restricted to sequences the component emphasises (its own top
#' sequences). Sequences that are not clearly owned by one component are left
#' unassigned, so a shared backbone sequence does not corrupt a pair's ratio.
#'
#' @param beta components x sequences exemplar matrix.
#' @param top_n Consider each component's top-N sequences as candidates
#'   (default 15).
#' @param owner_frac A sequence is "owned" by a component only if that component
#'   holds at least this fraction of the sequence's total weight across
#'   components (default 0.5 = clear majority).
#' @return A named character vector: sequence -> owning component (or dropped).
#' @keywords internal
.defining_seqs <- function(beta, top_n = 15, owner_frac = 0.5) {
  K <- nrow(beta)
  # Candidate sequences: union of each component's top_n
  cand <- unique(unlist(lapply(seq_len(K), function(k)
    names(utils::head(sort(beta[k, ], decreasing = TRUE), top_n)))))
  bc <- beta[, cand, drop = FALSE]
  col_tot <- colSums(bc)
  owner_idx <- apply(bc, 2, which.max)
  owner_share <- vapply(seq_along(cand), function(j)
    bc[owner_idx[j], j] / col_tot[j], numeric(1))
  keep <- owner_share >= owner_frac
  setNames(rownames(beta)[owner_idx[keep]], cand[keep])
}


#' Depth-aware overdispersion test (rebuilt on real read counts)
#'
#' Rebuilt to use **actual sequence read counts** rather than `theta * depth`.
#' Each component is represented by its *defining sequences* (those it clearly
#' owns in its exemplar). For a pair (A, B), `k_i` = sample i's summed reads of
#' A's defining sequences and `n_i` = reads of A's + B's defining sequences --
#' genuine integer read counts, so the binomial sampling model is legitimate.
#'
#' Overdispersion is estimated as the beta-binomial intra-class correlation
#' \eqn{\rho} (method of moments), which bounds in \eqn{[0, 1)} and properly
#' separates true cross-sample ratio variance from within-sample sampling
#' variance. \eqn{\rho \approx 0} = fixed ratio (consistent with intragenomic
#' variation); \eqn{\rho} clearly above 0 = variable ratio (consistent with
#' community mixing). A dispersion index `phi` (Pearson X^2 / df) is reported
#' alongside as a cross-check -- now well-behaved because `n_i` are real counts.
#'
#' @param sp A `symbayes` object with the chosen model fitted.
#' @param model One of `"lda"` or `"hdp"`.
#' @param min_frac Component presence threshold to count a sample as co-occurring
#'   (default 0.05).
#' @param min_samples Minimum co-occurring samples to test a pair (default 5).
#' @param min_pair_reads Minimum `n_i` (pair reads) to include a sample
#'   (default 100).
#' @param rho_mixing Beta-binomial rho above which a pair is called `"mixing"`
#'   (default 0.05). Between a small floor and this, `"ambiguous"`.
#' @param top_n,owner_frac Passed to the internal defining-sequence assignment.
#' @return Invisibly, a data frame: `comp_a`, `comp_b`, `n_samples`,
#'   `n_seqs_a`, `n_seqs_b`, `mean_ratio`, `rho`, `phi`, `p_value`, `call`.
#' @examples
#' \dontrun{
#' od <- test_overdispersion(sp, model = "lda")
#' subset(od, call == "mixing")
#' }
#' @export
test_overdispersion <- function(sp, model = c("lda", "hdp"),
                                min_frac = 0.05, min_samples = 5,
                                min_pair_reads = 100, rho_mixing = 0.05,
                                top_n = 15, owner_frac = 0.5) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta
  cm <- sp$count_mat[rownames(theta), , drop = FALSE]
  K <- ncol(theta)
  grp <- colnames(theta)

  # Assign defining sequences to components
  owner <- .defining_seqs(m$beta, top_n = top_n, owner_frac = owner_frac)
  seqs_by_comp <- split(names(owner), owner)

  # Per-sample summed reads of each component's defining sequences
  comp_reads <- sapply(grp, function(g) {
    s <- seqs_by_comp[[g]]
    if (is.null(s) || length(s) == 0) return(rep(0, nrow(cm)))
    rowSums(cm[, intersect(s, colnames(cm)), drop = FALSE])
  })
  rownames(comp_reads) <- rownames(cm)

  # Beta-binomial rho by method of moments for a vector of (k, n)
  bb_rho <- function(k, n) {
    m_s <- length(n)
    p_hat <- sum(k) / sum(n)
    if (p_hat <= 0 || p_hat >= 1) return(c(rho = NA, phi = NA, chisq = NA))
    # Pearson chi-square (binomial null)
    ev <- n * p_hat * (1 - p_hat)
    chisq <- sum((k - n * p_hat)^2 / ev)
    phi <- chisq / (m_s - 1)
    # Method-of-moments rho (Crowder / Williams): S = sum (k - n p)^2 / (p(1-p))
    # E[S] under bb = sum n (1 + (n-1) rho). Solve for rho.
    S <- sum((k - n * p_hat)^2) / (p_hat * (1 - p_hat))
    denom <- sum(n * (n - 1))
    rho <- if (denom > 0) max(0, (S - sum(n)) / denom) else NA
    c(rho = rho, phi = phi, chisq = chisq)
  }

  rows <- list()
  for (a in seq_len(K - 1)) for (b in (a + 1):K) {
    both <- theta[, a] >= min_frac & theta[, b] >= min_frac
    if (sum(both) < min_samples) next
    ka <- comp_reads[both, grp[a]]
    kb <- comp_reads[both, grp[b]]
    n_i <- ka + kb
    keep <- n_i >= min_pair_reads
    if (sum(keep) < min_samples) next
    ka <- ka[keep]; n_i <- n_i[keep]

    est <- bb_rho(ka, n_i)
    if (is.na(est["rho"])) next
    df <- length(n_i) - 1
    p_value <- stats::pchisq(est["chisq"], df, lower.tail = FALSE)

    call <- if (est["rho"] >= rho_mixing) "mixing"
    else if (est["rho"] > 0.01) "ambiguous"
    else "intragenomic"

    rows[[length(rows) + 1]] <- data.frame(
      comp_a = grp[a], comp_b = grp[b],
      n_samples = length(n_i),
      n_seqs_a = length(seqs_by_comp[[grp[a]]] %||% character(0)),
      n_seqs_b = length(seqs_by_comp[[grp[b]]] %||% character(0)),
      mean_ratio = round(sum(ka) / sum(n_i), 3),
      rho = round(est["rho"], 4),
      phi = round(est["phi"], 1),
      p_value = signif(p_value, 3),
      call = call, stringsAsFactors = FALSE)
  }

  if (length(rows) == 0) {
    message("No component pairs co-occur in enough samples with owned sequences.")
    return(invisible(data.frame()))
  }
  out <- do.call(rbind, rows)
  out <- out[order(-out$rho), ]
  rownames(out) <- NULL

  cat(sprintf("Overdispersion test (%s), real read counts: %d pairs tested\n",
              toupper(model), nrow(out)))
  cat(sprintf("  mixing (rho >= %.2f):        %d\n", rho_mixing,
              sum(out$call == "mixing")))
  cat(sprintf("  ambiguous (0.01 < rho < %.2f): %d\n", rho_mixing,
              sum(out$call == "ambiguous")))
  cat(sprintf("  intragenomic (rho <= 0.01):  %d\n",
              sum(out$call == "intragenomic")))
  cat("  rho = beta-binomial overdispersion (0 = fixed ratio; higher = mixing).\n")
  cat("  psbA (single-copy) remains the definitive arbiter.\n")
  invisible(out)
}

# Null-coalescing helper
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Plot the overdispersion result for one component pair (real read counts)
#'
#' Each co-occurring sample's observed A-share from real defining-sequence
#' reads, vs its pair read depth, with the binomial fixed-ratio null envelope
#' (+/-2 SE around the pooled ratio). Points tracking the band at all depths =
#' fixed ratio (intragenomic); systematic spread beyond it, especially at high
#' depth where the band is tight, = genuine mixing.
#'
#' @param sp A `symbayes` object with the model fitted.
#' @param model One of `"lda"` or `"hdp"`.
#' @param comp_a,comp_b Component labels (e.g. `"Topic_2"`, `"Topic_9"`).
#' @param min_frac Presence threshold (default 0.05).
#' @param min_pair_reads Minimum pair reads to include a sample (default 100).
#' @param top_n,owner_frac Defining-sequence assignment controls.
#' @return A \pkg{ggplot2} object.
#' @examples
#' \dontrun{ plot_overdispersion(sp, "lda", "Topic_2", "Topic_9") }
#' @export
plot_overdispersion <- function(sp, model = c("lda", "hdp"),
                                comp_a, comp_b, min_frac = 0.05,
                                min_pair_reads = 100,
                                top_n = 15, owner_frac = 0.5) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta
  cm <- sp$count_mat[rownames(theta), , drop = FALSE]

  owner <- .defining_seqs(m$beta, top_n = top_n, owner_frac = owner_frac)
  seqs_by_comp <- split(names(owner), owner)
  comp_reads_v <- function(g) {
    s <- intersect(seqs_by_comp[[g]], colnames(cm))
    if (length(s) == 0) rep(0, nrow(cm)) else rowSums(cm[, s, drop = FALSE])
  }
  ka_all <- comp_reads_v(comp_a); kb_all <- comp_reads_v(comp_b)
  names(ka_all) <- names(kb_all) <- rownames(cm)

  both <- theta[, comp_a] >= min_frac & theta[, comp_b] >= min_frac
  ka <- ka_all[both]; kb <- kb_all[both]
  n_i <- ka + kb
  keep <- n_i >= min_pair_reads
  ka <- ka[keep]; n_i <- n_i[keep]
  share <- ka / n_i
  p_hat <- sum(ka) / sum(n_i)

  df <- data.frame(depth = n_i, share = share)
  env <- data.frame(depth = seq(min(n_i), max(n_i), length.out = 100))
  env$se <- sqrt(p_hat * (1 - p_hat) / env$depth)
  env$lo <- pmax(0, p_hat - 2 * env$se)
  env$hi <- pmin(1, p_hat + 2 * env$se)

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = env, ggplot2::aes(
      x = .data$depth, ymin = .data$lo, ymax = .data$hi),
      fill = "grey80", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = p_hat, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::geom_point(data = df, ggplot2::aes(
      x = .data$depth, y = .data$share), size = 2, alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = sprintf("%s vs %s: ratio vs depth (real reads)", comp_a, comp_b),
      subtitle = "Grey = fixed-ratio (intragenomic) null +/-2 SE. Spread beyond = mixing.",
      x = "Reads in defining-sequence pair (log scale)",
      y = sprintf("%s share of pair", comp_a)) +
    ggplot2::theme_minimal(base_size = 11)
}


#' Combined pairwise topic comparison: redundancy, mixing, and ratio stability
#'
#' One row per co-occurring component pair, synthesising the three diagnostics
#' needed to decide what a pair *is*:
#'
#' \describe{
#'   \item{`cosine`}{Exemplar (beta) cosine similarity. High (> `redundant_thresh`)
#'     = redundant topics (the same signature split in two).}
#'   \item{`n_cooccur`, `mean_ratio`, `ratio_sd`}{How often the pair co-occurs
#'     and how stable their relative proportion is across those samples.}
#'   \item{`rho`}{Beta-binomial overdispersion from [test_overdispersion] on real
#'     reads. Low (~0) = fixed ratio (intragenomic); high = variable (mixing).}
#'   \item{`verdict`}{A synthesised call (see below).}
#' }
#'
#' The verdict combines all three:
#' \itemize{
#'   \item `"redundant"` -- cosine >= `redundant_thresh` (merge candidates).
#'   \item `"intragenomic"` -- co-occur often at a fixed ratio (rho <=
#'     `rho_intra`): likely one symbiont's intragenomic variants split across
#'     topics.
#'   \item `"mixing"` -- co-occur with a variable ratio (rho >= `rho_mixing`):
#'     genuine community admixture.
#'   \item `"ambiguous"` -- co-occur but rho between the thresholds.
#'   \item `"distinct"` -- rarely co-occur (< `min_samples`): separable types.
#' }
#'
#' Compositional evidence alone cannot prove intragenomic vs a fixed-proportion
#' obligate association; a single-copy marker (psbA) is the definitive arbiter.
#'
#' @param sp A `symbayes` object with the model fitted.
#' @param model One of `"lda"` or `"hdp"`.
#' @param redundant_thresh Cosine at/above which a pair is `"redundant"`
#'   (default 0.85).
#' @param rho_intra Rho at/below which a co-occurring pair is `"intragenomic"`
#'   (default 0.03).
#' @param rho_mixing Rho at/above which a co-occurring pair is `"mixing"`
#'   (default 0.05).
#' @param min_frac,min_samples,min_pair_reads Passed to the overdispersion test.
#' @return Invisibly, a data frame of pairs with all diagnostics and `verdict`,
#'   ordered by verdict priority then strength.
#' @examples
#' \dontrun{
#' ct <- compare_topics(sp, model = "lda")
#' subset(ct, verdict == "intragenomic")
#' }
#' @export
compare_topics <- function(sp, model = c("lda", "hdp"),
                           redundant_thresh = 0.85,
                           rho_intra = 0.03, rho_mixing = 0.05,
                           min_frac = 0.05, min_samples = 5,
                           min_pair_reads = 100) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta
  grp <- colnames(theta)
  K <- ncol(theta)

  # Cosine similarity among exemplars
  sim <- .cosine_matrix(m$beta)

  # Overdispersion for co-occurring pairs (real reads)
  od <- suppressMessages(
    test_overdispersion(sp, model = model, min_frac = min_frac,
                        min_samples = min_samples,
                        min_pair_reads = min_pair_reads,
                        rho_mixing = rho_mixing))
  od_key <- if (nrow(od) > 0) paste(od$comp_a, od$comp_b) else character(0)

  rows <- list()
  for (a in seq_len(K - 1)) for (b in (a + 1):K) {
    ga <- grp[a]; gb <- grp[b]
    both <- theta[, a] >= min_frac & theta[, b] >= min_frac
    n_co <- sum(both)

    # Ratio stats across co-occurring samples (theta-level, for description)
    if (n_co >= 1) {
      ra <- theta[both, a] / (theta[both, a] + theta[both, b])
      mean_ratio <- mean(ra); ratio_sd <- stats::sd(ra)
    } else { mean_ratio <- NA; ratio_sd <- NA }

    key <- paste(ga, gb)
    rho <- if (key %in% od_key) od$rho[match(key, od_key)] else NA_real_

    verdict <-
      if (sim[ga, gb] >= redundant_thresh) "redundant"
    else if (n_co < min_samples) "distinct"
    else if (!is.na(rho) && rho <= rho_intra) "intragenomic"
    else if (!is.na(rho) && rho >= rho_mixing) "mixing"
    else "ambiguous"

    rows[[length(rows) + 1]] <- data.frame(
      comp_a = ga, comp_b = gb,
      cosine = round(sim[ga, gb], 3),
      n_cooccur = n_co,
      mean_ratio = round(mean_ratio, 3),
      ratio_sd = round(ratio_sd, 3),
      rho = if (is.na(rho)) NA else round(rho, 4),
      verdict = verdict, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)

  # Order: redundant, intragenomic, mixing, ambiguous, distinct
  pri <- c(redundant = 1, intragenomic = 2, mixing = 3,
           ambiguous = 4, distinct = 5)
  out <- out[order(pri[out$verdict], -out$n_cooccur), ]
  rownames(out) <- NULL

  tab <- table(factor(out$verdict, levels = names(pri)))
  cat(sprintf("Topic pair comparison (%s): %d pairs\n", toupper(model), nrow(out)))
  for (v in names(pri)) cat(sprintf("  %-13s %d\n", v, tab[v]))
  cat("  (psbA remains the definitive arbiter for intragenomic vs mixing)\n")
  invisible(out)
}


#' Screen for intragenomic over-splitting and suggest a biological grouping
#'
#' Second-stage screen built on [compare_topics]. Identifies fixed-ratio
#' ("intragenomic") and redundant pairs, forms their transitive closure into
#' groups (if 2~9 and 9~5 are both fixed-ratio, then {2,5,9} is one group), and
#' reports how many biological units the topics reduce to -- as a **suggestion
#' only**. It does NOT modify `theta`: the intragenomic call is provisional
#' (psbA is the arbiter), and merging on a compositional threshold would
#' reintroduce exactly the threshold-dependent clustering the package avoids.
#'
#' For each suggested group it also flags **deviant samples** -- those whose
#' within-group ratio departs from the group's fixed ratio by more than
#' `deviant_z` robust SDs. These are the samples that break the intragenomic
#' pattern: either genuine mixing exceptions or low-depth noise, surfaced for
#' inspection rather than absorbed into the group verdict.
#'
#' To actually merge a suggested grouping into a coarser membership, call
#' [merge_topics] explicitly with the grouping you choose.
#'
#' @param sp A `symbayes` object with the model fitted.
#' @param model One of `"lda"` or `"hdp"`.
#' @param deviant_z Robust z-score above which a sample is flagged as deviating
#'   from its group's fixed ratio (default 3).
#' @param ... Passed to [compare_topics] (thresholds).
#' @return Invisibly, a list: `pairs` (the compare_topics table), `groups`
#'   (data frame mapping each component to a suggested group id), `n_units`
#'   (suggested biological unit count), and `deviants` (samples breaking their
#'   group's fixed ratio).
#' @examples
#' \dontrun{
#' scr <- screen_intragenomic(sp, model = "lda")
#' scr$n_units
#' scr$groups
#' scr$deviants
#' }
#' @export
screen_intragenomic <- function(sp, model = c("lda", "hdp"),
                                deviant_z = 3, ...) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta
  grp <- colnames(theta)

  pairs <- compare_topics(sp, model = model, ...)

  # Union-find over components linked by intragenomic/redundant verdicts.
  # Operate in integer index space (1..K); idx maps a component label to index.
  idx <- setNames(seq_along(grp), grp)
  parent <- seq_along(grp)
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  union <- function(a, b) { parent[find(a)] <<- find(b) }

  link <- pairs[pairs$verdict %in% c("intragenomic", "redundant"), ]
  if (nrow(link) > 0)
    for (i in seq_len(nrow(link)))
      union(idx[[link$comp_a[i]]], idx[[link$comp_b[i]]])

  root <- vapply(seq_along(grp), find, integer(1))
  group_id <- as.integer(factor(root))
  groups <- data.frame(component = grp, group = group_id,
                       stringsAsFactors = FALSE)
  n_units <- length(unique(group_id))

  # Deviant samples: within each multi-member group, flag samples whose
  # component split departs from the group's pooled fixed ratio.
  deviants <- list()
  for (g in unique(group_id)) {
    members <- grp[group_id == g]
    if (length(members) < 2) next
    sub <- theta[, members, drop = FALSE]
    present <- rowSums(sub >= 0.05) >= 2
    if (sum(present) < 3) next
    subp <- sub[present, , drop = FALSE]
    # Ratio of first member to group total, robust centre/spread
    r <- subp[, 1] / rowSums(subp)
    med <- stats::median(r)
    mad <- stats::mad(r)
    if (mad == 0) next
    z <- abs(r - med) / mad
    dev <- names(which(z > deviant_z))
    if (length(dev) > 0)
      deviants[[as.character(g)]] <- data.frame(
        group = g, sample_uid = dev,
        ratio = round(r[dev], 3), group_median = round(med, 3),
        robust_z = round(z[dev], 2), stringsAsFactors = FALSE)
  }
  deviants <- if (length(deviants)) do.call(rbind, deviants) else NULL

  multi <- table(group_id); multi <- sum(multi > 1)
  cat(sprintf("\nIntragenomic screen (%s):\n", toupper(model)))
  cat(sprintf("  %d topics -> %d suggested biological units (%d multi-topic groups)\n",
              length(grp), n_units, multi))
  if (!is.null(deviants))
    cat(sprintf("  %d sample(s) deviate from their group's fixed ratio (flagged)\n",
                nrow(deviants)))
  cat("  Suggestion only; theta unchanged. Use merge_topics() to apply a grouping.\n")

  invisible(list(pairs = pairs, groups = groups, n_units = n_units,
                 deviants = deviants))
}


#' Merge topics into a coarser membership (explicit, user-directed)
#'
#' Applies a grouping (e.g. from [screen_intragenomic]) by summing the theta
#' columns of components in the same group, producing a coarser membership
#' matrix. This is deliberately a separate, explicit call -- the package never
#' merges automatically, because merging on a compositional threshold is a
#' judgement the user must own (and ideally confirm with psbA).
#'
#' @param sp A `symbayes` object with the model fitted.
#' @param model One of `"lda"` or `"hdp"`.
#' @param groups Either a data frame with `component` and `group` columns (as
#'   returned by `screen_intragenomic()$groups`), or a named integer vector
#'   mapping component label -> group id.
#' @return A matrix (samples x merged groups) of summed memberships. Column
#'   names list the merged components. Returned, not written into `sp`.
#' @examples
#' \dontrun{
#' scr <- screen_intragenomic(sp, model = "lda")
#' merged_theta <- merge_topics(sp, model = "lda", groups = scr$groups)
#' }
#' @export
merge_topics <- function(sp, model = c("lda", "hdp"), groups) {
  model <- match.arg(model)
  m <- .get_membership(sp, model)
  theta <- m$theta
  grp <- colnames(theta)

  if (is.data.frame(groups))
    gmap <- setNames(groups$group, groups$component)
  else gmap <- groups
  gmap <- gmap[grp]                       # align to theta column order

  ug <- sort(unique(gmap))
  merged <- vapply(ug, function(g) {
    cols <- grp[gmap == g]
    rowSums(theta[, cols, drop = FALSE])
  }, numeric(nrow(theta)))
  # Name each merged column by its constituents
  colnames(merged) <- vapply(ug, function(g)
    paste(grp[gmap == g], collapse = "+"), character(1))
  rownames(merged) <- rownames(theta)
  merged
}

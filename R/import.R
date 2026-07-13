# ==============================================================================
# import.R — read SymPortal output into a symbayes object
# ==============================================================================

#' Import SymPortal output into a `symbayes` object
#'
#' Reads a SymPortal post-MED sequence count table (the only required input)
#' and, optionally, ITS2 type profiles, sample metadata, and the SymPortal
#' colour dictionary. Returns a `symbayes` object that flows through the rest
#' of the pipeline.
#'
#' @param seqs_abund Path to `*.seqs.absolute.abund_only.txt` (samples x
#'   post-MED sequence variants). Required.
#' @param profs_abund Path to `*.profiles.absolute.abund_only.txt`
#'   (samples x ITS2 type profiles). Optional; enables profile comparison.
#' @param profs_meta Path to `*.profiles.meta_only.txt` (profile UID, name,
#'   clade). Optional.
#' @param seqs_meta Path to `*.seqs.absolute.meta_only.txt` (per-sample QC).
#'   Optional.
#' @param sample_sheet Path to a SymPortal datasheet `.xlsx`. Optional; needs
#'   the \pkg{readxl} package.
#' @param colour_dict Path to `color_dict_post_med.json`. Optional; used for
#'   SymPortal-consistent sequence colours in bar plots.
#' @param prof_colour_dict Path to `prof_color_dict.json` (SymPortal profile
#'   colours, keyed by profile UID). Optional; used for SymPortal-consistent
#'   profile colours in the comparison plots.
#' @param sheet_skip Header rows to skip in the datasheet (default 1).
#'
#' @return A `symbayes` object: a list with `count_mat`, `prof_mat`,
#'   `prof_meta`, `metadata`, `col_dict`, `seq_clades`, and `filtered = FALSE`.
#'
#' @examples
#' \dontrun{
#' sp <- import(
#'   seqs_abund  = "run.seqs.absolute.abund_only.txt",
#'   profs_abund = "run.profiles.absolute.abund_only.txt",
#'   profs_meta  = "run.profiles.meta_only.txt",
#'   colour_dict = "color_dict_post_med.json"
#' )
#' sp
#' }
#' @export
import <- function(seqs_abund,
                   profs_abund  = NULL,
                   profs_meta   = NULL,
                   seqs_meta    = NULL,
                   sample_sheet = NULL,
                   colour_dict  = NULL,
                   prof_colour_dict = NULL,
                   sheet_skip   = 1) {

  message("Importing SymPortal output ...")

  # --- Sequence count matrix (required) ---------------------------------------
  seqs_raw  <- utils::read.delim(seqs_abund, row.names = 1, check.names = FALSE)
  count_mat <- as.matrix(seqs_raw)
  storage.mode(count_mat) <- "integer"
  message(sprintf("  sequences: %d samples x %d post-MED variants",
                  nrow(count_mat), ncol(count_mat)))

  # --- Profiles (optional) ----------------------------------------------------
  prof_mat <- NULL
  prof_meta_df <- NULL
  if (!is.null(profs_abund) && file.exists(profs_abund)) {
    prof_raw <- utils::read.delim(profs_abund, row.names = 1, check.names = FALSE)
    prof_mat <- as.matrix(prof_raw)
    storage.mode(prof_mat) <- "integer"
    message(sprintf("  profiles:  %d samples x %d ITS2 type profiles",
                    nrow(prof_mat), ncol(prof_mat)))
  }
  if (!is.null(profs_meta) && file.exists(profs_meta)) {
    prof_meta_df <- utils::read.delim(profs_meta, check.names = FALSE)
  }

  # --- Sample metadata --------------------------------------------------------
  metadata <- data.frame(sample_uid = rownames(count_mat),
                         stringsAsFactors = FALSE)

  if (!is.null(seqs_meta) && file.exists(seqs_meta)) {
    sm <- utils::read.delim(seqs_meta, check.names = FALSE)
    sm$sample_uid <- as.character(sm$sample_uid)
    keep_cols <- intersect(
      c("sample_uid", "sample_name", "raw_contigs",
        "post_qc_absolute_seqs", "post_qc_unique_seqs",
        "post_taxa_id_absolute_symbiodiniaceae_seqs"),
      colnames(sm))
    metadata <- dplyr::left_join(metadata, sm[, keep_cols], by = "sample_uid")
  }

  if (!is.null(sample_sheet) && file.exists(sample_sheet)) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("Install 'readxl' to read .xlsx datasheets", call. = FALSE)
    ds <- readxl::read_excel(sample_sheet, skip = sheet_skip)
    if ("sample_name" %in% colnames(metadata) &&
        "sample_name" %in% colnames(ds)) {
      ds$sample_name <- as.character(ds$sample_name)
      metadata$sample_name <- as.character(metadata$sample_name)
      metadata <- dplyr::left_join(metadata, ds, by = "sample_name",
                                   suffix = c("", ".ds"))
    }
  }
  rownames(metadata) <- metadata$sample_uid

  # --- Colour dictionaries ----------------------------------------------------
  col_dict <- NULL
  if (!is.null(colour_dict) && file.exists(colour_dict)) {
    col_dict <- jsonlite::fromJSON(colour_dict)
  }
  # SymPortal profile colours (prof_color_dict.json), keyed by profile UID.
  prof_col_dict <- NULL
  if (!is.null(prof_colour_dict) && file.exists(prof_colour_dict)) {
    prof_col_dict <- unlist(jsonlite::fromJSON(prof_colour_dict))
  }

  sp <- structure(
    list(
      count_mat  = count_mat,
      prof_mat   = prof_mat,
      prof_meta  = prof_meta_df,
      metadata   = metadata,
      col_dict   = col_dict,
      prof_col_dict = prof_col_dict,
      seq_clades = .infer_clade(colnames(count_mat)),
      filtered   = FALSE
    ),
    class = "symbayes"
  )
  message("  done.")
  sp
}


#' Filter a `symbayes` object by prevalence and sequencing depth
#'
#' Soft-filters the count matrix: removes low-prevalence sequences and
#' low-depth samples, and re-orders sequences by abundance. Unlike SymPortal's
#' hard minimum-sequence cutoff, `min_reads` is a floor applied once, not a
#' rule that reshapes profiles.
#'
#' @param sp A `symbayes` object from [import].
#' @param min_reads Minimum total reads per sample (default 100).
#' @param min_prev Minimum number of samples a sequence must appear in
#'   (default 2).
#' @return The filtered `symbayes` object (`filtered = TRUE`).
#' @examples
#' \dontrun{
#' sp <- filter_samples(sp, min_reads = 100, min_prev = 2)
#' }
#' @export
filter_samples <- function(sp, min_reads = 100, min_prev = 2) {
  .check_sp(sp)
  mat <- sp$count_mat

  prev <- colSums(mat > 0)
  mat  <- mat[, prev >= min_prev, drop = FALSE]

  depth <- rowSums(mat)
  keep  <- depth >= max(min_reads, 1)   # never keep a zero-read sample
  mat   <- mat[keep, , drop = FALSE]
  mat   <- mat[, colSums(mat) > 0, drop = FALSE]
  mat   <- mat[, order(colSums(mat), decreasing = TRUE), drop = FALSE]

  # Final safety: drop any sample left with zero reads after the column re-cull
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]

  sp$metadata   <- sp$metadata[sp$metadata$sample_uid %in% rownames(mat), ]
  sp$count_mat  <- mat
  sp$seq_clades <- sp$seq_clades[names(sp$seq_clades) %in% colnames(mat)]
  sp$filtered   <- TRUE

  depth <- rowSums(mat)
  message(sprintf("Filtered: %d samples x %d sequences (depth %.0f-%.0f, median %.0f)",
                  nrow(mat), ncol(mat), min(depth), max(depth), median(depth)))
  sp
}


#' Stable cross-run signature for a SymPortal profile
#'
#' SymPortal profile UIDs are run-specific, but the defining-sequence
#' *accessions* are database-level identifiers stable across runs. The leading
#' accessions (the majority / co-dominant sequences) are the intrinsic identity;
#' the trailing minor DIVs are run-context-dependent (a minor sequence is
#' appended only if it co-occurred often enough in that run's samples). This
#' builds a signature from the first `n_lead` accessions, so the same biological
#' profile matches across runs even when its minor-DIV tail differs.
#'
#' @param accession_str The `Sequence accession / SymPortal UID` string, e.g.
#'   `"62825/742546-32-1147953"`.
#' @param n_lead Number of leading accessions to use (default 3).
#' @return A signature string (leading accessions, sorted, joined by `_`).
#' @keywords internal
.profile_signature <- function(accession_str, n_lead = 3) {
  vapply(accession_str, function(s) {
    if (is.na(s) || s == "") return(NA_character_)
    toks <- unlist(strsplit(s, "[-/]"))        # split on both separators
    toks <- toks[toks != ""]
    lead <- utils::head(toks, n_lead)
    paste(sort(lead), collapse = "_")          # sort: order-independent identity
  }, character(1), USE.NAMES = FALSE)
}


#' Merge multiple SymPortal runs into one `symbayes` object
#'
#' Combines several imported runs so models can be fit on a larger dataset.
#' Sequences are aligned by name and zero-filled where absent; samples are
#' stacked. **Profiles are retained and matched across runs** by their stable
#' accession signature (see [.profile_signature]) rather than by run-specific
#' UID -- so the same biological profile is recognised as one entity even when
#' its run-dependent minor-DIV tail differs.
#'
#' Why signature-matching works, and its limits: a SymPortal profile's identity
#' lives in its majority / co-dominant sequences, whose database accessions are
#' stable across runs. Its trailing minor DIVs are contextual -- appended only
#' if they co-occurred enough in that run -- so the full profile *name* can
#' differ between runs for the same community. Matching on the leading
#' accessions captures the stable identity while tolerating the variable tail.
#' The tail variability is itself a manifestation of SymPortal's core issue
#' (co-occurrence-based, context-dependent profile definition); signature
#' matching sidesteps it.
#'
#' @param ... Two or more `symbayes` objects (from [import]), or a single list.
#' @param run_ids Optional run names (provenance and collision suffixes);
#'   defaults to `run1`, `run2`, ...
#' @param on_collision Duplicate sample names across runs: `"suffix"`
#'   (append `.<run_id>`, default), `"error"`, or `"keep_first"`.
#' @param numeric_warn_frac Warn if a run has more than this fraction of its
#'   abundance in numeric-only sequences that may not align (default 0.2).
#' @param n_lead Leading accessions used to build the profile signature
#'   (default 3). Higher = stricter matching (more distinct profiles); lower =
#'   looser (more merging across runs).
#' @param merge_profiles If `TRUE` (default) retain and cross-match profiles by
#'   signature; if `FALSE` drop profiles entirely (`prof_mat = NULL`).
#' @return A merged `symbayes` object with a `run` column in metadata,
#'   `filtered = FALSE`, and (if `merge_profiles`) a `prof_mat`/`prof_meta`
#'   whose columns are signature-matched across runs. `prof_meta` gains a
#'   `signature` column and a `runs` column listing which runs each profile
#'   appeared in.
#' @examples
#' \dontrun{
#' spA <- import("runA.seqs.absolute.abund_only.txt",
#'               profs_abund = "runA.profiles.absolute.abund_only.txt",
#'               profs_meta  = "runA.profiles.meta_only.txt")
#' spB <- import("runB.seqs.absolute.abund_only.txt",
#'               profs_abund = "runB.profiles.absolute.abund_only.txt",
#'               profs_meta  = "runB.profiles.meta_only.txt")
#' sp  <- combine_symportal(spA, spB, run_ids = c("2021", "2023"))
#' sp$prof_meta[, c("ITS2 type profile", "signature", "runs")]
#' }
#' @export
combine_symportal <- function(..., run_ids = NULL,
                              on_collision = c("suffix", "error", "keep_first"),
                              numeric_warn_frac = 0.2, n_lead = 3,
                              merge_profiles = TRUE) {
  on_collision <- match.arg(on_collision)

  objs <- list(...)
  if (length(objs) == 1 && is.list(objs[[1]]) &&
      !inherits(objs[[1]], "symbayes")) objs <- objs[[1]]
  if (length(objs) < 2)
    stop("Provide at least two symbayes objects to combine.", call. = FALSE)
  for (o in objs) if (!inherits(o, "symbayes"))
    stop("All inputs must be symbayes objects from import().", call. = FALSE)

  if (is.null(run_ids)) run_ids <- paste0("run", seq_along(objs))
  if (length(run_ids) != length(objs))
    stop("run_ids must have one entry per object.", call. = FALSE)

  # --- Diagnose numeric-only (non-reference) sequence load per run ------------
  is_numeric_seq <- function(nm) grepl("^[0-9]", nm)
  for (i in seq_along(objs)) {
    m <- objs[[i]]$count_mat
    num_frac <- sum(m[, is_numeric_seq(colnames(m)), drop = FALSE]) / sum(m)
    if (num_frac > numeric_warn_frac)
      warning(sprintf(paste0(
        "Run '%s': %.0f%% of abundance is in numeric-only sequences that may ",
        "not align across runs. Cross-run merging of these is approximate."),
        run_ids[i], 100 * num_frac), call. = FALSE)
  }

  # --- Union of sequences, aligned by name ------------------------------------
  all_seqs <- Reduce(union, lapply(objs, function(o) colnames(o$count_mat)))
  aligned <- lapply(seq_along(objs), function(i) {
    m <- objs[[i]]$count_mat
    out <- matrix(0L, nrow(m), length(all_seqs),
                  dimnames = list(rownames(m), all_seqs))
    out[, colnames(m)] <- m
    out
  })

  # --- Resolve sample-name collisions -----------------------------------------
  all_names <- unlist(lapply(aligned, rownames))
  dups <- unique(all_names[duplicated(all_names)])
  suffixed <- rep(FALSE, length(objs))
  if (length(dups) > 0) {
    if (on_collision == "error") {
      stop(sprintf("Duplicate sample names across runs: %s%s",
                   paste(utils::head(dups, 5), collapse = ", "),
                   if (length(dups) > 5) ", ..." else ""), call. = FALSE)
    } else if (on_collision == "suffix") {
      for (i in seq_along(aligned)) {
        rn <- rownames(aligned[[i]]); hit <- rn %in% dups
        rn[hit] <- paste0(rn[hit], ".", run_ids[i])
        rownames(aligned[[i]]) <- rn
        if (any(hit)) suffixed[i] <- TRUE
      }
    } else { # keep_first
      seen <- character(0)
      for (i in seq_along(aligned)) {
        rn <- rownames(aligned[[i]]); drop <- rn %in% seen
        aligned[[i]] <- aligned[[i]][!drop, , drop = FALSE]
        seen <- c(seen, rownames(aligned[[i]]))
      }
    }
  }
  count_mat <- do.call(rbind, aligned)

  # --- Provenance + shared metadata -------------------------------------------
  run_vec <- unlist(lapply(seq_along(aligned), function(i)
    rep(run_ids[i], nrow(aligned[[i]]))))
  metadata <- data.frame(sample_uid = rownames(count_mat),
                         run = run_vec, stringsAsFactors = FALSE)
  shared_cols <- Reduce(intersect, lapply(objs, function(o)
    setdiff(colnames(o$metadata), "sample_uid")))
  if (length(shared_cols) > 0) {
    md_all <- do.call(rbind, lapply(objs, function(o)
      o$metadata[, c("sample_uid", shared_cols), drop = FALSE]))
    metadata <- dplyr::left_join(metadata, md_all,
                                 by = "sample_uid", multiple = "first")
  }
  rownames(metadata) <- metadata$sample_uid

  # --- Profiles: retain and cross-match by accession signature ----------------
  prof_mat <- NULL; prof_meta <- NULL
  have_profs <- all(vapply(objs, function(o)
    !is.null(o$prof_mat) && !is.null(o$prof_meta), logical(1)))

  if (merge_profiles && have_profs) {
    acc_col <- "Sequence accession / SymPortal UID"
    uid_col <- "ITS2 type profile UID"
    name_col <- "ITS2 type profile"

    # Build a per-run tidy table: sample_uid x signature abundance
    prof_long <- list()
    sig_meta  <- list()
    for (i in seq_along(objs)) {
      pm <- objs[[i]]$prof_meta
      pmat <- objs[[i]]$prof_mat
      if (!acc_col %in% colnames(pm)) {
        warning(sprintf("Run '%s' lacks accession column; profiles skipped.",
                        run_ids[i]), call. = FALSE)
        next
      }
      sig <- .profile_signature(pm[[acc_col]], n_lead = n_lead)
      names(sig) <- as.character(pm[[uid_col]])

      # Signature per profile-matrix column (columns are profile UIDs)
      col_sig <- sig[colnames(pmat)]
      # Row names may have been suffixed for collisions
      rn <- rownames(pmat)
      if (on_collision == "suffix" && suffixed[i]) {
        hit <- rn %in% dups
        rn[hit] <- paste0(rn[hit], ".", run_ids[i])
        rownames(pmat) <- rn
      }
      # Collapse columns sharing a signature within this run (sum abundance)
      keep_cols <- !is.na(col_sig)
      pmat <- pmat[, keep_cols, drop = FALSE]
      col_sig <- col_sig[keep_cols]
      agg <- t(rowsum(t(pmat), group = col_sig))   # signatures x samples -> back
      prof_long[[i]] <- agg                         # samples x signatures

      # Keep a representative metadata row per signature (first occurrence)
      rep_meta <- pm[!duplicated(sig) & !is.na(sig),
                     c(uid_col, name_col, acc_col, "Clade",
                       "Majority ITS2 sequence"), drop = FALSE]
      rep_meta$signature <- sig[!duplicated(sig) & !is.na(sig)]
      rep_meta$run <- run_ids[i]
      sig_meta[[i]] <- rep_meta
    }

    if (length(prof_long) > 0) {
      all_sigs <- Reduce(union, lapply(prof_long, colnames))
      prof_mat <- matrix(0L, nrow(count_mat), length(all_sigs),
                         dimnames = list(rownames(count_mat), all_sigs))
      for (agg in prof_long) {
        if (is.null(agg)) next
        prof_mat[rownames(agg), colnames(agg)] <- agg
      }

      meta_all <- do.call(rbind, sig_meta)
      # One row per signature; record which runs it appeared in
      runs_by_sig <- tapply(meta_all$run, meta_all$signature,
                            function(x) paste(sort(unique(x)), collapse = ","))
      prof_meta <- meta_all[!duplicated(meta_all$signature), , drop = FALSE]
      prof_meta$runs <- runs_by_sig[prof_meta$signature]
      rownames(prof_meta) <- NULL

      n_shared <- sum(grepl(",", prof_meta$runs))
      message(sprintf(
        "  Profiles: %d unique signatures across runs (%d shared by >1 run).",
        nrow(prof_meta), n_shared))
    }
  } else if (merge_profiles && !have_profs) {
    warning("Not all runs have profiles; merged object will have prof_mat = NULL.",
            call. = FALSE)
  }

  col_dict <- NULL
  for (o in objs) if (!is.null(o$col_dict)) { col_dict <- o$col_dict; break }

  message(sprintf("Combined %d runs -> %d samples x %d sequences.",
                  length(objs), nrow(count_mat), ncol(count_mat)))

  structure(
    list(
      count_mat  = count_mat,
      prof_mat   = prof_mat,
      prof_meta  = prof_meta,
      metadata   = metadata,
      col_dict   = col_dict,
      seq_clades = .infer_clade(colnames(count_mat)),
      filtered   = FALSE,
      merged     = TRUE,
      run_ids    = run_ids
    ),
    class = "symbayes")
}

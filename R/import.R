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

  # --- Colour dictionary ------------------------------------------------------
  col_dict <- NULL
  if (!is.null(colour_dict) && file.exists(colour_dict)) {
    col_dict <- jsonlite::fromJSON(colour_dict)
  }

  sp <- structure(
    list(
      count_mat  = count_mat,
      prof_mat   = prof_mat,
      prof_meta  = prof_meta_df,
      metadata   = metadata,
      col_dict   = col_dict,
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
  keep  <- depth >= min_reads
  mat   <- mat[keep, , drop = FALSE]
  mat   <- mat[, colSums(mat) > 0, drop = FALSE]
  mat   <- mat[, order(colSums(mat), decreasing = TRUE), drop = FALSE]

  sp$metadata   <- sp$metadata[sp$metadata$sample_uid %in% rownames(mat), ]
  sp$count_mat  <- mat
  sp$seq_clades <- sp$seq_clades[names(sp$seq_clades) %in% colnames(mat)]
  sp$filtered   <- TRUE

  depth <- rowSums(mat)
  message(sprintf("Filtered: %d samples x %d sequences (depth %.0f-%.0f, median %.0f)",
                  nrow(mat), ncol(mat), min(depth), max(depth), median(depth)))
  sp
}

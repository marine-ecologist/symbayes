# Tests for the membership schema and core object plumbing.
# These use a small synthetic count matrix so they run without SymPortal files.

make_toy <- function(n = 30, s = 12, seed = 1) {
  set.seed(seed)
  # Two latent groups with distinct dominant sequences
  mat <- matrix(0L, n, s,
                dimnames = list(paste0("samp", seq_len(n)),
                                c(paste0("C", seq_len(s / 2)),
                                  paste0("D", seq_len(s / 2)))))
  for (i in seq_len(n)) {
    if (i <= n / 2) {
      probs <- c(rep(5, s / 2), rep(1, s / 2))
    } else {
      probs <- c(rep(1, s / 2), rep(5, s / 2))
    }
    mat[i, ] <- stats::rmultinom(1, 2000, probs / sum(probs))[, 1]
  }
  structure(
    list(count_mat = mat, prof_mat = NULL, prof_meta = NULL,
         metadata = data.frame(sample_uid = rownames(mat),
                               row.names = rownames(mat)),
         col_dict = NULL,
         seq_clades = symbayes:::.infer_clade(colnames(mat)),
         filtered = TRUE),
    class = "symbayes")
}

test_that("import produces a valid symbayes object structure", {
  sp <- make_toy()
  expect_s3_class(sp, "symbayes")
  expect_true(is.matrix(sp$count_mat))
  expect_equal(nrow(sp$count_mat), 30)
})

test_that(".infer_clade reads leading and trailing clade letters", {
  cl <- symbayes:::.infer_clade(c("C3", "D1", "15443_A", "B7", "unknownseq"))
  expect_equal(unname(cl["C3"]), "C")
  expect_equal(unname(cl["D1"]), "D")
  expect_equal(unname(cl["15443_A"]), "A")
  expect_equal(unname(cl["B7"]), "B")
  expect_true(is.na(cl["unknownseq"]))
})

test_that("membership constructor standardises labels", {
  theta <- matrix(runif(30), 10, 3)
  theta <- theta / rowSums(theta)
  beta  <- matrix(runif(15), 3, 5)
  m <- symbayes:::.membership(theta, beta, dominant = factor(rep(1:3, length.out = 10)),
                              k = 3, model = "lda", label_prefix = "Topic")
  expect_s3_class(m, "symbayes_membership")
  expect_equal(colnames(m$theta), c("Topic_1", "Topic_2", "Topic_3"))
  expect_equal(rownames(m$beta),  c("Topic_1", "Topic_2", "Topic_3"))
})

test_that(".get_membership errors clearly when model absent", {
  sp <- make_toy()
  expect_error(symbayes:::.get_membership(sp, "hdp"), "No HDP fit")
})

test_that("fit_dmm populates the schema", {
  skip_if_not_installed("DirichletMultinomial")
  sp <- make_toy()
  sp <- fit_dmm(sp, max_k = 3, seed = 1)
  expect_false(is.null(sp$dmm))
  expect_true(sp$dmm$k >= 1 && sp$dmm$k <= 3)
  expect_equal(nrow(sp$dmm$theta), nrow(sp$count_mat))
  expect_equal(ncol(sp$dmm$beta), ncol(sp$count_mat))
  expect_true(all(c("certainty", "model_fit", "best_fit") %in%
                    names(sp$dmm$meta)))
})

test_that("cosine matrix is symmetric with unit diagonal", {
  m <- matrix(runif(20), 4, 5)
  s <- symbayes:::.cosine_matrix(m)
  expect_equal(diag(s), rep(1, 4), tolerance = 1e-8)
  expect_equal(s, t(s), tolerance = 1e-8)
})

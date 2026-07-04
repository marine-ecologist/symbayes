#' symbayes: Threshold-Free, Mixture-Aware Symbiodiniaceae ITS2 Profiling
#'
#' `symbayes` provides a probabilistic alternative to SymPortal's deterministic
#' ITS2 type profiles. It reads SymPortal output and fits three complementary
#' mixture models, each answering a different question:
#'
#' \describe{
#'   \item{**DMM** ([fit_dmm])}{How many discrete community *types* does the
#'     dataset support? Dirichlet-Multinomial Mixture with formal model
#'     selection (Laplace/BIC). Hard assignment: each sample maps to one
#'     component.}
#'   \item{**LDA** ([fit_lda])}{What is each sample *made of*? Latent Dirichlet
#'     Allocation with redundancy-based topic-number selection. Fractional
#'     membership: a sample can be a mixture of topics.}
#'   \item{**HDP** ([fit_hdp])}{Threshold-free and mixture-aware: the number of
#'     community types is inferred from the data via a Hierarchical Dirichlet
#'     Process, and samples may be admixed.}
#' }
#'
#' # The `sp` object and the membership schema
#'
#' Every function takes and returns a single list of class `symbayes` (created
#' by [import]). Fitting a model writes a standardised *membership* object into
#' `sp$dmm`, `sp$lda`, or `sp$hdp`, each with the same schema:
#'
#' \describe{
#'   \item{`theta`}{samples x groups membership matrix (near one-hot for DMM,
#'     fractional for LDA/HDP)}
#'   \item{`beta`}{groups x sequences exemplar compositions}
#'   \item{`dominant`}{per-sample dominant group (factor)}
#'   \item{`k`}{number of groups}
#'   \item{`meta`}{model-specific extras (certainty, entropy, residual, ...)}
#' }
#'
#' Because all three models expose the same schema, the *generalist* plotting
#' functions ([plot_pcoa], [plot_barplot], [plot_exemplars],
#' [plot_contingency], [plot_membership]) work on any fitted model via a
#' `model = c("dmm", "lda", "hdp")` argument.
#'
#' # Function tiers
#'
#' \describe{
#'   \item{Entry points}{[import], [fit_dmm], [fit_lda], [fit_hdp]}
#'   \item{Generalist plots (`model=`)}{[plot_pcoa], [plot_barplot],
#'     [plot_exemplars], [plot_contingency], [plot_membership],
#'     [plot_top_seqs]}
#'   \item{DMM-specific}{[dmm_plot_selection], [dmm_loo],
#'     [dmm_subsample_cv], [dmm_simulate]}
#'   \item{LDA-specific}{[lda_plot_ksearch], [lda_plot_similarity],
#'     [lda_topic_summary]}
#'   \item{HDP-specific}{[hdp_plot_numcomp], [hdp_vs_symportal]}
#'   \item{Cross-model comparison}{[compare_models], [compare_mixing]}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom stats median setNames
#' @importFrom utils head read.delim
## usethis namespace: end
NULL

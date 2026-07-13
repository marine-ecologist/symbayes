# Package index

## Entry points

Read a SymPortal run and prepare it for modelling. Everything flows
through a single `symbayes` object.

- [`import()`](https://marine-ecologist.github.io/symbayes/reference/import.md)
  :

  Import SymPortal output into a `symbayes` object

- [`filter_samples()`](https://marine-ecologist.github.io/symbayes/reference/filter_samples.md)
  :

  Filter a `symbayes` object by prevalence and sequencing depth

## Fit a model

The three mixture models. Each fits, selects (or infers) the number of
groups, and writes a standardised membership object into the `symbayes`
object.

- [`fit_dmm()`](https://marine-ecologist.github.io/symbayes/reference/fit_dmm.md)
  : Fit a Dirichlet-Multinomial Mixture model
- [`fit_lda()`](https://marine-ecologist.github.io/symbayes/reference/fit_lda.md)
  : Fit a Latent Dirichlet Allocation topic model
- [`fit_hdp()`](https://marine-ecologist.github.io/symbayes/reference/fit_hdp.md)
  : Fit a Hierarchical Dirichlet Process

## Generalist plots

Visualisations that apply to any fitted model via a
`model = c("dmm", "lda", "hdp")` argument.

- [`plot_pcoa()`](https://marine-ecologist.github.io/symbayes/reference/plot_pcoa.md)
  : PCoA ordination coloured by a fitted model's grouping (or any
  variable)
- [`plot_barplot()`](https://marine-ecologist.github.io/symbayes/reference/plot_barplot.md)
  : Stacked sequence-composition barplot faceted by a model's groups
- [`plot_exemplars()`](https://marine-ecologist.github.io/symbayes/reference/plot_exemplars.md)
  : Exemplar composition per group (what a "pure" sample looks like)
- [`plot_membership()`](https://marine-ecologist.github.io/symbayes/reference/plot_membership.md)
  : Membership heatmap (per-sample group proportions)
- [`plot_contingency()`](https://marine-ecologist.github.io/symbayes/reference/plot_contingency.md)
  : Contingency heatmap: a model's groups vs SymPortal profiles
- [`plot_top_seqs()`](https://marine-ecologist.github.io/symbayes/reference/plot_top_seqs.md)
  : Top characterising sequences per group
- [`plot_overdispersion()`](https://marine-ecologist.github.io/symbayes/reference/plot_overdispersion.md)
  : Plot the overdispersion result for one component pair (real read
  counts)

## DMM: diagnostics, prediction, validation

Functions specific to the Dirichlet-Multinomial Mixture (hard
assignment).

- [`dmm_plot_selection()`](https://marine-ecologist.github.io/symbayes/reference/dmm_plot_selection.md)
  : Plot DMM model selection (Laplace, AIC, BIC vs K)
- [`predict_dmm()`](https://marine-ecologist.github.io/symbayes/reference/predict_dmm.md)
  : Predict DMM component membership for new samples
- [`dmm_loo()`](https://marine-ecologist.github.io/symbayes/reference/dmm_loo.md)
  : Leave-one-out cross-validation for DMM
- [`dmm_subsample_cv()`](https://marine-ecologist.github.io/symbayes/reference/dmm_subsample_cv.md)
  : Repeated subsampling cross-validation for DMM
- [`dmm_simulate()`](https://marine-ecologist.github.io/symbayes/reference/dmm_simulate.md)
  : Simulate samples and test DMM assignment (blend / perturb /
  from_alpha)

## LDA: topic diagnostics

Functions specific to Latent Dirichlet Allocation (fractional
membership).

- [`lda_plot_ksearch()`](https://marine-ecologist.github.io/symbayes/reference/lda_plot_ksearch.md)
  : Plot the LDA K-selection search
- [`lda_plot_similarity()`](https://marine-ecologist.github.io/symbayes/reference/lda_plot_similarity.md)
  : Plot the LDA topic cosine-similarity matrix
- [`lda_topic_summary()`](https://marine-ecologist.github.io/symbayes/reference/lda_topic_summary.md)
  : Text summary of LDA topics with SymPortal / host cross-reference

## HDP: threshold-free diagnostics

Functions specific to the Hierarchical Dirichlet Process, including
over-classification quantification.

- [`hdp_plot_numcomp()`](https://marine-ecologist.github.io/symbayes/reference/hdp_plot_numcomp.md)
  : Plot the HDP posterior distribution of component number
- [`hdp_vs_symportal()`](https://marine-ecologist.github.io/symbayes/reference/hdp_vs_symportal.md)
  : Quantify SymPortal over-classification relative to HDP

## Cross-model comparison

Compare group counts and within-sample mixing across models.

- [`combine_symportal()`](https://marine-ecologist.github.io/symbayes/reference/combine_symportal.md)
  :

  Merge multiple SymPortal runs into one `symbayes` object

- [`compare_models()`](https://marine-ecologist.github.io/symbayes/reference/compare_models.md)
  : Compare group counts and selection strategy across fitted models

- [`compare_mixing()`](https://marine-ecologist.github.io/symbayes/reference/compare_mixing.md)
  : Compare within-sample mixing across models and SymPortal

- [`compare_soft_hard()`](https://marine-ecologist.github.io/symbayes/reference/compare_soft_hard.md)
  : Compare a soft (fractional) model against SymPortal's profiles via
  the fuzzy Rand index

- [`compare_topics()`](https://marine-ecologist.github.io/symbayes/reference/compare_topics.md)
  : Combined pairwise topic comparison: redundancy, mixing, and ratio
  stability

- [`match_profiles()`](https://marine-ecologist.github.io/symbayes/reference/match_profiles.md)
  : Match model topics/components to SymPortal profiles by composition

- [`merge_topics()`](https://marine-ecologist.github.io/symbayes/reference/merge_topics.md)
  : Merge topics into a coarser membership (explicit, user-directed)

- [`name_topics()`](https://marine-ecologist.github.io/symbayes/reference/name_topics.md)
  : Name a fitted model's topics/components in SymPortal style

- [`plot_admixtures()`](https://marine-ecologist.github.io/symbayes/reference/plot_admixtures.md)
  : Admixture plot: mixed samples, sequence mix over boxed component
  structure

- [`plot_soft_hard()`](https://marine-ecologist.github.io/symbayes/reference/plot_soft_hard.md)
  : Visualise soft-vs-hard co-membership agreement as heatmaps

- [`plot_symportal_comparison()`](https://marine-ecologist.github.io/symbayes/reference/plot_symportal_comparison.md)
  : Four-row SymPortal-vs-model comparison over all samples

- [`set_dominance()`](https://marine-ecologist.github.io/symbayes/reference/set_dominance.md)
  : Relabel dominant assignment with a minimum-dominance threshold

- [`screen_intragenomic()`](https://marine-ecologist.github.io/symbayes/reference/screen_intragenomic.md)
  : Screen for intragenomic over-splitting and suggest a biological
  grouping

- [`test_overdispersion()`](https://marine-ecologist.github.io/symbayes/reference/test_overdispersion.md)
  : Depth-aware overdispersion test (rebuilt on real read counts)

## Package

- [`symbayes`](https://marine-ecologist.github.io/symbayes/reference/symbayes-package.md)
  [`symbayes-package`](https://marine-ecologist.github.io/symbayes/reference/symbayes-package.md)
  : symbayes: Threshold-Free, Mixture-Aware Symbiodiniaceae ITS2
  Profiling
- [`print(`*`<symbayes>`*`)`](https://marine-ecologist.github.io/symbayes/reference/print.symbayes.md)
  : Print method for symbayes objects
- [`print(`*`<symbayes_membership>`*`)`](https://marine-ecologist.github.io/symbayes/reference/print.symbayes_membership.md)
  : Print method for membership objects

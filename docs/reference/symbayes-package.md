# symbayes: Threshold-Free, Mixture-Aware Symbiodiniaceae ITS2 Profiling

`symbayes` provides a probabilistic alternative to SymPortal's
deterministic ITS2 type profiles. It reads SymPortal output and fits
three complementary mixture models, each answering a different question:

## Details

- **DMM**
  ([fit_dmm](https://marine-ecologist.github.io/symbayes/reference/fit_dmm.md)):

  How many discrete community *types* does the dataset support?
  Dirichlet-Multinomial Mixture with formal model selection
  (Laplace/BIC). Hard assignment: each sample maps to one component.

- **LDA**
  ([fit_lda](https://marine-ecologist.github.io/symbayes/reference/fit_lda.md)):

  What is each sample *made of*? Latent Dirichlet Allocation with
  redundancy-based topic-number selection. Fractional membership: a
  sample can be a mixture of topics.

- **HDP**
  ([fit_hdp](https://marine-ecologist.github.io/symbayes/reference/fit_hdp.md)):

  Threshold-free and mixture-aware: the number of community types is
  inferred from the data via a Hierarchical Dirichlet Process, and
  samples may be admixed.

## The `sp` object and the membership schema

Every function takes and returns a single list of class `symbayes`
(created by
[import](https://marine-ecologist.github.io/symbayes/reference/import.md)).
Fitting a model writes a standardised *membership* object into `sp$dmm`,
`sp$lda`, or `sp$hdp`, each with the same schema:

- `theta`:

  samples x groups membership matrix (near one-hot for DMM, fractional
  for LDA/HDP)

- `beta`:

  groups x sequences exemplar compositions

- `dominant`:

  per-sample dominant group (factor)

- `k`:

  number of groups

- `meta`:

  model-specific extras (certainty, entropy, residual, ...)

Because all three models expose the same schema, the *generalist*
plotting functions
([plot_pcoa](https://marine-ecologist.github.io/symbayes/reference/plot_pcoa.md),
[plot_barplot](https://marine-ecologist.github.io/symbayes/reference/plot_barplot.md),
[plot_exemplars](https://marine-ecologist.github.io/symbayes/reference/plot_exemplars.md),
[plot_contingency](https://marine-ecologist.github.io/symbayes/reference/plot_contingency.md),
[plot_membership](https://marine-ecologist.github.io/symbayes/reference/plot_membership.md))
work on any fitted model via a `model = c("dmm", "lda", "hdp")`
argument.

## Function tiers

- Entry points:

  [import](https://marine-ecologist.github.io/symbayes/reference/import.md),
  [fit_dmm](https://marine-ecologist.github.io/symbayes/reference/fit_dmm.md),
  [fit_lda](https://marine-ecologist.github.io/symbayes/reference/fit_lda.md),
  [fit_hdp](https://marine-ecologist.github.io/symbayes/reference/fit_hdp.md)

- Generalist plots (`model=`):

  [plot_pcoa](https://marine-ecologist.github.io/symbayes/reference/plot_pcoa.md),
  [plot_barplot](https://marine-ecologist.github.io/symbayes/reference/plot_barplot.md),
  [plot_exemplars](https://marine-ecologist.github.io/symbayes/reference/plot_exemplars.md),
  [plot_contingency](https://marine-ecologist.github.io/symbayes/reference/plot_contingency.md),
  [plot_membership](https://marine-ecologist.github.io/symbayes/reference/plot_membership.md),
  [plot_top_seqs](https://marine-ecologist.github.io/symbayes/reference/plot_top_seqs.md)

- DMM-specific:

  [dmm_plot_selection](https://marine-ecologist.github.io/symbayes/reference/dmm_plot_selection.md),
  [dmm_loo](https://marine-ecologist.github.io/symbayes/reference/dmm_loo.md),
  [dmm_subsample_cv](https://marine-ecologist.github.io/symbayes/reference/dmm_subsample_cv.md),
  [dmm_simulate](https://marine-ecologist.github.io/symbayes/reference/dmm_simulate.md)

- LDA-specific:

  [lda_plot_ksearch](https://marine-ecologist.github.io/symbayes/reference/lda_plot_ksearch.md),
  [lda_plot_similarity](https://marine-ecologist.github.io/symbayes/reference/lda_plot_similarity.md),
  [lda_topic_summary](https://marine-ecologist.github.io/symbayes/reference/lda_topic_summary.md)

- HDP-specific:

  [hdp_plot_numcomp](https://marine-ecologist.github.io/symbayes/reference/hdp_plot_numcomp.md),
  [hdp_vs_symportal](https://marine-ecologist.github.io/symbayes/reference/hdp_vs_symportal.md)

- Cross-model comparison:

  [compare_models](https://marine-ecologist.github.io/symbayes/reference/compare_models.md),
  [compare_mixing](https://marine-ecologist.github.io/symbayes/reference/compare_mixing.md)

## See also

Useful links:

- <https://github.com/yourname/symbayes>

- Report bugs at <https://github.com/yourname/symbayes/issues>

## Author

**Maintainer**: Jez Roff <you@example.org>

Authors:

- Jez Roff <you@example.org>

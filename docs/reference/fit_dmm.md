# Fit a Dirichlet-Multinomial Mixture model

Fits DMMs for `k = 1:max_k`, selects the number of components by the
Laplace approximation to the model evidence (with BIC reported
alongside), extracts posterior membership, and, if SymPortal profiles
are present, maps each sample's dominant profile and clade.

## Usage

``` r
fit_dmm(sp, max_k = 10, seed = 42)
```

## Arguments

- sp:

  A filtered `symbayes` object.

- max_k:

  Maximum number of components to evaluate (default 10).

- seed:

  Random seed (default 42).

## Value

The `symbayes` object with `sp$dmm` populated.

## Details

DMM is a *hard-assignment* model: each sample is assigned to one
component. Its posterior certainty is often near 1.0 even for mixed
samples, because it models between-sample variation, not within-sample
mixing. Use
[fit_lda](https://marine-ecologist.github.io/symbayes/reference/fit_lda.md)
or
[fit_hdp](https://marine-ecologist.github.io/symbayes/reference/fit_hdp.md)
to represent admixture. See
[dmm_loo](https://marine-ecologist.github.io/symbayes/reference/dmm_loo.md)
and
[dmm_simulate](https://marine-ecologist.github.io/symbayes/reference/dmm_simulate.md)
to probe assignment stability.

Writes a membership object to `sp$dmm` with the standard schema
(`theta`, `beta`, `dominant`, `k`, `meta`). The `meta` list contains
`certainty` (per-sample max posterior), `model_fit` (the K-selection
table), and `best_fit` (the fitted `DMN` object, used by dmm_predict
etc.).

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- fit_dmm(sp, max_k = 10)
sp$dmm$k
dmm_plot_selection(sp)
} # }
```

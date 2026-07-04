# Construct a standardised membership object

Internal constructor used by the `fit_*()` functions. All three models
emit this same structure so generalist plots can consume any of them.

## Usage

``` r
.membership(
  theta,
  beta,
  dominant,
  k,
  model,
  label_prefix = "Group",
  meta = list()
)
```

## Arguments

- theta:

  samples x groups membership matrix (rows sum to ~1)

- beta:

  groups x sequences exemplar matrix

- dominant:

  named integer/factor of per-sample dominant group

- k:

  number of groups

- model:

  one of "dmm", "lda", "hdp"

- label_prefix:

  prefix for group labels (e.g. "Comp", "Topic", "HDP")

- meta:

  named list of model-specific extras

## Value

A list of class `symbayes_membership`

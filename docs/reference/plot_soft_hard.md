# Visualise soft-vs-hard co-membership agreement as heatmaps

Renders the pairwise co-membership structure behind
[compare_soft_hard](https://marine-ecologist.github.io/symbayes/reference/compare_soft_hard.md):
the model's soft co-membership matrix, SymPortal's, and their
difference. Samples are ordered by the model's dominant group so block
structure is visible. In the difference panel, warm cells are pairs the
model co-clusters more than SymPortal (mixing SymPortal splits apart);
cool cells are the reverse (SymPortal groups pairs the model separates).

## Usage

``` r
plot_soft_hard(
  sp,
  model = c("lda", "hdp"),
  panels = c("model", "sp", "diff"),
  mixed_thresh = 0.7
)
```

## Arguments

- sp:

  A `symbayes` object with the chosen model fitted and SymPortal
  profiles present.

- model:

  One of `"lda"` or `"hdp"`.

- panels:

  Which panels to show: any of `"model"`, `"sp"`, `"diff"` (default all
  three).

- mixed_thresh:

  Passed to
  [compare_soft_hard](https://marine-ecologist.github.io/symbayes/reference/compare_soft_hard.md)
  (default 0.7).

## Value

A patchwork object (or a single ggplot if one panel).

## Examples

``` r
if (FALSE)  plot_soft_hard(sp, model = "lda")  # \dontrun{}
```

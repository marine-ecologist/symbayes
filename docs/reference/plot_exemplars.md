# Exemplar composition per group (what a "pure" sample looks like)

Stacked bar of each group's exemplar sequence distribution (`beta`),
with a clade annotation and dominant-sample count per group.

## Usage

``` r
plot_exemplars(sp, model = "dmm", min_pct = 1)
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- min_pct:

  Minimum percent in any group to show a sequence (default 1).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_exemplars(sp, model = "lda")  # \dontrun{}
```

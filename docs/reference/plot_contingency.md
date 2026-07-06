# Contingency heatmap: a model's groups vs SymPortal profiles

Tile grid of SymPortal ITS2 type profiles (rows, grouped by clade)
against a model's groups (columns); tile colour = mean assignment
confidence, number = sample count.

## Usage

``` r
plot_contingency(sp, model = "dmm")
```

## Arguments

- sp:

  A `symbayes` object with SymPortal profiles present.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_contingency(sp, model = "hdp")  # \dontrun{}
```

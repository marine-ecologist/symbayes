# Membership heatmap (per-sample group proportions)

Heatmap of `theta`: near-binary for DMM, gradient for LDA/HDP.

## Usage

``` r
plot_membership(sp, model = "dmm")
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_membership(sp, model = "lda")  # \dontrun{}
```

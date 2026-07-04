# Plot the LDA K-selection search

Four panels across the searched K: maximum pairwise topic cosine
similarity (with the redundancy threshold marked), number of redundant
pairs, percent of topics with pure samples, and perplexity. The selected
K (largest with no redundant pairs) is marked. Requires `fit_lda(sp)` to
have been run with adaptive selection (i.e. `k = NULL`).

## Usage

``` r
lda_plot_ksearch(sp)
```

## Arguments

- sp:

  A `symbayes` object with LDA fitted adaptively.

## Value

A patchwork composition of four ggplot2 panels.

## Examples

``` r
if (FALSE)  lda_plot_ksearch(sp)  # \dontrun{}
```

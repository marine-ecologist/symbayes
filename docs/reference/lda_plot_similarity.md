# Plot the LDA topic cosine-similarity matrix

Heatmap of pairwise topic cosine similarity, annotated with values.
Pairs above ~0.85 indicate redundant (over-split) topics.

## Usage

``` r
lda_plot_similarity(sp)
```

## Arguments

- sp:

  A `symbayes` object with LDA fitted.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  lda_plot_similarity(sp)  # \dontrun{}
```

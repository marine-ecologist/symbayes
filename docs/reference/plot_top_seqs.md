# Top characterising sequences per group

Faceted bar chart of the highest-weight sequences in each group's
exemplar, coloured by clade.

## Usage

``` r
plot_top_seqs(sp, model = "dmm", top_n = 8)
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- top_n:

  Sequences per group (default 8).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_top_seqs(sp, model = "dmm")  # \dontrun{}
```

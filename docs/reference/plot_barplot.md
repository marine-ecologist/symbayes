# Stacked sequence-composition barplot faceted by a model's groups

ITS2 sequence composition per sample, faceted by the dominant group of a
fitted model, using SymPortal colours where available. Optionally
annotates each bar with the model's assignment confidence (DMM certainty
/ LDA-HDP dominant proportion).

## Usage

``` r
plot_barplot(
  sp,
  model = "dmm",
  min_pct = 1,
  sample_label = "sample_name",
  show_conf = TRUE
)
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- min_pct:

  Minimum percent in any sample to show a sequence (default 1).

- sample_label:

  Metadata column for x labels (default `"sample_name"`).

- show_conf:

  Annotate bars with confidence (default `TRUE`).

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_barplot(sp, model = "hdp")  # \dontrun{}
```

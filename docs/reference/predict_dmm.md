# Predict DMM component membership for new samples

Scores new count vectors against a fitted DMM's components using the
closed-form Dirichlet-Multinomial predictive, returning posterior
membership. Columns are aligned to the fitted sequences (missing
zero-filled, extras dropped). No refitting is required.

## Usage

``` r
predict_dmm(sp, new_counts)
```

## Arguments

- sp:

  A `symbayes` object with DMM fitted.

- new_counts:

  Matrix (samples x sequences) or a single named vector.

## Value

A data frame with `assigned`, `certainty`, `total_reads`, and per
component posterior probabilities.

## Examples

``` r
if (FALSE)  predict_dmm(sp, new_counts)  # \dontrun{}
```

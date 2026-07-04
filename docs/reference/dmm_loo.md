# Leave-one-out cross-validation for DMM

For each sample: refit the DMM on the other N-1 samples and predict the
held-out sample. Reports concordance with the full-model assignment and
flags samples that switch component (candidate boundary/mixed samples).

## Usage

``` r
dmm_loo(sp, seed = 42)
```

## Arguments

- sp:

  A `symbayes` object with DMM fitted.

- seed:

  Random seed (default 42).

## Value

A data frame of per-sample LOO results, returned invisibly.

## Examples

``` r
if (FALSE)  loo <- dmm_loo(sp)  # \dontrun{}
```

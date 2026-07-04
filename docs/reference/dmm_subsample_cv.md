# Repeated subsampling cross-validation for DMM

Holds out a fraction of samples, refits, predicts held-out; repeats.
Reports overall concordance and per-sample stability.

## Usage

``` r
dmm_subsample_cv(sp, frac_test = 0.2, n_reps = 50, seed = 42)
```

## Arguments

- sp:

  A `symbayes` object with DMM fitted.

- frac_test:

  Fraction held out per replicate (default 0.2).

- n_reps:

  Replicates (default 50).

- seed:

  Random seed (default 42).

## Value

A data frame of per-sample stability, returned invisibly.

## Examples

``` r
if (FALSE)  dmm_subsample_cv(sp)  # \dontrun{}
```

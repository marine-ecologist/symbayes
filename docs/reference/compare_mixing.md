# Compare within-sample mixing across models and SymPortal

Computes per-sample "mixedness" (1 - dominant proportion) and the
effective number of groups (Hill number, q = 1) for LDA, HDP, DMM, and
SymPortal (multi-profile proportions), and reports correlations between
them. DMM mixedness is ~0 by construction (hard assignment).

## Usage

``` r
compare_mixing(sp)
```

## Arguments

- sp:

  A `symbayes` object.

## Value

A data frame of per-sample mixing metrics, returned invisibly.

## Examples

``` r
if (FALSE)  mix <- compare_mixing(sp)  # \dontrun{}
```

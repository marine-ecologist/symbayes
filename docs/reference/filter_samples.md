# Filter a `symbayes` object by prevalence and sequencing depth

Soft-filters the count matrix: removes low-prevalence sequences and
low-depth samples, and re-orders sequences by abundance. Unlike
SymPortal's hard minimum-sequence cutoff, `min_reads` is a floor applied
once, not a rule that reshapes profiles.

## Usage

``` r
filter_samples(sp, min_reads = 100, min_prev = 2)
```

## Arguments

- sp:

  A `symbayes` object from
  [import](https://marine-ecologist.github.io/symbayes/reference/import.md).

- min_reads:

  Minimum total reads per sample (default 100).

- min_prev:

  Minimum number of samples a sequence must appear in (default 2).

## Value

The filtered `symbayes` object (`filtered = TRUE`).

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- filter_samples(sp, min_reads = 100, min_prev = 2)
} # }
```

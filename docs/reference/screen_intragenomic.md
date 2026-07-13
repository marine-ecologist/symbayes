# Screen for intragenomic over-splitting and suggest a biological grouping

Second-stage screen built on
[compare_topics](https://marine-ecologist.github.io/symbayes/reference/compare_topics.md).
Identifies fixed-ratio ("intragenomic") and redundant pairs, forms their
transitive closure into groups (if 2~9 and 9~5 are both fixed-ratio,
then 2,5,9 is one group), and reports how many biological units the
topics reduce to – as a **suggestion only**. It does NOT modify `theta`:
the intragenomic call is provisional (psbA is the arbiter), and merging
on a compositional threshold would reintroduce exactly the
threshold-dependent clustering the package avoids.

## Usage

``` r
screen_intragenomic(sp, model = c("lda", "hdp"), deviant_z = 3, ...)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted.

- model:

  One of `"lda"` or `"hdp"`.

- deviant_z:

  Robust z-score above which a sample is flagged as deviating from its
  group's fixed ratio (default 3).

- ...:

  Passed to
  [compare_topics](https://marine-ecologist.github.io/symbayes/reference/compare_topics.md)
  (thresholds).

## Value

Invisibly, a list: `pairs` (the compare_topics table), `groups` (data
frame mapping each component to a suggested group id), `n_units`
(suggested biological unit count), and `deviants` (samples breaking
their group's fixed ratio).

## Details

For each suggested group it also flags **deviant samples** – those whose
within-group ratio departs from the group's fixed ratio by more than
`deviant_z` robust SDs. These are the samples that break the
intragenomic pattern: either genuine mixing exceptions or low-depth
noise, surfaced for inspection rather than absorbed into the group
verdict.

To actually merge a suggested grouping into a coarser membership, call
[merge_topics](https://marine-ecologist.github.io/symbayes/reference/merge_topics.md)
explicitly with the grouping you choose.

## Examples

``` r
if (FALSE) { # \dontrun{
scr <- screen_intragenomic(sp, model = "lda")
scr$n_units
scr$groups
scr$deviants
} # }
```

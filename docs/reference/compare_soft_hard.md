# Compare a soft (fractional) model against SymPortal's profiles via the fuzzy Rand index

Quantifies how much a fractional model (LDA or HDP) agrees with
SymPortal's ITS2 type-profile assignments at the level of *sample
pairs*, using the fuzzy / probabilistic Rand index (Campello 2007;
Huellermeier et al. 2011). This is a principled upgrade over a scalar
mixedness correlation: it compares the full relational structure (which
samples cluster with which, fractionally) rather than a per-sample
summary.

## Usage

``` r
compare_soft_hard(sp, model = c("lda", "hdp"), mixed_thresh = 0.7)
```

## Arguments

- sp:

  A `symbayes` object with the chosen model fitted and SymPortal
  profiles present.

- model:

  One of `"lda"` or `"hdp"` (fractional models). DMM is hard- assignment
  and not meaningful here.

- mixed_thresh:

  A sample is "mixed" if its dominant proportion is below this (default
  0.7). Applied to BOTH the model's `theta` and SymPortal's
  multi-profile proportions, so "SymPortal mixed" means SymPortal
  assigned the sample multiple profiles with non-trivial mass – not
  whether the profile name contains a separator (SymPortal names list
  intragenomic DIVs and almost always contain '-' or '/').

## Value

Invisibly, a list with `fuzzy_rand`, `adjusted_fuzzy_rand`,
`sample_divergence` (data frame, most-divergent first, with
`model_mixed` and `sp_mixed` flags), and `disagreement` (counts of the
two disagreement types: SymPortal missed mixing vs SymPortal
over-split).

## Details

For each pair of samples the probability they are co-assigned is the dot
product of their membership vectors (`theta %*% t(theta)`). SymPortal's
multi-profile proportions are treated the same way, giving two soft
co-membership matrices. The fuzzy Rand index is
`1 - mean|EQ_model - EQ_sp|` over all pairs; the adjusted version
corrects for chance agreement (the raw index is inflated by the many
pairs that are clearly separate under both).

Also returns per-sample co-membership divergence (where the two methods
disagree most — typically the admixed samples) and a disagreement
breakdown separating "SymPortal missed mixing the model found" from
"SymPortal invented a compound profile the model sees as pure".

## Examples

``` r
if (FALSE) { # \dontrun{
res <- compare_soft_hard(sp, model = "lda")
res$fuzzy_rand
head(res$sample_divergence)
} # }
```

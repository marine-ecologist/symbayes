# Depth-aware overdispersion test (rebuilt on real read counts)

Rebuilt to use **actual sequence read counts** rather than
`theta * depth`. Each component is represented by its *defining
sequences* (those it clearly owns in its exemplar). For a pair (A, B),
`k_i` = sample i's summed reads of A's defining sequences and `n_i` =
reads of A's + B's defining sequences – genuine integer read counts, so
the binomial sampling model is legitimate.

## Usage

``` r
test_overdispersion(
  sp,
  model = c("lda", "hdp"),
  min_frac = 0.05,
  min_samples = 5,
  min_pair_reads = 100,
  rho_mixing = 0.05,
  top_n = 15,
  owner_frac = 0.5
)
```

## Arguments

- sp:

  A `symbayes` object with the chosen model fitted.

- model:

  One of `"lda"` or `"hdp"`.

- min_frac:

  Component presence threshold to count a sample as co-occurring
  (default 0.05).

- min_samples:

  Minimum co-occurring samples to test a pair (default 5).

- min_pair_reads:

  Minimum `n_i` (pair reads) to include a sample (default 100).

- rho_mixing:

  Beta-binomial rho above which a pair is called `"mixing"` (default
  0.05). Between a small floor and this, `"ambiguous"`.

- top_n, owner_frac:

  Passed to the internal defining-sequence assignment.

## Value

Invisibly, a data frame: `comp_a`, `comp_b`, `n_samples`, `n_seqs_a`,
`n_seqs_b`, `mean_ratio`, `rho`, `phi`, `p_value`, `call`.

## Details

Overdispersion is estimated as the beta-binomial intra-class correlation
\\\rho\\ (method of moments), which bounds in \\\[0, 1)\\ and properly
separates true cross-sample ratio variance from within-sample sampling
variance. \\\rho \approx 0\\ = fixed ratio (consistent with intragenomic
variation); \\\rho\\ clearly above 0 = variable ratio (consistent with
community mixing). A dispersion index `phi` (Pearson X^2 / df) is
reported alongside as a cross-check – now well-behaved because `n_i` are
real counts.

## Examples

``` r
if (FALSE) { # \dontrun{
od <- test_overdispersion(sp, model = "lda")
subset(od, call == "mixing")
} # }
```

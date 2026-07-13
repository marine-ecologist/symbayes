# Combined pairwise topic comparison: redundancy, mixing, and ratio stability

One row per co-occurring component pair, synthesising the three
diagnostics needed to decide what a pair *is*:

## Usage

``` r
compare_topics(
  sp,
  model = c("lda", "hdp"),
  redundant_thresh = 0.85,
  rho_intra = 0.03,
  rho_mixing = 0.05,
  min_frac = 0.05,
  min_samples = 5,
  min_pair_reads = 100
)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted.

- model:

  One of `"lda"` or `"hdp"`.

- redundant_thresh:

  Cosine at/above which a pair is `"redundant"` (default 0.85).

- rho_intra:

  Rho at/below which a co-occurring pair is `"intragenomic"` (default
  0.03).

- rho_mixing:

  Rho at/above which a co-occurring pair is `"mixing"` (default 0.05).

- min_frac, min_samples, min_pair_reads:

  Passed to the overdispersion test.

## Value

Invisibly, a data frame of pairs with all diagnostics and `verdict`,
ordered by verdict priority then strength.

## Details

- `cosine`:

  Exemplar (beta) cosine similarity. High (\> `redundant_thresh`) =
  redundant topics (the same signature split in two).

- `n_cooccur`, `mean_ratio`, `ratio_sd`:

  How often the pair co-occurs and how stable their relative proportion
  is across those samples.

- `rho`:

  Beta-binomial overdispersion from
  [test_overdispersion](https://marine-ecologist.github.io/symbayes/reference/test_overdispersion.md)
  on real reads. Low (~0) = fixed ratio (intragenomic); high = variable
  (mixing).

- `verdict`:

  A synthesised call (see below).

The verdict combines all three:

- `"redundant"` – cosine \>= `redundant_thresh` (merge candidates).

- `"intragenomic"` – co-occur often at a fixed ratio (rho \<=
  `rho_intra`): likely one symbiont's intragenomic variants split across
  topics.

- `"mixing"` – co-occur with a variable ratio (rho \>= `rho_mixing`):
  genuine community admixture.

- `"ambiguous"` – co-occur but rho between the thresholds.

- `"distinct"` – rarely co-occur (\< `min_samples`): separable types.

Compositional evidence alone cannot prove intragenomic vs a
fixed-proportion obligate association; a single-copy marker (psbA) is
the definitive arbiter.

## Examples

``` r
if (FALSE) { # \dontrun{
ct <- compare_topics(sp, model = "lda")
subset(ct, verdict == "intragenomic")
} # }
```

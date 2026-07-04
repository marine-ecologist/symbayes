# Fit a Latent Dirichlet Allocation topic model

Fits LDA by Gibbs sampling. Unlike DMM, LDA gives *fractional*
membership: each sample is a mixture of topics, so within-sample
admixture is represented directly.

## Usage

``` r
fit_lda(
  sp,
  k = NULL,
  k_range = 2:15,
  redundant_thresh = 0.85,
  burnin = 1000,
  iter = 2000,
  thin = 5,
  seed = 42
)
```

## Arguments

- sp:

  A filtered `symbayes` object.

- k:

  Fixed number of topics, or `NULL` (default) to select adaptively.

- k_range:

  Integer vector of K to search when `k = NULL` (default 2:15).

- redundant_thresh:

  Cosine similarity above which two topics are deemed redundant (default
  0.85).

- burnin, iter, thin:

  Gibbs sampler controls (defaults 1000, 2000, 5).

- seed:

  Random seed (default 42).

## Value

The `symbayes` object with `sp$lda` populated.

## Details

Topic number selection is threshold-free-oriented: with `k = NULL`
(default) the function searches `k_range` and picks the **largest k with
no redundant topic pairs** (cosine similarity `< redundant_thresh`).
This avoids the "shadow topics" that perplexity-optimal K tends to
produce. Pass an integer `k` to fix it. The search table is stored in
`sp$lda$meta$selection` and can be plotted with
[lda_plot_ksearch](https://marine-ecologist.github.io/symbayes/reference/lda_plot_ksearch.md).

Writes a membership object to `sp$lda`. `meta` contains `entropy`,
`max_prop`, `n_topics` (per sample), `similarity` (topic cosine matrix),
`redundant` (flagged pairs), and `selection` (the K-search table).

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- fit_lda(sp)              # adaptive K
sp <- fit_lda(sp, k = 9)       # fixed K
lda_plot_ksearch(sp)
lda_plot_similarity(sp)
} # }
```

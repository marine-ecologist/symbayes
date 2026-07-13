# Plot the overdispersion result for one component pair (real read counts)

Each co-occurring sample's observed A-share from real defining-sequence
reads, vs its pair read depth, with the binomial fixed-ratio null
envelope (+/-2 SE around the pooled ratio). Points tracking the band at
all depths = fixed ratio (intragenomic); systematic spread beyond it,
especially at high depth where the band is tight, = genuine mixing.

## Usage

``` r
plot_overdispersion(
  sp,
  model = c("lda", "hdp"),
  comp_a,
  comp_b,
  min_frac = 0.05,
  min_pair_reads = 100,
  top_n = 15,
  owner_frac = 0.5
)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted.

- model:

  One of `"lda"` or `"hdp"`.

- comp_a, comp_b:

  Component labels (e.g. `"Topic_2"`, `"Topic_9"`).

- min_frac:

  Presence threshold (default 0.05).

- min_pair_reads:

  Minimum pair reads to include a sample (default 100).

- top_n, owner_frac:

  Defining-sequence assignment controls.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  plot_overdispersion(sp, "lda", "Topic_2", "Topic_9")  # \dontrun{}
```

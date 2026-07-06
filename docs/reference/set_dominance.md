# Relabel dominant assignment with a minimum-dominance threshold

The hard `dominant` label is `which.max(theta)` – a plurality that can
be as low as ~1/k. For near-boundary samples this label is fragile and
overstates confidence (a 0.35/0.34 split is labelled definitively by its
0.35 winner). This function recomputes the dominant label so that
samples whose top component is below `min_dominance` are labelled
`"mixed"` instead of forced to a plurality winner. It changes only the
*label* (the `<model>_dominant` metadata column and the membership
`dominant` factor); `theta` – the model's actual fractional output – is
untouched.

## Usage

``` r
set_dominance(sp, model = c("dmm", "lda", "hdp"), min_dominance = 0.5)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- min_dominance:

  Minimum top-component proportion to keep a hard label; below this the
  sample is labelled `"mixed"` (default 0.5).

## Value

The updated `sp` object. Adds a `<model>_dominant_thr` factor to
metadata (levels = group labels plus `"mixed"`) and stores the threshold
in `sp$<model>$meta$min_dominance`. The original `<model>_dominant` is
left intact for backwards compatibility.

## Details

There is no "correct" threshold: this is a reporting-honesty choice, not
a model parameter. `0.5` (call a sample mixed unless one component holds
an outright majority) is a natural, explainable default. Note SymPortal
has no equivalent control, because it never produces a proportion to
threshold – it assigns each sample one profile (possibly a compound name
that already absorbs the mixing). This function is symbayes
acknowledging admixture explicitly where SymPortal cannot.

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- set_dominance(sp, model = "hdp", min_dominance = 0.5)
table(sp$metadata$hdp_dominant_thr)   # how many samples are "mixed"
} # }
```

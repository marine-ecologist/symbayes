# Admixture plot: mixed samples, sequence mix over boxed component structure

Filters to samples that are genuinely *mixed* under a fitted model (at
least `min_groups` groups above `min_frac`) and returns a two-panel
figure sharing one x-axis:

## Usage

``` r
plot_admixtures(
  sp,
  model = "hdp",
  min_frac = 0.05,
  min_groups = 2,
  seq_min_pct = 1,
  sample_label = "sample_name",
  box_colour = "grey15",
  label_profile = TRUE,
  combine = TRUE
)
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- min_frac:

  Minimum group proportion to count as present (default 0.05).

- min_groups:

  Minimum number of present groups for a sample to be shown (default 2;
  i.e. only mixed samples).

- seq_min_pct:

  Minimum percent in any admixed sample to show a sequence individually
  in the top panel (default 1).

- sample_label:

  Metadata column for x labels (default `"sample_name"`).

- box_colour:

  Outline colour for each within-sample component segment (default
  `"grey15"`).

- label_profile:

  If `TRUE` (default) and SymPortal profiles are present, annotate the
  bottom panel with each sample's dominant SymPortal profile.

- combine:

  If `TRUE` (default) return a patchwork of both panels; if `FALSE`
  return a list with `seq` and `admix` ggplots.

## Value

A patchwork object (or a list if `combine = FALSE`).

## Details

- Top:

  ITS2 **sequence** composition of the admixed samples (a single stacked
  barplot, no facet) – the raw mixture SymPortal sees.

- Bottom:

  **Component** proportions (`theta`) for the same samples, with every
  contributing segment outlined so the within-sample admixture is
  demarcated.

Together they show what SymPortal would collapse to a single compound
profile name (top) versus the fractional community structure the mixture
model recovers (bottom). DMM is near-degenerate here (hard assignment),
so this is most informative for `"lda"` and `"hdp"`.

## Examples

``` r
if (FALSE)  plot_admixtures(sp, model = "lda")  # \dontrun{}
```

# Four-row SymPortal-vs-model comparison over all samples

Stacks four aligned panels sharing one sample order (all samples,
ordered by the chosen model's dominant group then confidence):

## Usage

``` r
plot_symportal_comparison(
  sp,
  model = "hdp",
  sample_label = "sample_name",
  box_colour = "grey15",
  show_x = FALSE,
  screen = FALSE,
  majority = 0.6,
  max_legend = 15,
  n_seqs = 20,
  residual = c("renormalise", "show"),
  use_names = FALSE
)
```

## Arguments

- sp:

  A `symbayes` object with the chosen model fitted and SymPortal
  profiles present.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- sample_label:

  Metadata column for x labels (default `"sample_name"`).

- box_colour:

  Outline colour for admixture segments (default `"grey15"`).

- show_x:

  Draw sample labels on the bottom row (default `FALSE`; with many
  samples they are unreadable).

- screen:

  If `TRUE`, row 1 uses *screen-aware* assignment
  ([screen_intragenomic](https://marine-ecologist.github.io/symbayes/reference/screen_intragenomic.md)):
  intragenomic topic pairs are collapsed to one unit before labelling,
  and samples that genuinely mix distinct units (a "mixing"-verdict
  pair) are shown as a distinct "true mix" band rather than forced to a
  single plurality topic. Default `FALSE` (raw `which.max`).

- majority:

  Minimum merged-unit proportion for a sample to be called that unit
  under `screen = TRUE`; below it the sample is "true mix" (default
  0.6).

- max_legend:

  Maximum number of SymPortal profiles to list in the legend (the most
  frequent dominant profiles). All profiles are still coloured; only the
  legend is capped, to stay on-panel for large datasets (default 15).

- residual:

  How to show model mass assigned to no component in the admixture row:
  `"renormalise"` (default; rescale each bar to sum to 1) or `"show"`
  (grey Residual segment so bars reach 1.0).

- use_names:

  If `TRUE`, label topics/components in the model legend with their
  SymPortal-style composition names (from
  [name_topics](https://marine-ecologist.github.io/symbayes/reference/name_topics.md))
  instead of `Topic_N` / `HDP_N` (default `FALSE`).

## Value

A patchwork object (four rows).

## Details

1.  **Model assignment** – a single colour band per sample showing the
    model's dominant group.

2.  **Model admixture** – boxed component proportions (`theta`); each
    outlined segment is one component's share within the sample.

3.  **SymPortal assignment** – a single colour band per sample showing
    its dominant SymPortal profile.

4.  **SymPortal mix** – stacked SymPortal multi-profile proportions;
    samples with more than one profile show the admixture SymPortal
    records.

Rows 1-2 are the model's view (hard then soft); rows 3-4 are SymPortal's
(hard then soft). Reading a column top-to-bottom shows, for one sample,
how the model and SymPortal each assign and each mix it.

## Examples

``` r
if (FALSE)  plot_symportal_comparison(sp, model = "hdp")  # \dontrun{}
```

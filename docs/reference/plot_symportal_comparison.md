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
  show_x = FALSE
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

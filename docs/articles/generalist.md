# 5. Generalist plots across models

The generalist plotting functions read the standardised membership
schema and so work on **any** fitted model via a
`model = c("dmm", "lda", "hdp")` argument. This lets you render the same
view for each model and compare like with like.

``` r
library(symbayes)
sp <- import("run.seqs.absolute.abund_only.txt",
             profs_abund = "run.profiles.absolute.abund_only.txt",
             profs_meta  = "run.profiles.meta_only.txt",
             colour_dict = "color_dict_post_med.json")
sp <- filter_samples(sp)
sp <- fit_dmm(sp)
sp <- fit_lda(sp)
sp <- fit_hdp(sp)
```

## PCoA ordination

Colour a Bray-Curtis PCoA by any model’s dominant group, or by a
metadata variable:

``` r
plot_pcoa(sp, model = "dmm")
plot_pcoa(sp, model = "lda")
plot_pcoa(sp, model = "hdp")

plot_pcoa(sp, colour_by = "sp_clade",
          colours = c(A = "#FF8C00", C = "#20B2AA", D = "#9370DB"))
```

## Sequence-composition barplot

Faceted by a model’s dominant group, annotated with assignment
confidence:

``` r
plot_barplot(sp, model = "dmm")
plot_barplot(sp, model = "hdp", min_pct = 1)
```

## Exemplars

Each group’s archetypal composition (`beta`):

``` r
plot_exemplars(sp, model = "lda")
plot_exemplars(sp, model = "hdp")
```

## Membership heatmap

Near-binary for DMM, gradient for LDA/HDP — the visual signature of hard
vs fractional assignment:

``` r
plot_membership(sp, model = "dmm")   # sharp
plot_membership(sp, model = "lda")   # graded
plot_membership(sp, model = "hdp")
```

## Contingency vs SymPortal profiles

Each model’s groups against SymPortal profiles (rows grouped by clade);
tile colour is mean confidence, number is sample count:

``` r
plot_contingency(sp, model = "dmm")
plot_contingency(sp, model = "hdp")
```

## Top characterising sequences

``` r
plot_top_seqs(sp, model = "dmm", top_n = 8)
plot_top_seqs(sp, model = "lda", top_n = 8)
```

## Rendering all three side by side

Because the interface is uniform, you can loop:

``` r
library(patchwork)
plots <- lapply(c("dmm", "lda", "hdp"),
                function(m) plot_membership(sp, model = m))
wrap_plots(plots, ncol = 1)
```

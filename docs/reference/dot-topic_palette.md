# Stable, distinguishable colour palette for topics / components

Assigns each topic/component an evenly-spaced hue across the full colour
wheel so all topics are maximally distinguishable, while keeping the
assignment DETERMINISTIC by label (same label -\> same colour across
plots, models, and runs). Hue order is set by a deterministic hash of
the label, so colours are stable but spread out rather than clustered.

## Usage

``` r
.topic_palette(labels, names_map = NULL)
```

## Arguments

- labels:

  Character vector of topic/component labels.

- names_map:

  Accepted for API compatibility; unused here.

## Value

A named character vector of hex colours, one per label.

## Details

Clade is not used to anchor hue: in datasets dominated by one clade,
clade-anchoring collapses most topics onto near-identical colours. This
palette prioritises telling topics apart. (Clade grouping remains
available via
[.clade_palette](https://marine-ecologist.github.io/symbayes/reference/dot-clade_palette.md)
for clade-level fills such as
[plot_top_seqs](https://marine-ecologist.github.io/symbayes/reference/plot_top_seqs.md).)

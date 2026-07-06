# PCoA ordination coloured by a fitted model's grouping (or any variable)

Principal coordinates analysis of Bray-Curtis distances among samples,
coloured by the dominant group of a fitted model, or by any metadata
column.

## Usage

``` r
plot_pcoa(sp, model = "dmm", colour_by = NULL, colours = NULL)
```

## Arguments

- sp:

  A `symbayes` object.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`. Colours points by that model's
  dominant group. Ignored if `colour_by` is set.

- colour_by:

  Optional metadata column name to colour by instead (e.g.
  `"host_species"`, `"sp_clade"`).

- colours:

  Optional named colour vector.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_pcoa(sp, model = "hdp")
plot_pcoa(sp, colour_by = "sp_clade")
} # }
```

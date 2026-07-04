# Plot DMM model selection (Laplace, AIC, BIC vs K)

Shows the three information criteria across the candidate component
numbers, with the Laplace-selected K marked. Laplace (the model evidence
approximation) is the primary criterion.

## Usage

``` r
dmm_plot_selection(sp)
```

## Arguments

- sp:

  A `symbayes` object with DMM fitted.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  dmm_plot_selection(sp)  # \dontrun{}
```

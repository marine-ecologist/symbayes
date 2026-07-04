# Plot the HDP posterior distribution of component number

Shows two quantities honestly: the *raw* cluster count the Gibbs sampler
used at each posterior sample (a distribution), and the number of
*extracted* stable components after merging cosine-similar and dropping
transient clusters (a single line). The raw-to-extracted reduction is
expected, not a loss.

## Usage

``` r
hdp_plot_numcomp(sp)
```

## Arguments

- sp:

  A `symbayes` object with HDP fitted.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE)  hdp_plot_numcomp(sp)  # \dontrun{}
```

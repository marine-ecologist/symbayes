# Simulate samples and test DMM assignment (blend / perturb / from_alpha)

Three modes for probing DMM behaviour: `from_alpha` draws from each
component's Dirichlet then multinomial; `blend` mixes two components
across a gradient to trace decision boundaries; `perturb` reassigns a
fraction of reads in real samples to test robustness.

## Usage

``` r
dmm_simulate(
  sp,
  mode = c("from_alpha", "blend", "perturb"),
  n_sim = 100,
  lib_size = 5000,
  perturb_frac = 0.1,
  blend_comps = c(1, 2),
  blend_range = seq(0, 1, 0.1),
  seed = 42
)
```

## Arguments

- sp:

  A `symbayes` object with DMM fitted.

- mode:

  One of `"from_alpha"`, `"blend"`, `"perturb"`.

- n_sim:

  Number of simulated samples (default 100).

- lib_size:

  Reads per simulated sample (default 5000).

- perturb_frac:

  Fraction of reads to reassign for `"perturb"` (default 0.1).

- blend_comps:

  Two component indices to mix for `"blend"` (default `c(1, 2)`).

- blend_range:

  Mixing proportions for `"blend"` (default `seq(0, 1, 0.1)`).

- seed:

  Random seed (default 42).

## Value

A data frame of simulated samples with predicted assignments.

## Examples

``` r
if (FALSE) { # \dontrun{
dmm_simulate(sp, mode = "blend", blend_comps = c(1, 4))
dmm_simulate(sp, mode = "perturb", perturb_frac = 0.1)
} # }
```

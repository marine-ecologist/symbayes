# Fit a Hierarchical Dirichlet Process

The threshold-free, mixture-aware core of `symbayes`. Unlike DMM and
LDA, the number of community types is **inferred from the data** rather
than specified: the DP concentration parameters (with Gamma hyperpriors)
govern component number, and samples may be admixed. Uses the hdp
package (Roberts). See <https://github.com/nicolaroberts/hdp>.

## Usage

``` r
fit_hdp(
  sp,
  group_var = NULL,
  n_chains = 4,
  burnin = 5000,
  n = 100,
  space = 50,
  initcc = 10,
  alphaa = 1,
  alphab = 1,
  cos_merge = 0.9,
  min_sample = 1,
  seed = 42
)
```

## Arguments

- sp:

  A filtered `symbayes` object.

- group_var:

  Optional metadata column defining a middle DP layer, or `NULL`
  (default) for a flat, covariate-free structure.

- n_chains:

  Independent Gibbs chains (default 4).

- burnin:

  Burn-in iterations per chain (default 5000).

- n:

  Posterior samples collected per chain (default 100).

- space:

  Iterations between collected samples (default 50).

- initcc:

  Initial cluster count (default 10; a starting point only).

- alphaa, alphab:

  Gamma hyperprior (shape, rate) on DP concentration (default 1, 1).

- cos_merge:

  Cosine similarity above which raw clusters are merged during
  extraction (default 0.90; the hdp default). Higher retains more
  components. Do not set very low or everything collapses.

- min_sample:

  Minimum posterior samples a component must appear in to be retained
  (default 1).

- seed:

  Base random seed (chains use `seed + chain`).

## Value

The `symbayes` object with `sp$hdp` populated.

## Details

With `group_var = NULL` the structure is flat (one parent DP, one child
DP per sample) and assignment is covariate-free. Supplying a metadata
column (e.g. `"host_species"`) inserts a middle DP layer so samples
within a group share strength preferentially. **Caution:** grouping by a
covariate you later want to *test* (e.g. host effect on symbiont
community) introduces circularity; use flat structure for host-blind
claims, grouped only as an explicitly-labelled predictive enhancement.

Writes a membership object to `sp$hdp`. `theta` excludes the residual
(unassigned) mass, which is kept in `meta$residual`. `meta` also holds
`entropy`, `max_prop`, `n_comp` per sample, `raw_ncl` (sampler cluster
count), and the `multi` chain object.

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- fit_hdp(sp)                       # flat, threshold-free
hdp_plot_numcomp(sp)
hdp_vs_symportal(sp)
} # }
```

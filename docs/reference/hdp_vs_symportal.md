# Quantify SymPortal over-classification relative to HDP

The core method metric. Compares SymPortal's profile count to HDP's
inferred component count and, crucially, partitions SymPortal profiles
into three mutually exclusive classes so they do not contaminate one
another:

## Usage

``` r
hdp_vs_symportal(sp, over_split_min_frac = 0.7, mixture_thresh = 0.7)
```

## Arguments

- sp:

  A `symbayes` object with HDP fitted and SymPortal profiles present.

- over_split_min_frac:

  Minimum dominant-component share for a profile to count as over-split
  (default 0.7).

- mixture_thresh:

  Dominant-component share below which a profile is a mixture (default
  0.7).

## Value

Invisibly, a list with `n_sp_profiles`, `n_hdp_comp`, `oversplit_ratio`,
`mapping`, `redundant_splits`, `mixture_profiles` (profile-level),
`mixture_samples` (one row per sample whose dominant profile is a
mixture, with ranked `comp_i`/`frac_i` columns showing which HDP
components co-occur in that individual sample), and `orphan_profiles`.

## Details

- Over-splitting:

  Several profiles that *confidently* map to one HDP component
  (`dominant_frac >= over_split_min_frac`). SymPortal split one real
  community type into multiple named profiles.

- Mixtures:

  Profiles whose samples spread across components
  (`dominant_frac < mixture_thresh`). SymPortal called an admixture a
  unique type; these are not force-assigned.

- Orphans:

  Not confidently in any component and not clear mixtures (usually rare,
  n = 1). Reported separately and excluded from the over-splitting
  count.

The partition matters: a rare profile whose single sample has its
largest component at only a few percent must not be counted as
"over-split into" that component.

## Examples

``` r
if (FALSE)  res <- hdp_vs_symportal(sp)  # \dontrun{}
```

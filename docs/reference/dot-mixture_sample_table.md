# Build a per-sample HDP mixture breakdown table

One row per sample, listing the top contributing HDP components and
their proportions as ranked columns.

## Usage

``` r
.mixture_sample_table(
  sp,
  theta,
  dom_prof,
  name_map,
  mixture_uids,
  top_n = 3,
  min_frac = 0.05
)
```

## Arguments

- sp:

  A `symbayes` object.

- theta:

  HDP sample x component matrix.

- dom_prof:

  Named vector: sample_uid -\> dominant profile UID.

- name_map:

  Named vector: profile UID -\> profile name.

- mixture_uids:

  Sample UIDs whose dominant profile is a flagged mixture.

- top_n:

  Number of ranked component columns to emit (default 3).

- min_frac:

  Minimum proportion to include a component (default 0.05).

## Value

A data frame: sample_uid, sample_name, sp_profile, n_comp, then
comp_1/frac_1 ... comp_top_n/frac_top_n.

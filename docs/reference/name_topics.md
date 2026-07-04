# Name a fitted model's topics/components in SymPortal style

Renders each group's exemplar sequence distribution (`beta`) as a
SymPortal-style name, so model groups can be compared to SymPortal
profiles in the same convention. Returns a data frame mapping the group
label to its name, dominant clade, and the exemplar's Shannon entropy (a
diversity / potential-over-naming flag: high-entropy exemplars produce
long names).

## Usage

``` r
name_topics(
  sp,
  model = c("dmm", "lda", "hdp"),
  codominant_ratio = 0.5,
  min_weight = 0.05,
  max_seqs = 6,
  assign = FALSE
)
```

## Arguments

- sp:

  A `symbayes` object with the chosen model fitted.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- codominant_ratio:

  Co-dominance cutoff for `/` (default 0.5).

- min_weight:

  Minimum sequence weight to include (default 0.05).

- max_seqs:

  Maximum sequences per name (default 6).

- assign:

  If `TRUE`, also writes the names onto `sp` as `sp$<model>$meta$names`
  and returns the updated `sp`; if `FALSE` (default) returns the naming
  data frame.

## Value

A data frame (group, name, clade, entropy, n_samples), or the updated
`sp` object if `assign = TRUE`.

## Details

The names are a labelling convenience for display and comparison. A
group's identity is its full distribution and its samples' membership,
not this string. Groups sharing a backbone but differing in minor
variants can receive similar names – the same mechanism behind SymPortal
over-splitting – so do not treat similar names as evidence groups are
redundant (use the cosine similarity from
[lda_plot_similarity](https://marine-ecologist.github.io/symbayes/reference/lda_plot_similarity.md)
or the HDP component distances for that).

## Examples

``` r
if (FALSE) { # \dontrun{
name_topics(sp, model = "lda")
sp <- name_topics(sp, model = "hdp", assign = TRUE)
sp$hdp$meta$names
} # }
```

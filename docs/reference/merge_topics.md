# Merge topics into a coarser membership (explicit, user-directed)

Applies a grouping (e.g. from
[screen_intragenomic](https://marine-ecologist.github.io/symbayes/reference/screen_intragenomic.md))
by summing the theta columns of components in the same group, producing
a coarser membership matrix. This is deliberately a separate, explicit
call – the package never merges automatically, because merging on a
compositional threshold is a judgement the user must own (and ideally
confirm with psbA).

## Usage

``` r
merge_topics(sp, model = c("lda", "hdp"), groups)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted.

- model:

  One of `"lda"` or `"hdp"`.

- groups:

  Either a data frame with `component` and `group` columns (as returned
  by `screen_intragenomic()$groups`), or a named integer vector mapping
  component label -\> group id.

## Value

A matrix (samples x merged groups) of summed memberships. Column names
list the merged components. Returned, not written into `sp`.

## Examples

``` r
if (FALSE) { # \dontrun{
scr <- screen_intragenomic(sp, model = "lda")
merged_theta <- merge_topics(sp, model = "lda", groups = scr$groups)
} # }
```

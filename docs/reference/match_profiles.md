# Match model topics/components to SymPortal profiles by composition

String-matching topic names to SymPortal profile names fails: the same
community can be written `C22-C3-C22ad` (topic) and `C22-C22ad-C3-C22ah`
(SymPortal) – overlapping sequences, different order and set. This
function instead matches by *composition*: it represents each SymPortal
profile as a sequence vector (from the profile's defining sequences, via
`prof_meta` if available, else from the mean sequence composition of
samples dominant for that profile) and computes cosine similarity to
each topic's exemplar (`beta`). Each topic is matched to its most
similar profile.

## Usage

``` r
match_profiles(sp, model = c("dmm", "lda", "hdp"), top_n = 3)
```

## Arguments

- sp:

  A `symbayes` object with the model fitted and SymPortal profiles
  present.

- model:

  One of `"dmm"`, `"lda"`, `"hdp"`.

- top_n:

  Number of best-matching profiles to report per topic (default 3).

## Value

A data frame: one row per topic-profile match, with `group`,
`topic_name`, `sp_profile`, `cosine`, and match rank.

## Examples

``` r
if (FALSE)  match_profiles(sp, model = "lda")  # \dontrun{}
```

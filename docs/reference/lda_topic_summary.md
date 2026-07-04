# Text summary of LDA topics with SymPortal / host cross-reference

Prints, for each topic: dominant clade, number of samples, mean entropy
and purity, top host species (if present), top SymPortal profiles (if
present), and the exemplar sequence composition.

## Usage

``` r
lda_topic_summary(sp, top_n = 8)
```

## Arguments

- sp:

  A `symbayes` object with LDA fitted.

- top_n:

  Number of top sequences per topic (default 8).

## Value

Invisibly `NULL`; called for its printed output.

## Examples

``` r
if (FALSE)  lda_topic_summary(sp)  # \dontrun{}
```

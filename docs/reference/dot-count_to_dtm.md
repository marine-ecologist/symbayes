# Convert a count matrix to a topicmodels DocumentTermMatrix

Builds the DTM directly as a
[`slam::simple_triplet_matrix`](https://rdrr.io/pkg/slam/man/matrix.html)
with the `DocumentTermMatrix` class and weighting attribute, avoiding a
hard dependency on
[`tm::as.DocumentTermMatrix`](https://rdrr.io/pkg/tm/man/matrix.html).

## Usage

``` r
.count_to_dtm(count_mat)
```

## Arguments

- count_mat:

  Integer sample x sequence matrix

## Value

A `DocumentTermMatrix`

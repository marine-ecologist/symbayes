# Name a single exemplar vector in SymPortal style

`/` separates co-dominant sequences (within `codominant_ratio` of the
top); `-` appends less-abundant associates in descending order, down to
`min_weight`.

## Usage

``` r
.name_exemplar(
  beta_row,
  codominant_ratio = 0.5,
  min_weight = 0.05,
  max_seqs = 6
)
```

## Arguments

- beta_row:

  Named numeric vector (sequence -\> weight), sums to ~1.

- codominant_ratio:

  Sequences whose weight is \>= this fraction of the single largest are
  treated as co-dominant and joined with `/` (default 0.5).

- min_weight:

  Minimum weight to include a sequence at all (default 0.05).

- max_seqs:

  Maximum sequences in the name (default 6).

## Value

A SymPortal-style name string.

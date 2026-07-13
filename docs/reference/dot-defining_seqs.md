# Assign defining sequences to each component from its exemplar

Each component's defining sequences are those where it holds the largest
share of that sequence's weight across components (argmax over the beta
columns), restricted to sequences the component emphasises (its own top
sequences). Sequences that are not clearly owned by one component are
left unassigned, so a shared backbone sequence does not corrupt a pair's
ratio.

## Usage

``` r
.defining_seqs(beta, top_n = 15, owner_frac = 0.5)
```

## Arguments

- beta:

  components x sequences exemplar matrix.

- top_n:

  Consider each component's top-N sequences as candidates (default 15).

- owner_frac:

  A sequence is "owned" by a component only if that component holds at
  least this fraction of the sequence's total weight across components
  (default 0.5 = clear majority).

## Value

A named character vector: sequence -\> owning component (or dropped).

# Infer Symbiodiniaceae clade from sequence names

SymPortal names sequences by clade (leading letter, e.g. `C3`, `D1`,
`A6b`) or by clade suffix for named intragenomic variants (e.g.
`15443_A`).

## Usage

``` r
.infer_clade(seq_names)
```

## Arguments

- seq_names:

  Character vector of sequence names

## Value

Named character vector of clade letters (A-D, F, G) or NA

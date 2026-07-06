# Build a colour vector from a SymPortal colour dictionary plus fallback

Build a colour vector from a SymPortal colour dictionary plus fallback

## Usage

``` r
.build_colours(seq_names, col_dict = NULL)
```

## Arguments

- seq_names:

  Character vector of sequence names to colour

- col_dict:

  Named list/vector mapping sequence names to hex colours (from
  SymPortal's `color_dict_post_med.json`); may be `NULL`

## Value

Named character vector of hex colours, including a grey `Other`

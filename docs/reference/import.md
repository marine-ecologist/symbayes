# Import SymPortal output into a `symbayes` object

Reads a SymPortal post-MED sequence count table (the only required
input) and, optionally, ITS2 type profiles, sample metadata, and the
SymPortal colour dictionary. Returns a `symbayes` object that flows
through the rest of the pipeline.

## Usage

``` r
import(
  seqs_abund,
  profs_abund = NULL,
  profs_meta = NULL,
  seqs_meta = NULL,
  sample_sheet = NULL,
  colour_dict = NULL,
  sheet_skip = 1
)
```

## Arguments

- seqs_abund:

  Path to `*.seqs.absolute.abund_only.txt` (samples x post-MED sequence
  variants). Required.

- profs_abund:

  Path to `*.profiles.absolute.abund_only.txt` (samples x ITS2 type
  profiles). Optional; enables profile comparison.

- profs_meta:

  Path to `*.profiles.meta_only.txt` (profile UID, name, clade).
  Optional.

- seqs_meta:

  Path to `*.seqs.absolute.meta_only.txt` (per-sample QC). Optional.

- sample_sheet:

  Path to a SymPortal datasheet `.xlsx`. Optional; needs the readxl
  package.

- colour_dict:

  Path to `color_dict_post_med.json`. Optional; used for
  SymPortal-consistent sequence colours in bar plots.

- sheet_skip:

  Header rows to skip in the datasheet (default 1).

## Value

A `symbayes` object: a list with `count_mat`, `prof_mat`, `prof_meta`,
`metadata`, `col_dict`, `seq_clades`, and `filtered = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
sp <- import(
  seqs_abund  = "run.seqs.absolute.abund_only.txt",
  profs_abund = "run.profiles.absolute.abund_only.txt",
  profs_meta  = "run.profiles.meta_only.txt",
  colour_dict = "color_dict_post_med.json"
)
sp
} # }
```

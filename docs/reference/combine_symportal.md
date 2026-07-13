# Merge multiple SymPortal runs into one `symbayes` object

Combines several imported runs so models can be fit on a larger dataset.
Sequences are aligned by name and zero-filled where absent; samples are
stacked. **Profiles are retained and matched across runs** by their
stable accession signature (see
[.profile_signature](https://marine-ecologist.github.io/symbayes/reference/dot-profile_signature.md))
rather than by run-specific UID – so the same biological profile is
recognised as one entity even when its run-dependent minor-DIV tail
differs.

## Usage

``` r
combine_symportal(
  ...,
  run_ids = NULL,
  on_collision = c("suffix", "error", "keep_first"),
  numeric_warn_frac = 0.2,
  n_lead = 3,
  merge_profiles = TRUE
)
```

## Arguments

- ...:

  Two or more `symbayes` objects (from
  [import](https://marine-ecologist.github.io/symbayes/reference/import.md)),
  or a single list.

- run_ids:

  Optional run names (provenance and collision suffixes); defaults to
  `run1`, `run2`, ...

- on_collision:

  Duplicate sample names across runs: `"suffix"` (append `.<run_id>`,
  default), `"error"`, or `"keep_first"`.

- numeric_warn_frac:

  Warn if a run has more than this fraction of its abundance in
  numeric-only sequences that may not align (default 0.2).

- n_lead:

  Leading accessions used to build the profile signature (default 3).
  Higher = stricter matching (more distinct profiles); lower = looser
  (more merging across runs).

- merge_profiles:

  If `TRUE` (default) retain and cross-match profiles by signature; if
  `FALSE` drop profiles entirely (`prof_mat = NULL`).

## Value

A merged `symbayes` object with a `run` column in metadata,
`filtered = FALSE`, and (if `merge_profiles`) a `prof_mat`/`prof_meta`
whose columns are signature-matched across runs. `prof_meta` gains a
`signature` column and a `runs` column listing which runs each profile
appeared in.

## Details

Why signature-matching works, and its limits: a SymPortal profile's
identity lives in its majority / co-dominant sequences, whose database
accessions are stable across runs. Its trailing minor DIVs are
contextual – appended only if they co-occurred enough in that run – so
the full profile *name* can differ between runs for the same community.
Matching on the leading accessions captures the stable identity while
tolerating the variable tail. The tail variability is itself a
manifestation of SymPortal's core issue (co-occurrence-based,
context-dependent profile definition); signature matching sidesteps it.

## Examples

``` r
if (FALSE) { # \dontrun{
spA <- import("runA.seqs.absolute.abund_only.txt",
              profs_abund = "runA.profiles.absolute.abund_only.txt",
              profs_meta  = "runA.profiles.meta_only.txt")
spB <- import("runB.seqs.absolute.abund_only.txt",
              profs_abund = "runB.profiles.absolute.abund_only.txt",
              profs_meta  = "runB.profiles.meta_only.txt")
sp  <- combine_symportal(spA, spB, run_ids = c("2021", "2023"))
sp$prof_meta[, c("ITS2 type profile", "signature", "runs")]
} # }
```

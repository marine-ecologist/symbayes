# Stable cross-run signature for a SymPortal profile

SymPortal profile UIDs are run-specific, but the defining-sequence
*accessions* are database-level identifiers stable across runs. The
leading accessions (the majority / co-dominant sequences) are the
intrinsic identity; the trailing minor DIVs are run-context-dependent (a
minor sequence is appended only if it co-occurred often enough in that
run's samples). This builds a signature from the first `n_lead`
accessions, so the same biological profile matches across runs even when
its minor-DIV tail differs.

## Usage

``` r
.profile_signature(accession_str, n_lead = 3)
```

## Arguments

- accession_str:

  The `Sequence accession / SymPortal UID` string, e.g.
  `"62825/742546-32-1147953"`.

- n_lead:

  Number of leading accessions to use (default 3).

## Value

A signature string (leading accessions, sorted, joined by `_`).

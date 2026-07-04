# Build notes

The `man/*.Rd` files and the exact `NAMESPACE` are generated from the
roxygen comments in `R/*.R`. This package ships with a hand-maintained
`NAMESPACE` mirror, but you should regenerate both before building:

``` r
# from the package root
install.packages(c("devtools", "roxygen2"))
devtools::document()      # writes man/*.Rd and NAMESPACE from roxygen
```

## Dependencies

Bioconductor (not on CRAN):

``` r
install.packages("BiocManager")
BiocManager::install(c("DirichletMultinomial", "Biostrings"))
```

GitHub (for HDP, optional):

``` r
remotes::install_github("nicolaroberts/hdp")
```

CRAN imports install automatically with the package.

## Check and install

``` r
devtools::check()         # R CMD check
devtools::install()       # install locally
devtools::build_vignettes()
```

## Notes on `.data` and R CMD check

Generalist functions reference metadata columns dynamically via
`.data[[col]]` (imported from `rlang`). This is the tidy-eval idiom that
keeps `R CMD check` free of “no visible binding” NOTEs for column names.

## Test suite

``` r
devtools::test()
```

The tests in `tests/testthat/` use a small synthetic count matrix, so
the schema and DMM smoke tests run without any SymPortal files or the
`hdp` package. HDP-dependent paths are skipped when `hdp` is
unavailable.

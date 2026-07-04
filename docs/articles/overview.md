# 1. SymBayes Overview

## What `symbayes` does

SymPortal assigns ITS2 amplicon data to *type profiles* using
deterministic co-occurrence rules. This has two well-known failure
modes: it **over-splits** (creating many near-identical profiles that
differ only in minor intragenomic variants), and it **mistakes mixed
symbiont communities for unique types** (giving an admixed sample a
single compound profile name).

`symbayes` reframes profiling as probabilistic mixture modelling and
provides three complementary models, each answering a different
question:

| Model | Question | Assignment | K |
|----|----|----|----|
| **DMM** | How many discrete community *types*? | hard | selected (Laplace) |
| **LDA** | What is each sample *made of*? | fractional | selected (redundancy) |
| **HDP** | Both, threshold-free | fractional | **inferred** |

## The `sp` object

Everything flows through a single object created by the
[`import()`](https://marine-ecologist.github.io/symbayes/reference/import.md)
function which imports sections of the Symportal output:

``` r
library(symbayes)

sp <- import(
  seqs_abund  = "run.seqs.absolute.abund_only.txt",
  profs_abund = "run.profiles.absolute.abund_only.txt",
  profs_meta  = "run.profiles.meta_only.txt",
  colour_dict = "color_dict_post_med.json"
)
sp <- filter_samples(sp, min_reads = 100, min_prev = 2)
```

Once the data is imported into `sp`, each `fit_*()` adds a model under
`sp$dmm`, `sp$lda`, or `sp$hdp`, all sharing one schema, so the
generalist plots work on any of them.

``` r
sp <- fit_dmm(sp)          # sp$dmm
sp <- fit_lda(sp)          # sp$lda  (adaptive K)
sp <- fit_hdp(sp)          # sp$hdp  (threshold-free)

sp                          # print method summarises fitted models
```

## Two tiers of function

**Generalist** functions take a `model=` argument and work on any fitted
model:

``` r
plot_pcoa(sp, model = "hdp")
plot_barplot(sp, model = "lda")
plot_exemplars(sp, model = "dmm")
plot_membership(sp, model = "lda")
plot_contingency(sp, model = "hdp")
plot_top_seqs(sp, model = "dmm")
```

**Specialist** functions are prefixed by model and cover diagnostics
unique to that model. See the per-model vignettes:

- [`vignette("dmm")`](https://marine-ecologist.github.io/symbayes/articles/dmm.md)
  — model selection, prediction, LOO, simulation
- [`vignette("lda")`](https://marine-ecologist.github.io/symbayes/articles/lda.md)
  — K search, topic similarity, topic summary
- [`vignette("hdp")`](https://marine-ecologist.github.io/symbayes/articles/hdp.md)
  — inferred component number, over-splitting quantification

## The method argument

Once all three are fitted, compare them directly:

``` r
compare_models(sp)     # group counts + selection strategy + mixture-awareness
compare_mixing(sp)     # per-sample mixing agreement (LDA/HDP vs SymPortal)
```

## Model background

**Dirichlet-multinomial mixture models (DMMs)** treat each sample as
belonging to one discrete community type. The model assumes that samples
are generated from one of K underlying multinomial profiles, where each
profile represents a characteristic assemblage or community state. In
this sense, DMMs are useful when the goal is to classify samples into
relatively clear groups, such as dominant Symbiodiniaceae community
types. Assignment is usually hard, meaning each sample is allocated to a
single best-fitting component, although assignment uncertainty can still
be examined. The number of components, K, is typically selected by
comparing models using criteria such as the Laplace approximation.

**Latent Dirichlet allocation (LDA)** treats each sample as a mixture of
latent topics rather than assigning it to a single community type. Each
topic represents a recurring community profile, and each sample is
described by fractional membership across those topics. This makes LDA
useful when communities are not cleanly separated, or when samples
contain combinations of multiple underlying assemblages. Instead of
asking “which group does this sample belong to?”, LDA asks “what is this
sample made of?”. The number of topics, K, is usually chosen by fitting
models across a range of values and selecting a solution that balances
interpretability, redundancy, and model fit.

**Hierarchical Dirichlet processes (HDPs)** extend the mixture-model
framework by allowing the number of latent components to be inferred
from the data rather than fixed in advance. Like LDA, HDPs allow samples
to have fractional membership across multiple latent community profiles,
but they are threshold-free in the sense that the model can adapt the
number of components as needed. This makes HDPs useful when both the
number of community types and the composition of each sample are
uncertain. In practice, HDPs can reveal fine-scale structure, but they
may also infer fewer or more components depending on how much shared
structure exists across samples.

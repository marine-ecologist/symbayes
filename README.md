# symbayes

Threshold-free, mixture-aware profiling of Symbiodiniaceae ITS2 amplicon data —
a probabilistic alternative deterministic type profiles. `symbayes` reframes profiling as
probabilistic mixture modelling to determine intragenomic variants and admixtures using  
Dirichlet-Multinomial Mixture, Latent Dirichlet Allocation, and Hierarchical Dirichlet Process models. 

## Install

```r
# install.packages("remotes")
remotes::install_github("marine-ecologist/symbayes")

# HDP support
remotes::install_github("NickWilliamsSanger/hdp")

# DirichletMultinomial via Bioconductor:
# BiocManager::install("DirichletMultinomial")
```

## Quick start

```r
library(symbayes)

sp <- import(
  seqs_abund  = "run.seqs.absolute.abund_only.txt",
  profs_abund = "run.profiles.absolute.abund_only.txt",
  profs_meta  = "run.profiles.meta_only.txt",
  colour_dict = "color_dict_post_med.json"
)
sp <- filter_samples(sp, min_reads = 100, min_prev = 2)

sp <- fit_dmm(sp)     # how many discrete types?      (hard assignment)
sp <- fit_lda(sp)     # what is each sample made of?   (fractional, adaptive K)
sp <- fit_hdp(sp)     # both, threshold-free           (K inferred)

compare_models(sp)
hdp_vs_symportal(sp)  # over-splitting quantification
```

## The three models

| Model | Question | Assignment | K selection |
|-------|----------|------------|-------------|
| `fit_dmm` | How many discrete community types? | hard | Laplace/BIC |
| `fit_lda` | What is each sample made of? | fractional | cosine redundancy |
| `fit_hdp` | Both, without thresholds | fractional | **inferred** |

## Function tiers

- **Entry points**: `import`, `filter_samples`, `fit_dmm`, `fit_lda`, `fit_hdp`
- **Generalist plots** (`model=` argument): `plot_pcoa`, `plot_barplot`,
  `plot_exemplars`, `plot_membership`, `plot_contingency`, `plot_top_seqs`
- **DMM-specific**: `dmm_plot_selection`, `predict_dmm`, `dmm_loo`,
  `dmm_subsample_cv`, `dmm_simulate`
- **LDA-specific**: `lda_plot_ksearch`, `lda_plot_similarity`,
  `lda_topic_summary`
- **HDP-specific**: `hdp_plot_numcomp`, `hdp_vs_symportal`
- **Comparison**: `compare_models`, `compare_mixing`

## Documentation

```r
vignette("overview",   package = "symbayes")   # the pipeline
vignette("dmm",        package = "symbayes")   # all DMM functions
vignette("lda",        package = "symbayes")   # all LDA functions
vignette("hdp",        package = "symbayes")   # all HDP functions
vignette("generalist", package = "symbayes")   # model= plots side by side
```

## License

MIT

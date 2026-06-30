# bayes-ace-ordinal-supp

Supplementary materials for the manuscript "Bayesian Estimation of the Univariate ACE Model With Ordinal Data."

To ensure reproducibility, after cloning/downloading the repository, run the following in R to activate `renv`:

```r
renv::restore()
```

- Stan (Marginal Likelihood) vs blavaan: Timing Benchmark
    * [Source code](notebooks/benchmark-stan-blavaan.qmd)
    * [Markdown](notebooks/benchmark-stan-blavaan.md)
    * [PDF](notebooks/benchmark-stan-blavaan.pdf)

- ACE Model (Full Data)
    * [Source code](notebooks/ace-full.qmd)
    * [Markdown](notebooks/ace-full.md)
    * [PDF](notebooks/ace-full.pdf)

- ACE Model (Subset Data)
    * [Source code](notebooks/ace-subset.qmd)
    * [Markdown](notebooks/ace-subset.md)
    * [PDF](notebooks/ace-subset.pdf)

- ACE Model (Full Data) with Mx
    * [Source code](notebooks/ace_mx.qmd)
    * [Markdown](notebooks/ace_mx.md)
    * [PDF](notebooks/ace_mx.pdf)

- Simulation Results
    * [Source code](notebooks/sim_joreskog.qmd)
    * [Markdown](notebooks/sim_joreskog.md)
    * [PDF](notebooks/sim_joreskog.pdf)

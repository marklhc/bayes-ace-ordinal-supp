# bayes-ace-ordinal-supp

Supplementary materials for the manuscript "Bayesian Estimation of the Univariate ACE Model With Ordinal Data."

After cloning/downloading, start R in the supplementary repository root and restore
the package environment:

```r
renv::restore()
```

The packaged `.Rprofile` activates `renv` at startup, and `.here` identifies the
project root. Run the commands below from that root without `--vanilla` so the
startup profile is used. CmdStan must also be installed and configured for
`cmdstanr`; `renv` restores R packages, not CmdStan or the C++ toolchain.

After setup, start with the [ACE Model (Subset Data)](notebooks/ace-subset.qmd)
example. It walks through fitting the model on a smaller data set before moving
on to the full-data analysis.

## Other supplementary materials

- ACE Model (Subset Data)
    * [Source code](notebooks/ace-subset.qmd)
    * [Markdown](notebooks/ace-subset.md)
    * [PDF](notebooks/ace-subset.pdf)

- Stan (Marginal Likelihood) vs blavaan: Timing Benchmark
    * [Source code](notebooks/benchmark-stan-blavaan.qmd)
    * [Markdown](notebooks/benchmark-stan-blavaan.md)
    * [PDF](notebooks/benchmark-stan-blavaan.pdf)

- ACE Model (Full Data)
    * [Source code](notebooks/ace-full.qmd)
    * [Markdown](notebooks/ace-full.md)
    * [PDF](notebooks/ace-full.pdf)

- ACE Model (Full Data) with Mx
    * [Source code](notebooks/ace_mx.qmd)
    * [Markdown](notebooks/ace_mx.md)
    * [PDF](notebooks/ace_mx.pdf)

- Simulation Results
    * [Source code](notebooks/sim_joreskog.qmd)
    * [Markdown](notebooks/sim_joreskog.md)
    * [PDF](notebooks/sim_joreskog.pdf)

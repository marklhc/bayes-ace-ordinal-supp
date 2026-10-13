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

- Stan (Marginal Likelihood) vs blavaan: Illustrative Timing Benchmark
    * [Source code](notebooks/benchmark-stan-blavaan.qmd)
    * [Rendered report](notebooks/benchmark-stan-blavaan.html)

### Running the timing benchmark

This is a small illustrative comparison of Stan's marginal-likelihood approach
with data augmentation in *blavaan* for unrestricted MZ and DZ correlations. It
uses one run per method at each of 100, 500, and 1,552 twin pairs; it is not a
repeated performance study or a comparison of the constrained ACE model. Each
fit uses three chains, 500 warmup iterations, and 1,000 sampling iterations per
chain. The reported whole-fit times include warmup and internal package
post-processing, but exclude compilation.

From the project root, render the report with:

```sh
quarto render notebooks/benchmark-stan-blavaan.qmd
```

The notebook saves individual fit results as RDS files and reuses them on later
renders. Delete the relevant `notebooks/bench_stan_blavaan_<100|500|full>_<stan|blavaan>.rds`
file to rerun that fit. Rendering requires both CmdStan and the R package
*blavaan*, in addition to the packages restored by `renv`.

## Other supplementary materials

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

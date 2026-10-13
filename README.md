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

- Stan (Marginal Likelihood) vs blavaan: Corrected Repeated Benchmark
    * [Source code](notebooks/benchmark-stan-blavaan.qmd)
    * [Reusable R helper](notebooks/benchmark-stan-blavaan.R)
    * [Cache verification script](notebooks/verify-benchmark-stan-blavaan.R)
    * [Markdown](notebooks/benchmark-stan-blavaan.md)
    * [PDF](notebooks/benchmark-stan-blavaan.pdf)

### Running and validating the corrected benchmark

The corrected benchmark uses matching first-threshold/log-spacing priors, common
90% equal-tail intervals, and modern diagnostics for the same four structural
parameters. It includes three repetitions per method at each of 100, 500, and
1,552 twin pairs, with three parallel chains, 500 warmup and 1,000 sampling
iterations per chain. Reported times are whole-fit wall-clock times including
warmup, with compilation and explicit reporting measured separately; the
CmdStanR and rstan implementations perform different internal tasks.

The package includes 18 corrected caches named
`notebooks/bench_stan_blavaan_v2_<100|500|full>_<stan|blavaan>_rep<1|2|3>.rds`.
They preserve draws, diagnostics, conditions and run provenance. To validate them
against the **current environment**, run:

```sh
Rscript notebooks/verify-benchmark-stan-blavaan.R
```

Verification checks provenance, numeric summaries, common interval quantiles,
group/prior mappings, diagnostics and cache invalidation. It prepares the model
(compiling or checking its executable as needed), but does not run fitting or
write caches. Successful validation does not imply that every diagnostic passes;
inspect the printed diagnostic outcomes as well.

The distributed caches record the original run environment. Relocation, changed
library or CmdStan paths, compiler/executable changes, hardware/software changes,
or changes to input data, benchmark sources, priors or settings can invalidate
them even after `renv::restore()`. A provenance mismatch is not a reason to treat
the archived timings as measurements of the new environment. Recompute using:

```sh
Rscript -e 'source("notebooks/benchmark-stan-blavaan.R"); runs <- benchmark_execute()'
Rscript notebooks/verify-benchmark-stan-blavaan.R
```

`benchmark_execute()` reuses validated caches, archives invalid corrected caches
with an `.invalid-...` suffix, and fits missing or invalid runs. Failures are saved
as failures and stop execution; diagnostic flags remain visible. Recomputing the
full repeated benchmark can take tens of minutes. To generate the report, run:

```sh
quarto render notebooks/benchmark-stan-blavaan.qmd
```

Rendering also invokes cache validation and recomputes invalid runs. The corrected
workflow does not consume historical `bench_blavaan_*.rds` summaries.

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

# Corrected, versioned benchmark. Source from the project root (renv activated).
benchmark_settings <- function() {
    list(repetitions = 3L, chains = 3L, warmup = 500L, sampling = 1000L,
         seeds = c(42L, 1042L, 2042L), adapt_delta = 0.95,
         max_treedepth = 12L, cmdstan_sig_figs = 18L, interval_probs = c(0.05, 0.95),
         rhat_limit = 1.01, ess_limit = 400, ebfmi_limit = 0.3)
}

benchmark_model <- function() {
    paste("twin1 ~~ twin2", "twin1 ~~ 1*twin1", "twin2 ~~ 1*twin2",
          "twin1 ~ 0*1", "twin2 ~ 0*1",
          "twin1 | th1*t1 + th2*t2", "twin2 | th1*t1 + th2*t2", sep = "\n")
}

benchmark_data <- function() {
    make_group <- function(counts, group) {
        data.frame(zygosity = group,
                   twin1 = ordered(rep(rep(1:3, 3), counts), levels = 1:3),
                   twin2 = ordered(rep(rep(1:3, each = 3), counts), levels = 1:3))
    }
    full <- rbind(make_group(c(325, 60, 7, 78, 150, 40, 8, 38, 88), "MZ"),
                  make_group(c(202, 97, 19, 96, 112, 64, 46, 56, 66), "DZ"))
    set.seed(1234)
    small <- full[sample.int(nrow(full), 100), ]
    medium <- full[sample.int(nrow(full), 500), ]
    list(`100` = small, `500` = medium, full = full)
}

benchmark_hash <- function(x) digest::digest(x, algo = "sha256")
benchmark_file_hash <- function(path) digest::digest(file = path, algo = "sha256")

benchmark_prepare <- function(root = here::here(), settings = benchmark_settings()) {
    stopifnot(settings$repetitions == length(settings$seeds), settings$chains >= 2L,
              identical(settings$interval_probs, c(0.05, 0.95)))
    packages <- c("cmdstanr", "blavaan", "lavaan", "rstan", "StanHeaders",
                  "posterior", "digest", "renv", "Rcpp", "RcppParallel")
    start <- proc.time()[[3]]
    invisible(lapply(packages, loadNamespace))
    options(mc.cores = settings$chains)
    # Hold post-processing concurrency fixed too; MCMC concurrency is explicit below.
    future::plan(future::sequential)
    package_load_s <- proc.time()[[3]] - start
    start <- proc.time()[[3]]
    model <- cmdstanr::cmdstan_model(file.path(root, "stan", "ace_ordinal_symmetry_for_compare.stan"))
    compilation_check_s <- proc.time()[[3]] - start
    blav_model <- get("stanmodels", asNamespace("blavaan"))$stanmarg
    code <- blav_model@model_code
    # Refuse an unreviewed change to the ordinal prior parameterization.
    stopifnot(grepl("normal_lpdf(Tau_ufree", code, fixed = TRUE),
              grepl("+ exp(Tau_un", code, fixed = TRUE))
    versions <- setNames(lapply(packages, function(p) as.character(utils::packageVersion(p))), packages)
    descriptions <- lapply(packages, function(p) {
        d <- utils::packageDescription(p)
        as.list(d[intersect(names(d), c("Version", "Built", "RemoteSha", "RemoteRef"))])
    })
    names(descriptions) <- packages
    command <- function(cmd, args) paste(system2(cmd, args, stdout = TRUE, stderr = TRUE), collapse = "\n")
    cpu <- if (file.exists("/proc/cpuinfo")) {
        unique(grep("^model name", readLines("/proc/cpuinfo"), value = TRUE))
    } else "unavailable"
    memory <- if (file.exists("/proc/meminfo")) readLines("/proc/meminfo")[1] else "unavailable"
    source_paths <- c("notebooks/benchmark-stan-blavaan.qmd",
                      "notebooks/benchmark-stan-blavaan.R",
                      "stan/ace_ordinal_symmetry_for_compare.stan", "stan/ace_fun.stan", "renv.lock")
    source_hashes <- setNames(lapply(file.path(root, source_paths), benchmark_file_hash), source_paths)
    ns <- asNamespace("blavaan")
    function_hashes <- setNames(lapply(c("bcfa", "blavaan", "lav2stanmarg", "format_priors"),
                                     function(n) benchmark_hash(deparse(get(n, ns)))),
                               c("bcfa", "blavaan", "lav2stanmarg", "format_priors"))
    dll <- list.files(file.path(find.package("blavaan"), "libs"), full.names = TRUE)
    identity <- list(
        schema = "polycor-benchmark-v2", settings = settings,
        model = benchmark_model(), group_order = c("MZ", "DZ"),
        structural_parameters = c("rho_MZ", "rho_DZ", "tau1", "tau2"),
        effective_priors = list(rho = "independent Uniform(-1,1); blavaan LKJ(1) in dimension 2",
                                thresholds = "u1,u2 independent normal(0,1.5); tau1=u1; tau2=tau1+exp(u2)",
                                identification = "means zero; variances one; actual thresholds shared across twins/groups",
                                blavaan_augmentation_bounds = c(-30, 30)),
        interval = list(type = "equal-tail", probability = 0.90, probs = settings$interval_probs),
        timing = paste("External elapsed time from before method-specific input preparation to return of",
                       "sample()/bcfa(); warmup, initialization and internal package post-processing included;",
                       "explicit draw extraction, common summaries/diagnostics and cache I/O excluded;",
                       "libraries loaded and models precompiled before timer; backends differ."),
        source_hashes = source_hashes, blavaan_model_hash = benchmark_hash(code),
        blavaan_function_hashes = function_hashes,
        blavaan_dll_hashes = setNames(lapply(dll, benchmark_file_hash), basename(dll)),
        cmdstan_executable_hash = benchmark_file_hash(model$exe_file()),
        software = list(R = R.version, packages = versions, descriptions = descriptions,
                        cmdstan = as.character(cmdstanr::cmdstan_version()),
                        cmdstan_path = cmdstanr::cmdstan_path(), libraries = .libPaths()),
        compiler = list(cxx = command("g++", "--version"), make = command("make", "--version"),
                        R_CXX17 = command(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17")),
                        R_CXX17FLAGS = command(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17FLAGS"))),
        environment = list(system = Sys.info(), cpu = cpu, memory = memory,
                           detected_cores = parallel::detectCores(),
                           variables = Sys.getenv(c("STAN_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                    "MKL_NUM_THREADS", "CXX", "CXXFLAGS")),
                           future_plan = "sequential", parallel_chains = settings$chains))
    list(root = root, model = model, identity = identity,
         setup = list(package_load_s = package_load_s,
                      cmdstan_compilation_or_freshness_check_s = compilation_check_s,
                      blavaan_compilation = "Precompiled at package build; no per-fit compilation",
                      session = utils::sessionInfo(),
                      recorded_at = format(Sys.time(), tz = "UTC", usetz = TRUE)))
}

benchmark_expected <- function(context, data, repetition, method) {
    list(identity = context$identity, input_hash = benchmark_hash(data),
         input_rows = rownames(data), N = nrow(data),
         group_counts = table(factor(data$zygosity, levels = c("MZ", "DZ"))),
         repetition = repetition, method = method,
         seed = context$identity$settings$seeds[repetition])
}

benchmark_cache_valid <- function(cache, expected) {
    is.list(cache) && identical(cache$schema, "polycor-benchmark-v2") &&
        identical(cache$fingerprint, benchmark_hash(expected)) &&
        identical(cache$provenance, expected) && identical(cache$status, "ok") &&
        isTRUE(tryCatch({ benchmark_verify_run(cache); TRUE }, error = function(e) FALSE))
}

benchmark_blavaan_draws <- function(fit) {
    pt <- as.data.frame(fit@ParTable)
    labels <- fit@Data@group.label
    stopifnot(identical(labels, c("MZ", "DZ")))
    variances <- subset(pt, op == "~~" & lhs == rhs & lhs %in% c("twin1", "twin2"))
    means <- subset(pt, op == "~1" & lhs %in% c("twin1", "twin2"))
    stopifnot(nrow(variances) == 4, all(variances$free == 0), all(variances$est == 1),
              nrow(means) == 4, all(means$free == 0), all(means$est == 0))
    correlations <- subset(pt, op == "~~" & lhs != rhs &
                               lhs %in% c("twin1", "twin2") & rhs %in% c("twin1", "twin2"))
    stopifnot(nrow(correlations) == 2, !anyDuplicated(correlations$group),
              all(correlations$prior == "lkj_corr(1)"))
    source_names <- correlations$pxnames[match(c("MZ", "DZ"), labels[correlations$group])]
    thresholds <- subset(pt, op == "|" & lhs %in% c("twin1", "twin2"))
    for (j in 1:2) {
        rows <- subset(thresholds, rhs == paste0("t", j))
        stopifnot(nrow(rows) == 4, length(unique(rows$pxnames)) == 1,
                  sum(rows$prior != "") == 1, rows$prior[rows$prior != ""] == "normal(0,1.5)")
        source_names <- c(source_names, unique(rows$pxnames))
    }
    d <- fit@external$mcmcdata
    stopifnot(d$len_tau == 2, identical(as.numeric(d$tau_mn), c(0, 0)),
              identical(as.numeric(d$tau_sd), c(1.5, 1.5)),
              length(d$thetablkpri) == 2, all(d$thetablkpri == 1), d$Noent == 2 * d$Ntot)
    raw <- as.array(fit@external$mcmcout)
    stopifnot(!anyNA(source_names), !anyDuplicated(source_names), all(source_names %in% dimnames(raw)[[3]]))
    draws <- raw[, , source_names, drop = FALSE]
    dimnames(draws)[[3]] <- c("rho_MZ", "rho_DZ", "tau1", "tau2")
    list(draws = posterior::as_draws_array(draws),
         mapping = data.frame(parameter = dimnames(draws)[[3]], stan_name = source_names),
         parameter_table = pt,
         effective_prior_data = d[c("tau_mn", "tau_sd", "len_tau", "thetablkpri", "thetanblk", "Noent")])
}

benchmark_summary <- function(draws) {
    as.data.frame(posterior::summarise_draws(
        draws, mean = mean, sd = sd, q05 = function(x) unname(stats::quantile(x, 0.05)),
        q95 = function(x) unname(stats::quantile(x, 0.95)), rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail,
        mcse_mean = posterior::mcse_mean, mcse_sd = posterior::mcse_sd))
}

benchmark_sampler_diagnostics <- function(sampler, settings) {
    # All arrays below are post-warmup. Retain them, not just a success flag.
    sampler <- as.array(sampler)
    stopifnot(all(c("divergent__", "treedepth__", "energy__") %in% dimnames(sampler)[[3]]))
    chains <- lapply(seq_len(dim(sampler)[2]), function(j) {
        energy <- sampler[, j, "energy__"]
        data.frame(chain = j, divergences = sum(sampler[, j, "divergent__"]),
                   treedepth_hits = sum(sampler[, j, "treedepth__"] >= settings$max_treedepth),
                   ebfmi = mean(diff(energy)^2) / stats::var(energy))
    })
    list(scope = "post-warmup sampler states, all sampled dimensions (including augmentation)",
         chains = do.call(rbind, chains), states = sampler)
}

benchmark_run <- function(context, data, repetition, method) {
    stopifnot(method %in% c("stan", "blavaan"))
    # bcfa() rewrites its call to an unqualified blavaan(); provide its binding
    # without depending on which packages a notebook has attached.
    blavaan <- blavaan::blavaan
    expected <- benchmark_expected(context, data, repetition, method)
    s <- context$identity$settings
    warnings <- messages <- character()
    output <- character()
    started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    # Timer encloses input preparation and the full fitting call for BOTH methods.
    start <- proc.time()[[3]]
    fit <- NULL
    error <- NULL
    output <- capture.output(withCallingHandlers(tryCatch({
        if (method == "stan") {
            lower_counts <- function(group) {
                tab <- with(data[data$zygosity == group, ], table(twin1, twin2))
                symm <- tab + t(tab) - diag(diag(tab))
                as.integer(symm[lower.tri(symm, diag = TRUE)])
            }
            fit <- context$model$sample(
                data = list(k = 3L, ltri_ns1 = lower_counts("MZ"), ltri_ns2 = lower_counts("DZ")),
                chains = s$chains, parallel_chains = s$chains, iter_warmup = s$warmup,
                iter_sampling = s$sampling, seed = expected$seed, refresh = 0,
                sig_figs = s$cmdstan_sig_figs,
                adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth)
        } else {
            fit <- blavaan::bcfa(
                benchmark_model(), data = data, ordered = c("twin1", "twin2"),
                group = "zygosity", group.label = c("MZ", "DZ"), target = "stan",
                dp = blavaan::dpriors(tau = "normal(0,1.5)", rho = "beta(1,1)"),
                n.chains = s$chains, burnin = s$warmup, sample = s$sampling,
                seed = expected$seed, test = "none",
                bcontrol = list(cores = s$chains, refresh = 0,
                                control = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth)))
        }
    }, error = function(e) { error <<- conditionMessage(e) }),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }))
    elapsed <- proc.time()[[3]] - start
    result <- list(schema = "polycor-benchmark-v2", fingerprint = benchmark_hash(expected),
                   provenance = expected, setup = context$setup, started_at = started_at,
                   status = if (is.null(error)) "ok" else "error", error = error,
                   elapsed_s = unname(elapsed), conditions = list(warnings = warnings, messages = messages, stdout = output))
    if (!is.null(error)) return(result)
    post_start <- proc.time()[[3]]
    # Do not convert extraction/validation failures into apparently successful fits.
    extracted <- tryCatch({
        if (method == "stan") {
            draws <- fit$draws(c("rho1", "rho2", "thres"), format = "draws_array")
            names_map <- c(rho1 = "rho_MZ", rho2 = "rho_DZ", `thres[1]` = "tau1", `thres[2]` = "tau2")
            source_names <- names(names_map)
            draws <- draws[, , source_names, drop = FALSE]
            dimnames(draws)[[3]] <- unname(names_map)
            result$mapping <- data.frame(parameter = unname(names_map), stan_name = source_names)
            sampler <- fit$sampler_diagnostics(format = "draws_array")
            result$backend_timing <- fit$time()
            result$actual_sampler_settings <- fit$metadata()
            # Check the implemented threshold transform in actual draws.
            u <- as.array(fit$draws("thres_un", format = "draws_array"))
            stopifnot(max(abs(draws[, , "tau1"] - u[, , 1])) < 1e-10,
                      max(abs(draws[, , "tau2"] - draws[, , "tau1"] - exp(u[, , 2]))) < 1e-8)
        } else {
            mapped <- benchmark_blavaan_draws(fit)
            draws <- mapped$draws
            result$mapping <- mapped$mapping
            result$parameter_table <- mapped$parameter_table
            result$effective_prior_data <- mapped$effective_prior_data
            sampler <- simplify2array(rstan::get_sampler_params(fit@external$mcmcout, inc_warmup = FALSE))
            sampler <- aperm(sampler, c(1, 3, 2))
            dimnames(sampler)[[3]] <- colnames(rstan::get_sampler_params(fit@external$mcmcout, inc_warmup = FALSE)[[1]])
            result$backend_timing <- fit@timing
            result$actual_sampler_settings <- fit@external$mcmcout@stan_args
        }
        result$draws <- posterior::as_draws_array(draws)
        result$summary <- benchmark_summary(result$draws)
        result$sampler <- benchmark_sampler_diagnostics(sampler, s)
        z <- result$summary
        c <- result$sampler$chains
        result$diagnostics <- list(
            structural_scope = c("rho_MZ", "rho_DZ", "tau1", "tau2"),
            rhat_max = max(z$rhat), ess_bulk_min = min(z$ess_bulk), ess_tail_min = min(z$ess_tail),
            divergences = sum(c$divergences), treedepth_hits = sum(c$treedepth_hits), ebfmi_min = min(c$ebfmi),
            passed = all(is.finite(z$rhat)) && all(z$rhat <= s$rhat_limit) &&
                all(z$ess_bulk >= s$ess_limit) && all(z$ess_tail >= s$ess_limit) &&
                all(c$divergences == 0) && all(c$treedepth_hits == 0) && all(c$ebfmi >= s$ebfmi_limit))
        benchmark_verify_run(result)
        result
    }, error = function(e) {
        result$status <- "extraction_error"
        result$error <- conditionMessage(e)
        result
    })
    extracted$postprocessing_s <- unname(proc.time()[[3]] - post_start)
    extracted
}

benchmark_verify_run <- function(run) {
    stopifnot(identical(run$status, "ok"))
    s <- run$provenance$identity$settings
    d <- as.array(run$draws)
    stopifnot(identical(dim(d), c(s$sampling, s$chains, 4L)), all(is.finite(d)),
              identical(dimnames(d)[[3]], c("rho_MZ", "rho_DZ", "tau1", "tau2")),
              all(abs(d[, , 1:2]) < 1), all(d[, , "tau2"] > d[, , "tau1"]),
              is.finite(run$elapsed_s), run$elapsed_s > 0,
              identical(run$summary$variable, dimnames(d)[[3]]),
              isTRUE(all.equal(run$summary, benchmark_summary(run$draws), tolerance = 1e-12)),
              dim(run$sampler$states)[1] == s$sampling, dim(run$sampler$states)[2] == s$chains)
    if (run$provenance$method == "blavaan") {
        args <- run$actual_sampler_settings
        stopifnot(length(args) == s$chains,
                  all(vapply(args, function(x) x$warmup == s$warmup &&
                                 x$iter == s$warmup + s$sampling &&
                                 x$seed == run$provenance$seed &&
                                 x$control$adapt_delta == s$adapt_delta &&
                                 x$control$max_treedepth == s$max_treedepth, logical(1))))
    } else {
        args <- run$actual_sampler_settings
        stopifnot(args$iter_warmup == s$warmup, args$iter_sampling == s$sampling,
                  length(args$id) == s$chains, args$seed == run$provenance$seed,
                  args$adapt_delta == s$adapt_delta, args$max_treedepth == s$max_treedepth,
                  args$sig_figs == s$cmdstan_sig_figs)
    }
    invisible(TRUE)
}

benchmark_execute <- function(root = here::here(), settings = benchmark_settings(),
                              dataset_names = c("100", "500", "full")) {
    context <- benchmark_prepare(root, settings)
    datasets <- benchmark_data()
    stopifnot(length(dataset_names) > 0, !anyDuplicated(dataset_names), all(dataset_names %in% names(datasets)))
    datasets <- datasets[dataset_names]
    runs <- list()
    for (name in names(datasets)) {
        for (r in seq_len(settings$repetitions)) {
            # Alternate method order to reduce systematic first/second effects.
            methods <- if (r %% 2L) c("stan", "blavaan") else c("blavaan", "stan")
            for (method in methods) {
                data <- datasets[[name]]
                expected <- benchmark_expected(context, data, r, method)
                path <- file.path(root, "notebooks", sprintf("bench_stan_blavaan_v2_%s_%s_rep%d.rds", name, method, r))
                cached <- if (file.exists(path)) readRDS(path) else NULL
                if (benchmark_cache_valid(cached, expected)) {
                    message("Validated cache: ", basename(path))
                    run <- cached
                } else {
                    if (file.exists(path)) {
                        # Archive invalid provenance/results rather than overwriting evidence.
                        archive <- paste0(path, ".invalid-", format(Sys.time(), "%Y%m%dT%H%M%S"), "-", substr(benchmark_hash(cached), 1, 12))
                        stopifnot(file.rename(path, archive))
                    }
                    message("Running N=", nrow(data), ", ", method, ", repetition ", r)
                    run <- benchmark_run(context, data, r, method)
                    saveRDS(run, path)
                    if (!identical(run$status, "ok")) stop("Benchmark failed; evidence saved in ", path, ": ", run$error)
                    message(sprintf("Completed %.2f s; Rhat %.4f; min bulk/tail ESS %.0f/%.0f; diagnostics %s",
                                    run$elapsed_s, run$diagnostics$rhat_max,
                                    run$diagnostics$ess_bulk_min, run$diagnostics$ess_tail_min,
                                    if (run$diagnostics$passed) "PASS" else "FLAGGED"))
                }
                runs[[paste(name, method, r, sep = "_")]] <- run
            }
        }
    }
    runs
}

benchmark_tables <- function(runs) {
    timing <- estimates <- diagnostics <- list()
    for (i in seq_along(runs)) {
        run <- runs[[i]]
        benchmark_verify_run(run)
        p <- run$provenance
        method <- if (p$method == "stan") "Stan (marginal / CmdStanR)" else "blavaan (augmentation / rstan)"
        timing[[i]] <- data.frame(N = p$N, Method = method, repetition = p$repetition,
                                 elapsed_s = run$elapsed_s, postprocessing_s = run$postprocessing_s)
        z <- run$summary
        z$N <- p$N; z$Method <- method; z$repetition <- p$repetition
        z$ess_bulk_per_s <- z$ess_bulk / run$elapsed_s
        z$ess_tail_per_s <- z$ess_tail / run$elapsed_s
        estimates[[i]] <- z
        diagnostics[[i]] <- cbind(timing[[i]][, 1:3], as.data.frame(run$diagnostics[setdiff(names(run$diagnostics), "structural_scope")]),
                                  warning_count = length(run$conditions$warnings))
    }
    list(timing = do.call(rbind, timing), estimates = do.call(rbind, estimates),
         diagnostics = do.call(rbind, diagnostics))
}

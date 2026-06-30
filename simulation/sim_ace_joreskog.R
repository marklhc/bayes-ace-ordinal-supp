# sim_ace_joreskog.R
#
# Simulation comparing umx (OpenMx, frequentist) vs Stan (Bayesian) for
# estimating h² and c² from an ordinal ACE twin model.
#
# Population: Jöreskog BMI data (3 categories: normal / overweight / obese)
# Data generation: resample from empirical cell proportions
# Metrics: bias, RMSE, CI coverage, CI width, convergence rate,
#          admissible solution rate

library(here)
library(SimDesign)
library(cmdstanr)
library(R.utils) # for withTimeout()
library(umx)
library(posterior)

# ── Population data (Jöreskog BMI contingency tables) ─────────────────────────
# Rows = twin1, Cols = twin2; categories: 1=normal, 2=overweight, 3=obese
N_MZ_POP <- matrix(
    c(325, 60, 7, 78, 150, 40, 8, 38, 88),
    nrow = 3,
    byrow = TRUE
)
N_DZ_POP <- matrix(
    c(202, 97, 19, 96, 112, 64, 46, 56, 66),
    nrow = 3,
    byrow = TRUE
)
P_MZ <- N_MZ_POP / sum(N_MZ_POP)
P_DZ <- N_DZ_POP / sum(N_DZ_POP)

# Symmetrise: average upper and lower triangles so P[i,j] == P[j,i]
symmetrise <- function(p) (p + t(p)) / 2
P_MZ <- symmetrise(P_MZ)
P_DZ <- symmetrise(P_DZ)

K <- 3L
LABELS <- c("normal", "overweight", "obese")

# ── True parameter values ─────────────────────────────────────────────────────
# Polychoric MZ/DZ correlations → Falconer's formulas for h² and c²
r_mz_true <- polycor::polychor(N_MZ_POP)
r_dz_true <- polycor::polychor(N_DZ_POP)
H2_TRUE <- 2 * (r_mz_true - r_dz_true) # additive genetic
C2_TRUE <- 2 * r_dz_true - r_mz_true # shared environment
message(sprintf(
    "Population: r_MZ=%.3f, r_DZ=%.3f  →  h²=%.3f, c²=%.3f",
    r_mz_true,
    r_dz_true,
    H2_TRUE,
    C2_TRUE
))

# ── Helpers ───────────────────────────────────────────────────────────────────
# Symmetrised lower-triangular cell counts (input to the Stan model)
to_ltri <- function(tab) {
    symm <- tab + t(tab) - diag(diag(tab))
    symm[lower.tri(symm, diag = TRUE)]
}

# Expand integer count vector to long-format ordered data frame
counts_to_df <- function(counts, k = K, labels = LABELS) {
    data.frame(
        twin1 = ordered(rep(rep(seq_len(k), k), counts), labels = labels),
        twin2 = ordered(rep(rep(seq_len(k), each = k), counts), labels = labels)
    )
}

# ── Fixed simulation objects ──────────────────────────────────────────────────
FIXED <- list(
    p_mz = P_MZ,
    p_dz = P_DZ,
    h2_true = H2_TRUE,
    c2_true = C2_TRUE,
    ace_mod = cmdstan_model(
        here("stan", "ace_ordinal_symmetry_simplex_no_x2.stan"),
        include_paths = here("stan"),
        quiet = TRUE
    ),
    ace_reg_mod = cmdstan_model(
        here("stan", "ace_ordinal_symmetry_simplex_reg.stan"),
        include_paths = here("stan"),
        quiet = TRUE
    ),
    ace_nc_mod = cmdstan_model(
        here("stan", "ace_ordinal_symmetry_no_x2.stan"),
        include_paths = here("stan"),
        quiet = TRUE
    )
)

# ── Design conditions ─────────────────────────────────────────────────────────
DESIGN <- createDesign(n_pairs = c(50L, 100L))

# ── Data generation ───────────────────────────────────────────────────────────
gen <- function(condition, fixed_objects) {
    n <- condition$n_pairs
    list(
        dat_mz = counts_to_df(c(rmultinom(
            1L,
            n,
            prob = c(fixed_objects$p_mz)
        ))),
        dat_dz = counts_to_df(c(rmultinom(1L, n, prob = c(fixed_objects$p_dz))))
    )
}

# ── Stan helpers ──────────────────────────────────────────────────────────────
kd_mode <- function(x) {
    d <- density(x, bw = "SJ")
    d$x[which.max(d$y)]
}

.stan_na_return <- c(
    converged = 0L,
    admissible = NA_integer_,
    n_divergent = NA_integer_,
    h2_est = NA_real_,
    h2_est_med = NA_real_,
    h2_est_mode = NA_real_,
    h2_lo = NA_real_,
    h2_hi = NA_real_,
    h2_sd = NA_real_,
    c2_est = NA_real_,
    c2_est_med = NA_real_,
    c2_est_mode = NA_real_,
    c2_lo = NA_real_,
    c2_hi = NA_real_,
    c2_sd = NA_real_
)

# Escalating sampler settings for retries
.stan_retry_schedule <- list(
    list(adapt_delta = 0.90, max_treedepth = 10L),
    list(adapt_delta = 0.95, max_treedepth = 12L),
    list(adapt_delta = 0.99, max_treedepth = 15L)
)

# Worker: fit one Stan model and return the summary vector.
# `simplex = TRUE`  → admissibility guaranteed by parameterisation.
# `simplex = FALSE` → admissibility assessed from posterior means.
# Retries up to length(.stan_retry_schedule) times with escalating sampler
# settings if convergence fails; returns .stan_na_return if all attempts fail.
run_stan_model <- function(stan_mod, stan_data, simplex = TRUE) {
    summarise_fit <- function(fit) {
        diags <- fit$diagnostic_summary(quiet = TRUE)
        s <- fit$summary(c("h2", "c2"), mean, rhat, ess_bulk)
        n_div <- sum(diags$num_divergent)
        rhat_ok <- all(s$rhat < 1.01, na.rm = TRUE)
        ess_ok <- all(s$ess_bulk > 400, na.rm = TRUE)
        list(
            fit = fit,
            n_div = n_div,
            converged = rhat_ok && ess_ok && n_div == 0
        )
    }

    attempt <- NULL
    for (sched in .stan_retry_schedule) {
        attempt <- tryCatch(
            {
                fit <-
                    withTimeout(
                        {
                            stan_mod$sample(
                                data = stan_data,
                                chains = 3L,
                                iter_warmup = 1000L,
                                iter_sampling = 1500L,
                                adapt_delta = sched$adapt_delta,
                                max_treedepth = sched$max_treedepth,
                                parallel_chains = 1L,
                                refresh = 0L,
                                show_messages = FALSE,
                                show_exceptions = FALSE
                            )
                        },
                        timeout = 180
                    )
                summarise_fit(fit)
            },
            TimeoutException = function(ex) {
                message("Timeout. Skipping.")
            },
            error = function(e) NULL
        )
        if (!is.null(attempt) && attempt$converged) break
    }

    if (is.null(attempt)) {
        return(.stan_na_return)
    }

    fit <- attempt$fit
    n_div <- attempt$n_div

    h2_draws <- as.vector(fit$draws("h2", format = "matrix"))
    c2_draws <- as.vector(fit$draws("c2", format = "matrix"))
    h2_mean <- mean(h2_draws)
    c2_mean <- mean(c2_draws)

    admissible <- if (simplex) {
        1L # simplex parameterisation: always in [0,1]
    } else {
        as.integer(
            h2_mean >= -1e-4 &
                c2_mean >= -1e-4 &
                h2_mean + c2_mean <= 1 + 1e-4
        )
    }

    if (!attempt$converged) {
        return(.stan_na_return)
    }

    c(
        converged = 1L,
        admissible = admissible,
        n_divergent = n_div,
        h2_est = h2_mean,
        h2_est_med = median(h2_draws),
        h2_est_mode = kd_mode(h2_draws),
        h2_lo = quantile(h2_draws, 0.05, names = FALSE),
        h2_hi = quantile(h2_draws, 0.95, names = FALSE),
        h2_sd = sd(h2_draws),
        c2_est = c2_mean,
        c2_est_med = median(c2_draws),
        c2_est_mode = kd_mode(c2_draws),
        c2_lo = quantile(c2_draws, 0.05, names = FALSE),
        c2_hi = quantile(c2_draws, 0.95, names = FALSE),
        c2_sd = sd(c2_draws)
    )
}

# Build the contingency-table data list once per replicate.
make_stan_data <- function(dat) {
    list(
        k = K,
        ltri_ns1 = to_ltri(table(dat$dat_mz$twin1, dat$dat_mz$twin2)),
        ltri_ns2 = to_ltri(table(dat$dat_dz$twin1, dat$dat_dz$twin2))
    )
}

# ── Stan analysis (simplex) ───────────────────────────────────────────────────
analyse_stan <- function(condition, dat, fixed_objects) {
    run_stan_model(fixed_objects$ace_mod, make_stan_data(dat), simplex = TRUE)
}

# ── Stan analysis (simplex + Dirichlet regularising prior) ───────────────────
analyse_stan_reg <- function(condition, dat, fixed_objects) {
    run_stan_model(
        fixed_objects$ace_reg_mod,
        make_stan_data(dat),
        simplex = TRUE
    )
}

# ── Stan analysis (non-simplex / unconstrained) ───────────────────────────────
# rho1, rho2 free in (-1,1); h2 and c2 derived in generated quantities and
# can be outside [0,1].
analyse_stan_nc <- function(condition, dat, fixed_objects) {
    run_stan_model(
        fixed_objects$ace_nc_mod,
        make_stan_data(dat),
        simplex = FALSE
    )
}

# ── umx analysis ──────────────────────────────────────────────────────────────
analyse_umx <- function(condition, dat, fixed_objects) {
    tryCatch(
        {
            mz2 <- dat$dat_mz
            dz2 <- dat$dat_dz
            mz2[] <- lapply(mz2, umxFactor)
            dz2[] <- lapply(dz2, umxFactor)

            m <- suppressMessages(suppressWarnings(
                umxACE(
                    selDVs = "twin",
                    mzData = mz2,
                    dzData = dz2,
                    sep = "",
                    intervals = FALSE,
                    autoRun = TRUE
                )
            ))

            a_val <- m$top$a$values[1, 1]
            c_val <- m$top$c$values[1, 1]
            e_val <- m$top$e$values[1, 1]

            # For ordinal ACE, variance = 1 by constraint, so h² = a², c² = c²
            h2 <- a_val^2
            c2 <- c_val^2
            e2 <- e_val^2

            # Convergence: estimates finite + Hessian invertible (SEs available)
            par_df <- tryCatch(
                as.data.frame(summary(m)$parameters),
                error = function(e) NULL
            )

            # Extract SE by matrix label (top.a, top.c) and position (row=1, col=1)
            get_se <- function(df, mat_label) {
                if (is.null(df)) {
                    return(NA_real_)
                }
                idx <- which(df$matrix == mat_label & df$row == 1 & df$col == 1)
                if (length(idx) != 1) {
                    return(NA_real_)
                }
                se <- df$Std.Error[idx]
                if (!is.finite(se) || se <= 0) NA_real_ else se
            }
            se_a <- get_se(par_df, "top.a")
            se_c <- get_se(par_df, "top.c")

            estimates_ok <- all(is.finite(c(h2, c2, e2)))
            ses_ok <- !is.na(se_a) && !is.na(se_c)
            converged <- as.integer(estimates_ok && ses_ok)

            # Admissibility: h2 ≥ 0, c2 ≥ 0, e2 > 0, and h2+c2+e2 ≈ 1
            admissible <- as.integer(
                h2 >= -1e-4 &
                    c2 >= -1e-4 &
                    e2 > 1e-4 &
                    abs(h2 + c2 + e2 - 1) < 0.01
            )

            # 90% profile (likelihood-based) CI on h² and c² directly via MxAlgebra
            ci_bounds <- tryCatch(
                {
                    m_ci <- mxModel(
                        m,
                        mxAlgebra(top.a * top.a, name = "h2_ace"),
                        mxAlgebra(top.c * top.c, name = "c2_ace"),
                        mxCI(c("h2_ace[1,1]", "c2_ace[1,1]"), interval = 0.90)
                    )
                    m_ci <- suppressMessages(suppressWarnings(
                        mxRun(m_ci, intervals = TRUE, silent = TRUE)
                    ))
                    ci_mat <- m_ci$output$confidenceIntervals
                    get_ci <- function(mat, nm) {
                        idx <- which(rownames(mat) == nm)
                        if (length(idx) == 0) {
                            return(c(NA_real_, NA_real_))
                        }
                        c(mat[idx, "lbound"], mat[idx, "ubound"])
                    }
                    list(
                        h2 = get_ci(ci_mat, "h2_ace[1,1]"),
                        c2 = get_ci(ci_mat, "c2_ace[1,1]")
                    )
                },
                error = function(e) {
                    list(
                        h2 = c(NA_real_, NA_real_),
                        c2 = c(NA_real_, NA_real_)
                    )
                }
            )

            profile_ci_ok <- as.integer(all(is.finite(c(
                ci_bounds$h2,
                ci_bounds$c2
            ))))

            # Wald 90% CI via delta method: θ = x² → se(θ) = 2|x| × se(x)
            h2_se <- abs(2 * a_val) * se_a
            c2_se <- abs(2 * c_val) * se_c

            c(
                converged = converged,
                admissible = admissible,
                profile_ci_ok = profile_ci_ok,
                h2_est = h2,
                h2_lo = ci_bounds$h2[1], # profile CI lower
                h2_hi = ci_bounds$h2[2], # profile CI upper
                h2_lo_wald = h2 - 1.645 * h2_se, # Wald CI lower
                h2_hi_wald = h2 + 1.645 * h2_se, # Wald CI upper
                h2_sd = h2_se, # asymptotic SE
                c2_est = c2,
                c2_lo = ci_bounds$c2[1],
                c2_hi = ci_bounds$c2[2],
                c2_lo_wald = c2 - 1.645 * c2_se,
                c2_hi_wald = c2 + 1.645 * c2_se,
                c2_sd = c2_se
            )
        },
        error = function(e) {
            c(
                converged = 0L,
                admissible = NA_integer_,
                profile_ci_ok = NA_integer_,
                h2_est = NA_real_,
                h2_lo = NA_real_,
                h2_hi = NA_real_,
                h2_lo_wald = NA_real_,
                h2_hi_wald = NA_real_,
                h2_sd = NA_real_,
                c2_est = NA_real_,
                c2_lo = NA_real_,
                c2_hi = NA_real_,
                c2_lo_wald = NA_real_,
                c2_hi_wald = NA_real_,
                c2_sd = NA_real_
            )
        }
    )
}

# ── Evaluate ──────────────────────────────────────────────────────────────────
get_cols <- function(results, pattern) {
    results[, grep(pattern, colnames(results)), drop = FALSE]
}

# Strip trailing ".varname" to get clean method names (stan / umx) as col names
clean_names <- function(x, suffix) {
    setNames(as.data.frame(x), sub(paste0("\\.", suffix, "$"), "", colnames(x)))
}

evaluate <- function(condition, results, fixed_objects) {
    h2_true <- fixed_objects$h2_true
    c2_true <- fixed_objects$c2_true

    converged <- get_cols(results, "\\.converged$")
    admissible <- get_cols(results, "\\.admissible$")
    profile_ci_ok <- get_cols(results, "\\.profile_ci_ok$") # umx only

    h2_est <- get_cols(results, "\\.h2_est$")
    h2_est_med <- get_cols(results, "\\.h2_est_med$") # stan only
    h2_est_mode <- get_cols(results, "\\.h2_est_mode$") # stan only
    h2_lo <- get_cols(results, "\\.h2_lo$") # stan: posterior; umx: profile
    h2_hi <- get_cols(results, "\\.h2_hi$")
    h2_lo_wald <- get_cols(results, "\\.h2_lo_wald$") # umx Wald only
    h2_hi_wald <- get_cols(results, "\\.h2_hi_wald$")
    h2_sd <- get_cols(results, "\\.h2_sd$")

    c2_est <- get_cols(results, "\\.c2_est$")
    c2_est_med <- get_cols(results, "\\.c2_est_med$") # stan only
    c2_est_mode <- get_cols(results, "\\.c2_est_mode$") # stan only
    c2_lo <- get_cols(results, "\\.c2_lo$")
    c2_hi <- get_cols(results, "\\.c2_hi$")
    c2_lo_wald <- get_cols(results, "\\.c2_lo_wald$")
    c2_hi_wald <- get_cols(results, "\\.c2_hi_wald$")
    c2_sd <- get_cols(results, "\\.c2_sd$")

    # Helper: compute metrics (coverage uses primary CI: posterior for Stan, profile for umx)
    metrics <- function(est, lo, hi, sd_col, true) {
        err <- est - true
        m <- rbind(
            bias = colMeans(err, na.rm = TRUE),
            rmse = sqrt(colMeans(err^2, na.rm = TRUE)),
            emp_sd = apply(est, 2, sd, na.rm = TRUE),
            avg_se = colMeans(sd_col, na.rm = TRUE),
            coverage = colMeans(lo <= true & true <= hi, na.rm = TRUE),
            ci_width = colMeans(hi - lo, na.rm = TRUE)
        )
        colnames(m) <- sub("\\.[^.]+$", "", colnames(m))
        m
    }

    h2_m <- metrics(h2_est, h2_lo, h2_hi, h2_sd, h2_true)
    c2_m <- metrics(c2_est, c2_lo, c2_hi, c2_sd, c2_true)

    # Point-estimate-only metrics for median and mode (stan only; drop all-NA cols)
    pt_metrics <- function(est_df, true) {
        est_df <- est_df[, colSums(!is.na(est_df)) > 0, drop = FALSE]
        if (ncol(est_df) == 0) {
            return(NULL)
        }
        err <- est_df - true
        nms <- sub("\\.[^.]+$", "", colnames(est_df))
        rbind(
            bias = setNames(colMeans(err, na.rm = TRUE), nms),
            rmse = setNames(sqrt(colMeans(err^2, na.rm = TRUE)), nms),
            emp_sd = setNames(apply(est_df, 2, sd, na.rm = TRUE), nms)
        )
    }
    h2_med_m <- pt_metrics(h2_est_med, h2_true)
    h2_mode_m <- pt_metrics(h2_est_mode, h2_true)
    c2_med_m <- pt_metrics(c2_est_med, c2_true)
    c2_mode_m <- pt_metrics(c2_est_mode, c2_true)

    # Wald CI metrics (umx only)
    clean <- function(x) setNames(x, sub("\\.[^.]+$", "", names(x)))
    h2_cov_wald <- clean(colMeans(
        h2_lo_wald <= h2_true & h2_true <= h2_hi_wald,
        na.rm = TRUE
    ))
    h2_wid_wald <- clean(colMeans(h2_hi_wald - h2_lo_wald, na.rm = TRUE))
    c2_cov_wald <- clean(colMeans(
        c2_lo_wald <= c2_true & c2_true <= c2_hi_wald,
        na.rm = TRUE
    ))
    c2_wid_wald <- clean(colMeans(c2_hi_wald - c2_lo_wald, na.rm = TRUE))

    conv_rate <- clean(colMeans(converged == 1, na.rm = TRUE))
    admis_rate <- clean(colMeans(admissible == 1, na.rm = TRUE))
    prof_ci_rate <- clean(colMeans(profile_ci_ok == 1, na.rm = TRUE))

    # Helper to emit named metric rows for median/mode (NULL-safe)
    emit_pt <- function(m, param, estimator) {
        if (is.null(m)) {
            return(NULL)
        }
        c(
            setNames(
                m["bias", ],
                paste0(param, "_bias_", estimator, ".", names(m["bias", ]))
            ),
            setNames(
                m["rmse", ],
                paste0(param, "_rmse_", estimator, ".", names(m["rmse", ]))
            ),
            setNames(
                m["emp_sd", ],
                paste0(param, "_empsd_", estimator, ".", names(m["emp_sd", ]))
            )
        )
    }

    c(
        # h² metrics (stan: posterior CI; umx: profile CI)
        setNames(h2_m["bias", ], paste0("h2_bias.", names(h2_m["bias", ]))),
        setNames(h2_m["rmse", ], paste0("h2_rmse.", names(h2_m["rmse", ]))),
        setNames(
            h2_m["emp_sd", ],
            paste0("h2_empsd.", names(h2_m["emp_sd", ]))
        ),
        setNames(
            h2_m["avg_se", ],
            paste0("h2_avgse.", names(h2_m["avg_se", ]))
        ),
        setNames(
            h2_m["coverage", ],
            paste0("h2_coverage.", names(h2_m["coverage", ]))
        ),
        setNames(
            h2_m["ci_width", ],
            paste0("h2_ciwidth.", names(h2_m["ci_width", ]))
        ),
        setNames(h2_cov_wald, paste0("h2_coverage_wald.", names(h2_cov_wald))),
        setNames(h2_wid_wald, paste0("h2_ciwidth_wald.", names(h2_wid_wald))),
        # c² metrics
        setNames(c2_m["bias", ], paste0("c2_bias.", names(c2_m["bias", ]))),
        setNames(c2_m["rmse", ], paste0("c2_rmse.", names(c2_m["rmse", ]))),
        setNames(
            c2_m["emp_sd", ],
            paste0("c2_empsd.", names(c2_m["emp_sd", ]))
        ),
        setNames(
            c2_m["avg_se", ],
            paste0("c2_avgse.", names(c2_m["avg_se", ]))
        ),
        setNames(
            c2_m["coverage", ],
            paste0("c2_coverage.", names(c2_m["coverage", ]))
        ),
        setNames(
            c2_m["ci_width", ],
            paste0("c2_ciwidth.", names(c2_m["ci_width", ]))
        ),
        setNames(c2_cov_wald, paste0("c2_coverage_wald.", names(c2_cov_wald))),
        setNames(c2_wid_wald, paste0("c2_ciwidth_wald.", names(c2_wid_wald))),
        # h² median and mode (stan only)
        emit_pt(h2_med_m, "h2", "med"),
        emit_pt(h2_mode_m, "h2", "mode"),
        # c² median and mode (stan only)
        emit_pt(c2_med_m, "c2", "med"),
        emit_pt(c2_mode_m, "c2", "mode"),
        # Rates
        setNames(conv_rate, paste0("converge_rate.", names(conv_rate))),
        setNames(admis_rate, paste0("admissible_rate.", names(admis_rate))),
        setNames(prof_ci_rate, paste0("profile_ci_rate.", names(prof_ci_rate)))
    )
}

# ── Run simulation ─────────────────────────────────────────────────────────────
# Each Stan call uses 1 core (chains run sequentially); parallelism is across
# replicates via runSimulation().

out <- runSimulation(
    design = DESIGN,
    replications = 1000L,
    generate = gen,
    analyse = list(
        stan = analyse_stan,
        stan_reg = analyse_stan_reg,
        stan_nc = analyse_stan_nc,
        umx = analyse_umx
    ),
    summarise = evaluate,
    seed = rep(20250324L, nrow(DESIGN)),
    parallel = TRUE,
    ncores = min(parallel::detectCores() - 1L, 50L),
    packages = c("cmdstanr", "posterior", "umx", "polycor", "R.utils"),
    fixed_objects = FIXED,
    filename = "simulation/sim_ace_joreskog.rds",
    save_results = TRUE,
    save_details = list(
        save_results_dirname = "simulation/sim_ace_joreskog_results_"
    ),
    control = list(
        allow_na = TRUE,
        stop_on_fatal = FALSE
    )
)

# ── Quick summary table ───────────────────────────────────────────────────────
print(out)

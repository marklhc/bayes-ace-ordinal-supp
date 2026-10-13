# sim_ace_joreskog.R
#
# Simulation comparing umx (OpenMx, frequentist) vs Stan (Bayesian) for
# estimating h² and c² from an ordinal ACE twin model.
#
# Population: Jöreskog BMI data (3 categories: normal / overweight / obese)
# Data generation: resample from empirical cell proportions
# Metrics: bias, RMSE, CI coverage, CI width, convergence rate,
#          admissible solution rate
# This is empirical resampling, not generation under a known ACE model.
# Bias/coverage use descriptive polychoric reference values, not known truths.

library(here)
library(SimDesign)
library(cmdstanr)
library(R.utils) # for withTimeout()
library(umx)
library(posterior)
library(HDInterval)

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

# ── Reference parameter values ────────────────────────────────────────────────
# Two-step polychorics of the actual symmetrized generating tables.
# These tables need not satisfy a Gaussian ACE model with common thresholds.
# Keep the existing TRUE names for compatibility with the summary code.
r_mz_true <- polycor::polychor(symmetrise(N_MZ_POP), ML = FALSE)
r_dz_true <- polycor::polychor(symmetrise(N_DZ_POP), ML = FALSE)
H2_TRUE <- 2 * (r_mz_true - r_dz_true) # additive genetic reference
C2_TRUE <- 2 * r_dz_true - r_mz_true # shared environment reference
message(sprintf(
    "Reference: r_MZ=%.3f, r_DZ=%.3f  →  h²=%.3f, c²=%.3f",
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
        twin1 = ordered(rep(rep(seq_len(k), k), counts), levels = seq_len(k), labels = labels),
        twin2 = ordered(rep(rep(seq_len(k), each = k), counts), levels = seq_len(k), labels = labels)
    )
}

# ── Fixed simulation objects ──────────────────────────────────────────────────
# All Stan models use standard-normal threshold priors. Their ACE/correlation
# priors and parameter spaces still differ.
FIXED <- list(
    k = K,
    labels = LABELS,
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
        here("stan", "ace_ordinal_symmetry_no_x2_normal.stan"),
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
        )), k = fixed_objects$k, labels = fixed_objects$labels),
        dat_dz = counts_to_df(c(rmultinom(1L, n, prob = c(fixed_objects$p_dz))),
                             k = fixed_objects$k, labels = fixed_objects$labels)
    )
}

# ── Stan helpers ──────────────────────────────────────────────────────────────
# Per-method wall time in seconds, including preparation, retries, diagnostics
# and requested intervals, even for failed fits. Excludes compilation/generation.
timed_analysis <- function(expr) {
    started <- proc.time()[["elapsed"]]
    result <- force(expr)
    c(result, elapsed_sec = unname(proc.time()[["elapsed"]] - started))
}

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
    h2_lo_hdi = NA_real_,
    h2_hi_hdi = NA_real_,
    h2_sd = NA_real_,
    c2_est = NA_real_,
    c2_est_med = NA_real_,
    c2_est_mode = NA_real_,
    c2_lo = NA_real_,
    c2_hi = NA_real_,
    c2_lo_hdi = NA_real_,
    c2_hi_hdi = NA_real_,
    c2_sd = NA_real_
)

# Escalating sampler settings for retries
.stan_retry_schedule <- list(
    list(adapt_delta = 0.90, max_treedepth = 10L),
    list(adapt_delta = 0.95, max_treedepth = 12L),
    list(adapt_delta = 0.99, max_treedepth = 15L)
)
# Pass settings through fixed_objects so parallel workers do not require globals.
FIXED$stan_na_return <- .stan_na_return
FIXED$stan_retry_schedule <- .stan_retry_schedule

# Worker: fit one Stan model and return the summary vector.
# `simplex = TRUE`  → admissibility guaranteed by parameterisation.
# `simplex = FALSE` → admissibility assessed from posterior means.
# Retries up to length(.stan_retry_schedule) times with escalating sampler
# settings if convergence fails; returns .stan_na_return if all attempts fail.
run_stan_model <- function(stan_mod, stan_data, simplex = TRUE,
                           retry_schedule = .stan_retry_schedule,
                           na_return = .stan_na_return) {
    summarise_fit <- function(fit) {
        diags <- fit$diagnostic_summary(quiet = TRUE)
        s <- fit$summary(c("h2", "c2"), mean, rhat, ess_bulk)
        n_div <- sum(diags$num_divergent)
        rhat_ok <- nrow(s) > 0 && all(is.finite(s$rhat)) && all(s$rhat < 1.01)
        ess_ok <- nrow(s) > 0 && all(is.finite(s$ess_bulk)) && all(s$ess_bulk > 400)
        list(
            fit = fit,
            n_div = n_div,
            converged = rhat_ok && ess_ok && is.finite(n_div) && n_div == 0
        )
    }

    attempt <- NULL
    for (sched in retry_schedule) {
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
        return(na_return)
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
        return(na_return)
    }

    h2_hdi <- HDInterval::hdi(h2_draws, credMass = 0.9) |> unname()
    c2_hdi <- HDInterval::hdi(c2_draws, credMass = 0.9) |> unname()
    c(
        converged = 1L,
        admissible = admissible,
        n_divergent = n_div,
        h2_est = h2_mean,
        h2_est_med = median(h2_draws),
        h2_est_mode = kd_mode(h2_draws),
        h2_lo = quantile(h2_draws, 0.05, names = FALSE),
        h2_hi = quantile(h2_draws, 0.95, names = FALSE),
        h2_lo_hdi = h2_hdi[1],
        h2_hi_hdi = h2_hdi[2],
        h2_sd = sd(h2_draws),
        c2_est = c2_mean,
        c2_est_med = median(c2_draws),
        c2_est_mode = kd_mode(c2_draws),
        c2_lo = quantile(c2_draws, 0.05, names = FALSE),
        c2_hi = quantile(c2_draws, 0.95, names = FALSE),
        c2_lo_hdi = c2_hdi[1],
        c2_hi_hdi = c2_hdi[2],
        c2_sd = sd(c2_draws)
    )
}

# Build the contingency-table data list once per replicate.
make_stan_data <- function(dat) {
    list(
        k = nlevels(dat$dat_mz$twin1),
        ltri_ns1 = to_ltri(table(dat$dat_mz$twin1, dat$dat_mz$twin2)),
        ltri_ns2 = to_ltri(table(dat$dat_dz$twin1, dat$dat_dz$twin2))
    )
}

# ── Stan analysis (simplex) ───────────────────────────────────────────────────
analyse_stan <- function(condition, dat, fixed_objects) {
    timed_analysis(run_stan_model(fixed_objects$ace_mod, make_stan_data(dat), simplex = TRUE,
                   retry_schedule = fixed_objects$stan_retry_schedule,
                   na_return = fixed_objects$stan_na_return))
}

# ── Stan analysis (simplex + Dirichlet regularising prior) ───────────────────
analyse_stan_reg <- function(condition, dat, fixed_objects) {
    timed_analysis(run_stan_model(
        fixed_objects$ace_reg_mod,
        make_stan_data(dat),
        simplex = TRUE,
        retry_schedule = fixed_objects$stan_retry_schedule,
        na_return = fixed_objects$stan_na_return
    ))
}

# ── Stan analysis (non-simplex / unconstrained) ───────────────────────────────
# rho1, rho2 free in (-1,1); h2 and c2 derived in generated quantities and
# can be outside [0,1].
analyse_stan_nc <- function(condition, dat, fixed_objects) {
    timed_analysis(run_stan_model(
        fixed_objects$ace_nc_mod,
        make_stan_data(dat),
        simplex = FALSE,
        retry_schedule = fixed_objects$stan_retry_schedule,
        na_return = fixed_objects$stan_na_return
    ))
}

# OpenMx may prefix the target with the model name. Finite endpoints alone
# do not establish success: both endpoint status codes must also be zero.
get_profile_ci <- function(output, nm) {
    mat <- output$confidenceIntervals
    codes <- output$confidenceIntervalCodes
    find_row <- function(x) {
        if (is.null(x)) return(integer())
        which(rownames(x) == nm | endsWith(rownames(x), paste0(".", nm)))
    }
    idx <- find_row(mat)
    idx_code <- find_row(codes)
    if (length(idx) != 1L || length(idx_code) != 1L) return(c(NA_real_, NA_real_))
    bounds <- as.numeric(mat[idx, c("lbound", "ubound")])
    status <- as.numeric(codes[idx_code, c("lbound", "ubound")])
    if (!all(is.finite(bounds)) || bounds[1] > bounds[2] ||
        !all(is.finite(status)) || any(status != 0)) return(c(NA_real_, NA_real_))
    bounds
}

# ── umx analysis ──────────────────────────────────────────────────────────────
analyse_umx <- function(condition, dat, fixed_objects) {
    timed_analysis(tryCatch(
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
                    addCI = FALSE,
                    autoRun = FALSE
                )
            ))
            # Keep UMX's default ordinal identification; normalize its free
            # liability variance instead of assuming it equals one.
            m <- mxModel(
                m,
                mxAlgebra(top.A / top.Vtot, name = "h2_ace"),
                mxAlgebra(top.C / top.Vtot, name = "c2_ace")
            )
            m <- suppressMessages(suppressWarnings(mxRun(m, silent = TRUE)))

            a_val <- m$top$a$values[1, 1]
            c_val <- m$top$c$values[1, 1]
            e_val <- m$top$e$values[1, 1]

            total_var <- a_val^2 + c_val^2 + e_val^2
            h2 <- a_val^2 / total_var
            c2 <- c_val^2 / total_var
            e2 <- e_val^2 / total_var

            # Optimization success is distinct from SE/CI availability.
            estimates_ok <- all(is.finite(c(h2, c2, e2))) && total_var > 0
            status <- m$output$status$code
            converged <- as.integer(
                estimates_ok && length(status) == 1L && is.finite(status) &&
                    status %in% c(0L, 1L) && is.finite(m$output$fit)
            )

            # Zero components are admissible boundary estimates.
            admissible <- as.integer(
                estimates_ok && h2 >= -1e-4 && c2 >= -1e-4 && h2 + c2 <= 1 + 1e-4
            )
            boundary_h2 <- as.integer(estimates_ok && h2 <= 1e-4)
            boundary_c2 <- as.integer(estimates_ok && c2 <= 1e-4)
            boundary_e2 <- as.integer(estimates_ok && e2 <= 1e-4)

            # mxSE propagates the full parameter covariance through each ratio.
            get_se <- function(nm) {
                se <- tryCatch(as.numeric(mxSE(nm, m)), error = function(e) NA_real_)
                if (length(se) != 1L || !is.finite(se) || se <= 0) NA_real_ else se
            }
            h2_se <- get_se("h2_ace")
            c2_se <- get_se("c2_ace")
            se_ok <- as.integer(all(is.finite(c(h2_se, c2_se))))

            # 90% profile CI on the standardized proportions themselves.
            ci_bounds <- tryCatch(
                {
                    if (!converged) stop("Optimization did not converge")
                    m_ci <- mxModel(
                        m,
                        mxCI(c("h2_ace[1,1]", "c2_ace[1,1]"), interval = 0.90)
                    )
                    m_ci <- suppressMessages(suppressWarnings(
                        mxRun(m_ci, intervals = TRUE, silent = TRUE)
                    ))
                    list(
                        h2 = get_profile_ci(m_ci$output, "h2_ace[1,1]"),
                        c2 = get_profile_ci(m_ci$output, "c2_ace[1,1]")
                    )
                },
                error = function(e) {
                    list(
                        h2 = c(NA_real_, NA_real_),
                        c2 = c(NA_real_, NA_real_)
                    )
                }
            )

            h2_ci_ok <- as.integer(all(is.finite(ci_bounds$h2)))
            c2_ci_ok <- as.integer(all(is.finite(ci_bounds$c2)))
            profile_ci_ok <- as.integer(h2_ci_ok && c2_ci_ok)

            c(
                converged = converged,
                admissible = admissible,
                optimizer_status = if (length(status) == 1L) status else NA_real_,
                se_ok = se_ok,
                boundary_h2 = boundary_h2,
                boundary_c2 = boundary_c2,
                boundary_e2 = boundary_e2,
                h2_ci_ok = h2_ci_ok,
                c2_ci_ok = c2_ci_ok,
                profile_ci_ok = profile_ci_ok,
                h2_est = h2,
                h2_lo = ci_bounds$h2[1], # profile CI lower
                h2_hi = ci_bounds$h2[2], # profile CI upper
                h2_lo_wald = h2 - qnorm(0.95) * h2_se,
                h2_hi_wald = h2 + qnorm(0.95) * h2_se,
                h2_sd = h2_se, # asymptotic SE
                c2_est = c2,
                c2_lo = ci_bounds$c2[1],
                c2_hi = ci_bounds$c2[2],
                c2_lo_wald = c2 - qnorm(0.95) * c2_se,
                c2_hi_wald = c2 + qnorm(0.95) * c2_se,
                c2_sd = c2_se
            )
        },
        error = function(e) {
            c(
                converged = 0L,
                admissible = NA_integer_,
                optimizer_status = NA_integer_,
                se_ok = 0L,
                boundary_h2 = NA_integer_,
                boundary_c2 = NA_integer_,
                boundary_e2 = NA_integer_,
                h2_ci_ok = 0L,
                c2_ci_ok = 0L,
                profile_ci_ok = 0L,
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
    ))
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
    results <- as.matrix(results)
    h2_true <- fixed_objects$h2_true
    c2_true <- fixed_objects$c2_true

    converged <- get_cols(results, "\\.converged$")
    profile_ci_ok <- get_cols(results, "\\.profile_ci_ok$") # umx only
    se_ok <- get_cols(results, "\\.se_ok$") # umx only
    elapsed <- get_cols(results, "\\.elapsed_sec$") # all attempts, including failures

    # Performance/admissibility/boundary metrics use accepted fits consistently.
    # Failure and CI-success rates use all attempted replications.
    accepted_cols <- function(x) {
        for (j in seq_len(ncol(x))) {
            method <- sub("\\.[^.]+$", "", colnames(x)[j])
            ok <- converged[, paste0(method, ".converged")]
            x[!is.finite(ok) | ok != 1, j] <- NA_real_
        }
        x
    }
    admissible <- accepted_cols(get_cols(results, "\\.admissible$"))
    boundary <- accepted_cols(get_cols(results, "\\.boundary_(h2|c2|e2)$"))
    for (nm in colnames(results)[grepl("\\.(h2|c2)_", colnames(results))]) {
        method <- sub("\\.[^.]+$", "", nm)
        ok <- converged[, paste0(method, ".converged")]
        results[!is.finite(ok) | ok != 1, nm] <- NA_real_
    }

    h2_est <- get_cols(results, "\\.h2_est$")
    h2_est_med <- get_cols(results, "\\.h2_est_med$") # stan only
    h2_est_mode <- get_cols(results, "\\.h2_est_mode$") # stan only
    h2_lo <- get_cols(results, "\\.h2_lo$") # stan: posterior; umx: profile
    h2_hi <- get_cols(results, "\\.h2_hi$")
    h2_lo_hdi <- get_cols(results, "\\.h2_lo_hdi$") # stan HDI
    h2_hi_hdi <- get_cols(results, "\\.h2_hi_hdi$")
    h2_lo_wald <- get_cols(results, "\\.h2_lo_wald$") # umx Wald only
    h2_hi_wald <- get_cols(results, "\\.h2_hi_wald$")
    h2_sd <- get_cols(results, "\\.h2_sd$")

    c2_est <- get_cols(results, "\\.c2_est$")
    c2_est_med <- get_cols(results, "\\.c2_est_med$") # stan only
    c2_est_mode <- get_cols(results, "\\.c2_est_mode$") # stan only
    c2_lo <- get_cols(results, "\\.c2_lo$")
    c2_hi <- get_cols(results, "\\.c2_hi$")
    c2_lo_hdi <- get_cols(results, "\\.c2_lo_hdi$") # stan HDI
    c2_hi_hdi <- get_cols(results, "\\.c2_hi_hdi$")
    c2_lo_wald <- get_cols(results, "\\.c2_lo_wald$")
    c2_hi_wald <- get_cols(results, "\\.c2_hi_wald$")
    c2_sd <- get_cols(results, "\\.c2_sd$")

    # Coverage/width are conditional on successful intervals; report their counts.
    metrics <- function(est, lo, hi, sd_col, true) {
        err <- est - true
        m <- rbind(
            bias = colMeans(err, na.rm = TRUE),
            rmse = sqrt(colMeans(err^2, na.rm = TRUE)),
            emp_sd = apply(est, 2, sd, na.rm = TRUE),
            avg_se = colMeans(sd_col, na.rm = TRUE),
            coverage = colMeans(lo <= true & true <= hi, na.rm = TRUE),
            ci_width = colMeans(hi - lo, na.rm = TRUE),
            n_est = colSums(is.finite(est)),
            n_ci = colSums(is.finite(lo) & is.finite(hi)),
            n_sd = colSums(is.finite(sd_col))
        )
        m[is.nan(m)] <- NA_real_
        colnames(m) <- sub("\\.[^.]+$", "", colnames(m))
        m
    }

    h2_m <- metrics(h2_est, h2_lo, h2_hi, h2_sd, h2_true)
    c2_m <- metrics(c2_est, c2_lo, c2_hi, c2_sd, c2_true)

    # Point-estimate-only metrics for median and mode (stan only)
    pt_metrics <- function(est_df, true) {
        if (ncol(est_df) == 0) {
            return(NULL)
        }
        err <- est_df - true
        nms <- sub("\\.[^.]+$", "", colnames(est_df))
        m <- rbind(
            bias = setNames(colMeans(err, na.rm = TRUE), nms),
            rmse = setNames(sqrt(colMeans(err^2, na.rm = TRUE)), nms),
            emp_sd = setNames(apply(est_df, 2, sd, na.rm = TRUE), nms),
            n_est = setNames(colSums(is.finite(est_df)), nms)
        )
        m[is.nan(m)] <- NA_real_
        m
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
    h2_cov_hdi <- clean(colMeans(
        h2_lo_hdi <= h2_true & h2_true <= h2_hi_hdi,
        na.rm = TRUE
    ))
    h2_wid_hdi <- clean(colMeans(h2_hi_hdi - h2_lo_hdi, na.rm = TRUE))
    h2_n_ci_hdi <- clean(colSums(is.finite(h2_lo_hdi) & is.finite(h2_hi_hdi)))
    c2_cov_wald <- clean(colMeans(
        c2_lo_wald <= c2_true & c2_true <= c2_hi_wald,
        na.rm = TRUE
    ))
    c2_wid_wald <- clean(colMeans(c2_hi_wald - c2_lo_wald, na.rm = TRUE))
    c2_cov_hdi <- clean(colMeans(
        c2_lo_hdi <= c2_true & c2_true <= c2_hi_hdi,
        na.rm = TRUE
    ))
    c2_wid_hdi <- clean(colMeans(c2_hi_hdi - c2_lo_hdi, na.rm = TRUE))
    c2_n_ci_hdi <- clean(colSums(is.finite(c2_lo_hdi) & is.finite(c2_hi_hdi)))

    conv_rate <- clean(colMeans(is.finite(converged) & converged == 1))
    admis_rate <- clean(colMeans(admissible == 1, na.rm = TRUE))
    prof_ci_rate <- clean(colMeans(is.finite(profile_ci_ok) & profile_ci_ok == 1))
    se_rate <- clean(colMeans(is.finite(se_ok) & se_ok == 1))
    boundary_rate <- colMeans(boundary == 1, na.rm = TRUE)
    elapsed[!is.finite(elapsed)] <- NA_real_
    time_mean <- clean(colMeans(elapsed, na.rm = TRUE))
    time_median <- clean(apply(elapsed, 2, median, na.rm = TRUE))
    n_timed <- clean(colSums(is.finite(elapsed)))

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
            ),
            setNames(m["n_est", ], paste0(param, "_n_est_", estimator, ".", names(m["n_est", ])))
        )
    }
    emit_counts <- function(m, param) {
        unlist(lapply(c("n_est", "n_ci", "n_sd"), function(stat) {
            setNames(m[stat, ], paste0(param, "_", stat, ".", colnames(m)))
        }))
    }

    out <- c(
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
        setNames(h2_cov_hdi, paste0("h2_coverage_hdi.", names(h2_cov_hdi))),
        setNames(h2_wid_hdi, paste0("h2_ciwidth_hdi.", names(h2_wid_hdi))),
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
        setNames(c2_cov_hdi, paste0("c2_coverage_hdi.", names(c2_cov_hdi))),
        setNames(c2_wid_hdi, paste0("c2_ciwidth_hdi.", names(c2_wid_hdi))),
        # h² median and mode (stan only)
        emit_pt(h2_med_m, "h2", "med"),
        emit_pt(h2_mode_m, "h2", "mode"),
        # c² median and mode (stan only)
        emit_pt(c2_med_m, "c2", "med"),
        emit_pt(c2_mode_m, "c2", "mode"),
        # Rates
        setNames(conv_rate, paste0("converge_rate.", names(conv_rate))),
        setNames(admis_rate, paste0("admissible_rate.", names(admis_rate))),
        setNames(prof_ci_rate, paste0("profile_ci_rate.", names(prof_ci_rate))),
        setNames(se_rate, paste0("se_rate.", names(se_rate))),
        setNames(time_mean, paste0("time_mean_sec.", names(time_mean))),
        setNames(time_median, paste0("time_median_sec.", names(time_median))),
        setNames(n_timed, paste0("n_timed.", names(n_timed))),
        setNames(boundary_rate, sub("^([^.]+)\\.(.*)$", "\\2_rate.\\1", names(boundary_rate))),
        # Explicit denominators; interval counts also give success rates via n_attempted.
        emit_counts(h2_m, "h2"),
        emit_counts(c2_m, "c2"),
        setNames(clean(colSums(is.finite(h2_lo_wald) & is.finite(h2_hi_wald))), "h2_n_ci_wald.umx"),
        setNames(clean(colSums(is.finite(c2_lo_wald) & is.finite(c2_hi_wald))), "c2_n_ci_wald.umx"),
        setNames(h2_n_ci_hdi, paste0("h2_n_ci_hdi.", names(h2_n_ci_hdi))),
        setNames(c2_n_ci_hdi, paste0("c2_n_ci_hdi.", names(c2_n_ci_hdi))),
        setNames(rep(nrow(results), ncol(converged)), paste0("n_attempted.", names(conv_rate))),
        setNames(clean(colSums(is.finite(converged) & converged == 1)), paste0("n_accepted.", names(conv_rate))),
        setNames(clean(colSums(is.finite(admissible))), paste0("n_admissible_assessed.", names(admis_rate)))
    )
    out[is.nan(out)] <- NA_real_
    out
}

# ── Run simulation ─────────────────────────────────────────────────────────────
# Each Stan call uses 1 core (chains run sequentially); parallelism is across
# replicates via runSimulation().

out <- runSimulation(
    design = DESIGN,
    replications = 2500L,
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
    packages = c("cmdstanr", "posterior", "umx", "polycor", "R.utils", "HDInterval"),
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

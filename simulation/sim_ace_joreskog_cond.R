library(here)
library(cmdstanr)

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

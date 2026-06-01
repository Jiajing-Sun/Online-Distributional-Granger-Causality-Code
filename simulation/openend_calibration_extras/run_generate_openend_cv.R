
# ==============================================================
# run_generate_openend_cv.R
# ==============================================================
# Self-contained script to generate open-end (T = Inf) critical values
# for SSMS, RSMS, and HAC KS/CvM monitoring rules.

args <- commandArgs(trailingOnly = FALSE)
file.arg <- grep("--file=", args, value = TRUE)
if (length(file.arg) > 0L) {
  script.path <- normalizePath(sub("--file=", "", file.arg[1]))
  setwd(dirname(script.path))
}

source(file.path("R", "utils_openend.R"))
source(file.path("R", "openend_weights.R"))
source(file.path("R", "openend_limit_statistics.R"))
source(file.path("R", "openend_critical_values.R"))

# --------------------------------------------------------------
# User controls
# --------------------------------------------------------------
q_max <- 10L
gamma_vec <- c(0, 0.15)
weight_names <- c("U", "Early", "Late", "Mid")
alpha_levels <- c(0.05, 0.10)

nrep <- 5000L
n_train_grid <- 1500L
n_open_grid <- 2000L
ncores <- max(1L, parallel::detectCores() - 1L)
seed <- 13579L

# --------------------------------------------------------------
# Run simulation
# --------------------------------------------------------------
sim <- simulate_openend_critical_values(
  q_max = q_max,
  gamma_vec = gamma_vec,
  weight_names = weight_names,
  alpha_levels = alpha_levels,
  nrep = nrep,
  n_train_grid = n_train_grid,
  n_open_grid = n_open_grid,
  ncores = ncores,
  seed = seed,
  verbose = TRUE
)

paths <- write_openend_critical_values(sim, out_dir = "outputs", prefix = "openend_critical_values")

message("Done. Files written to:")
message("  ", normalizePath(paths$ks_path))
message("  ", normalizePath(paths$cvm_path))

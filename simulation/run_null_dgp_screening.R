rm(list = ls())

.bootstrap_get_script_path <- function() {
  frames <- sys.frames()
  ofiles <- vapply(frames, function(fr) {
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      val <- get("ofile", envir = fr, inherits = FALSE)
      if (length(val) && !is.null(val) && nzchar(val[1])) return(as.character(val[1]))
    }
    NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0L) return(normalizePath(ofiles[length(ofiles)], winslash = "/", mustWork = FALSE))

  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  if (length(file.arg) > 0L) return(normalizePath(sub("^--file=", "", file.arg[1]), winslash = "/", mustWork = FALSE))
  ""
}

.bootstrap_find_root <- function(start_dir, max_up = 12L) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) cur <- dirname(cur)
  for (ii in seq_len(max_up + 1L)) {
    if (file.exists(file.path(cur, "R", "utils.R"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

.env_int <- function(name, default = NULL) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  out <- suppressWarnings(as.integer(val))
  if (is.na(out)) default else out
}

.env_flag <- function(name, default = FALSE) {
  val <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(val)) return(default)
  val %in% c("1", "true", "yes", "y")
}

.env_chr <- function(name, default = "") {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) default else val
}

.this_file <- .bootstrap_get_script_path()
.this_dir <- if (nzchar(.this_file)) dirname(.this_file) else getwd()
setwd(.bootstrap_find_root(.this_dir))
rm(.bootstrap_get_script_path, .bootstrap_find_root, .this_file, .this_dir)

source(file.path("R", "utils.R"))
set_script_wd()
source_project("R", "joe_redesign_runner.R")

make_screen_tau <- function(tau_grid_id) {
  switch(
    tau_grid_id,
    Q5_tail = c(0.05, 0.10, 0.50, 0.90, 0.95),
    Q3_midtail = c(0.10, 0.50, 0.90),
    Q3_tail = c(0.05, 0.50, 0.95),
    stop("Unknown screening tau_grid_id: ", tau_grid_id)
  )
}

make_screen_variants <- function(scope = c("triage", "focused")) {
  scope <- match.arg(scope)
  base <- data.frame(
    dgp_variant_id = c("baseline_ar05", "low_ar03_mild_garch", "low_ar02_milder_garch", "iid_mild_garch"),
    phi_y = c(0.50, 0.30, 0.20, 0.00),
    phi_z = c(0.50, 0.30, 0.20, 0.00),
    n3_df = c(5, 7, 8, 8),
    n4_phi_z = c(0.80, 0.65, 0.60, 0.50),
    garch_alpha = c(0.05, 0.04, 0.03, 0.03),
    garch_beta = c(0.85, 0.80, 0.75, 0.70),
    stringsAsFactors = FALSE
  )
  base$garch_omega <- 1 - base$garch_alpha - base$garch_beta
  if (scope == "focused") {
    base <- base[base$dgp_variant_id %in% c("low_ar03_mild_garch", "low_ar02_milder_garch"), , drop = FALSE]
  }
  base
}

build_null_dgp_screen_grid <- function(scope = c("triage", "focused")) {
  scope <- match.arg(scope)
  variants <- make_screen_variants(scope)
  design_vec <- if (scope == "triage") c("N1", "N2") else c("N1", "N2", "N3", "N4")
  tau_vec <- if (scope == "triage") c("Q5_tail", "Q3_midtail") else c("Q3_midtail", "Q3_tail")
  inst_vec <- c("z", "asym", "scale")
  m_vec <- c(200, 500)
  T_vec <- if (scope == "triage") 5 else c(1, 2, 5)

  rows <- list()
  ii <- 1L
  for (vv in seq_len(nrow(variants))) {
    for (design_id in design_vec) {
      for (tau_grid_id0 in tau_vec) {
        tau_grid0 <- make_screen_tau(tau_grid_id0)
        inst_use <- if (design_id == "N4") "scale" else inst_vec
        for (inst in inst_use) for (m0 in m_vec) for (T0 in T_vec) {
          rows[[ii]] <- data.frame(
            run_id = ii,
            design_id = design_id,
            design_family = "null",
            null_match_id = design_id,
            dgp_variant_id = variants$dgp_variant_id[vv],
            phi_y = variants$phi_y[vv],
            phi_z = variants$phi_z[vv],
            n3_df = variants$n3_df[vv],
            n4_phi_z = variants$n4_phi_z[vv],
            garch_omega = variants$garch_omega[vv],
            garch_alpha = variants$garch_alpha[vv],
            garch_beta = variants$garch_beta[vv],
            model_type = "quantile",
            tau_grid_id = tau_grid_id0,
            tau_grid = paste(tau_grid0, collapse = ";"),
            tau_weight_scheme = "equal",
            instrument_id = inst,
            instrument_description = instrument_description(inst),
            q_expected = length(tau_grid0) * ifelse(inst == "z", 1L, 2L),
            m = m0,
            T = T0,
            n_monitor = as.integer(round(m0 * T0)),
            c_val = 0,
            break_frac = NA_real_,
            k_star = NA_integer_,
            include_eprocess = FALSE,
            eprocess_tau_grid_id = NA_character_,
            eprocess_feature_id = ifelse(inst == "asym", "zminus", "z"),
            notes = "small-B null DGP size screening",
            stringsAsFactors = FALSE
          )
          ii <- ii + 1L
        }
      }
    }
  }
  do.call(rbind, rows)
}

summarise_screen_configs <- function(size_summary, alpha = 0.05, size_cap = 0.15) {
  ss <- size_summary[size_summary$threshold_type == "brownian", , drop = FALSE]
  ss$size_distortion <- ss$empirical_size - alpha
  ss$abs_size_distortion <- abs(ss$size_distortion)
  ss$passes_size_cap <- ss$empirical_size <= size_cap
  write_csv(ss, file.path(output_dir, "null_dgp_screen_size_summary.csv"))

  key_cols <- c("dgp_variant_id", "design_id", "tau_grid_id", "instrument_id", "m", "T", "q")
  key <- make_group_key(ss, key_cols)
  groups <- split(seq_len(nrow(ss)), key)
  cfg <- lapply(groups, function(idx) {
    x <- ss[idx, , drop = FALSE]
    primary <- x[x$method_id %in% c("SSMS_KS_g0", "SSMS_CvM_Late"), , drop = FALSE]
    data.frame(
      dgp_variant_id = x$dgp_variant_id[1],
      design_id = x$design_id[1],
      tau_grid_id = x$tau_grid_id[1],
      instrument_id = x$instrument_id[1],
      m = x$m[1],
      T = x$T[1],
      q = x$q[1],
      max_brownian_size_all_methods = max(x$empirical_size, na.rm = TRUE),
      max_brownian_size_ssms = max(primary$empirical_size, na.rm = TRUE),
      mean_brownian_size_all_methods = mean(x$empirical_size, na.rm = TRUE),
      all_methods_pass_015 = all(x$empirical_size <= size_cap, na.rm = TRUE),
      ssms_pass_015 = all(primary$empirical_size <= size_cap, na.rm = TRUE),
      worst_method = x$method_id[which.max(x$empirical_size)],
      stringsAsFactors = FALSE
    )
  })
  cfg <- do.call(rbind, cfg)
  cfg <- cfg[order(cfg$max_brownian_size_ssms, cfg$max_brownian_size_all_methods), , drop = FALSE]
  write_csv(cfg, file.path(output_dir, "null_dgp_screen_config_summary.csv"))

  winners <- cfg[cfg$ssms_pass_015, , drop = FALSE]
  if (nrow(winners) > 0L) {
    winners <- winners[order(winners$max_brownian_size_ssms, winners$max_brownian_size_all_methods), , drop = FALSE]
  }
  write_csv(utils::head(winners, 40), file.path(output_dir, "null_dgp_screen_recommended_configs.csv"))
  invisible(list(size_summary = ss, config_summary = cfg, recommended = winners))
}

B <- .env_int("JOE_SCREEN_B", default = 200L)
ncores <- .env_int("JOE_NCORES", default = NULL)
reserve <- .env_int("JOE_RESERVE_CORES", default = 3L)
resume <- .env_flag("JOE_RESUME", default = TRUE)
checkpoint_batch_size <- .env_int("JOE_CHECKPOINT_BATCH_SIZE", default = NULL)
scope <- match.arg(.env_chr("JOE_SCREEN_SCOPE", default = "triage"), c("triage", "focused"))
alpha <- 0.05
size_cap <- as.numeric(.env_chr("JOE_SIZE_CAP", default = "0.15"))
if (!is.finite(size_cap)) size_cap <- 0.15

output_dir <- .env_chr("JOE_OUTPUT_DIR", default = "")
if (!nzchar(output_dir)) {
  output_dir <- file.path("output", paste0("null_dgp_screen_", scope, "_B", B, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}
if (!grepl("^(/|[A-Za-z]:[/\\])", output_dir)) output_dir <- project_path(output_dir)
dir_create(output_dir)

ncores <- reserve_cores_n(reserve = reserve, override = ncores)
message("Detected logical cores: ", parallel::detectCores(logical = TRUE))
message(sprintf("Null DGP screen scope=%s, B=%s, ncores=%s (reserve=%s), size cap=%0.3f, resume=%s",
                scope, B, ncores, reserve, size_cap, resume))
message("Output directory: ", output_dir)

cv <- load_critical_values(
  path_base = project_path("critical_values", "critical_values_all.csv"),
  path_weights = project_path("critical_values", "critical_values_all_weights.csv")
)
method_grid <- make_main_method_grid(include_eprocess = FALSE)
design_grid <- build_null_dgp_screen_grid(scope = scope)

write_csv(design_grid, file.path(output_dir, "null_dgp_screen_design_grid.csv"))
write_csv(method_grid, file.path(output_dir, "method_grid.csv"))

long_path <- file.path(output_dir, "null_dgp_screen_replication_long.csv")
screen_long <- run_replication_grid(
  grid = design_grid,
  B = B,
  method_grid = method_grid,
  cv = cv,
  cv_map = NULL,
  alpha = alpha,
  seed_base = 20260529,
  phase = "evaluation",
  ncores = ncores,
  prefer_quantreg = TRUE,
  output_path = long_path,
  resume = resume,
  checkpoint_batch_size = checkpoint_batch_size
)

screen_size <- summarise_size(screen_long)
screen_out <- summarise_screen_configs(screen_size, alpha = alpha, size_cap = size_cap)

manifest <- data.frame(
  scope = scope,
  B = B,
  ncores = ncores,
  reserved_cores = reserve,
  resume = resume,
  checkpoint_batch_size = ifelse(is.null(checkpoint_batch_size), NA_integer_, checkpoint_batch_size),
  alpha = alpha,
  size_cap = size_cap,
  started_or_completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  output_dir = output_dir,
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(output_dir, "run_manifest.csv"))

message("Best SSMS-size configurations:")
print(utils::head(screen_out$config_summary[, c(
  "dgp_variant_id", "design_id", "tau_grid_id", "instrument_id", "m", "T", "q",
  "max_brownian_size_ssms", "max_brownian_size_all_methods", "ssms_pass_015",
  "all_methods_pass_015", "worst_method"
)], 20), row.names = FALSE)


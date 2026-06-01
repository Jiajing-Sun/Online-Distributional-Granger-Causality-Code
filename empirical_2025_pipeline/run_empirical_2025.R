#!/usr/bin/env Rscript

# April 2025--March 2026 BTC--Deribit empirical pipeline.
# This script writes results under empirical_2025_pipeline/output_apr2025_mar2026_ncv20000 by default.

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
                           winslash = "/", mustWork = FALSE)
if (!nzchar(this_file) || is.na(this_file)) this_file <- normalizePath("empirical_2025_pipeline/run_empirical_2025.R", winslash = "/", mustWork = FALSE)
root <- dirname(this_file)

source(file.path(root, "R", "00_config.R"))
source(file.path(root, "R", "01_build_panel.R"))
source(file.path(root, "R", "02_monitoring_core.R"))
source(file.path(root, "R", "03_strategy_robustness.R"))

args <- emp2025_parse_args()
cfg <- emp2025_default_config(args)
emp2025_prepare_dirs(cfg)
set.seed(cfg$seed)
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(cfg$ncores)
}

message("Output directory: ", cfg$output_dir)
message("Mode: ", if (isTRUE(cfg$smoke)) "smoke" else "full")
message("Thread cap: ", cfg$ncores)
if (isTRUE(cfg$smoke) && isTRUE(cfg$smoke_synthetic)) {
  message("Smoke data source: synthetic panel (use --smoke_synthetic=false to exercise raw ingest).")
}
if (nzchar(cfg$panel_csv)) message("Panel override: ", cfg$panel_csv)

panel <- emp2025_build_panel(cfg)
panel_path <- file.path(cfg$output_dir, "panel", "panel_hourly_2025.csv")
message("Panel rows: ", nrow(panel), " (", min(panel$hour_end), " to ", max(panel$hour_end), ")")

main_alarm_path <- file.path(cfg$output_dir, "main", "alarm_table.csv")
if (file.exists(main_alarm_path) && !isTRUE(cfg$force)) {
  alarm_table <- data.table::fread(main_alarm_path)
  alarm_table[, `:=`(
    train_start_time = as.POSIXct(train_start_time, tz = "UTC"),
    train_end_time = as.POSIXct(train_end_time, tz = "UTC"),
    monitor_start_time = as.POSIXct(monitor_start_time, tz = "UTC"),
    monitor_end_time = as.POSIXct(monitor_end_time, tz = "UTC"),
    alarm_time = as.POSIXct(alarm_time, tz = "UTC"),
    trade_start_time = as.POSIXct(trade_start_time, tz = "UTC")
  )]
} else {
  alarm_table <- emp2025_run_windows(panel, cfg, label = "option_block")
  emp2025_write_csv(alarm_table, main_alarm_path)
}

strategy <- emp2025_make_strategy_outputs(panel, alarm_table, cfg)
quality <- emp2025_alarm_quality(panel, alarm_table, cfg)
sparsity <- emp2025_alarm_sparsity(panel, alarm_table, strategy$strategy_summary, quality)

emp2025_write_csv(strategy$strategy_summary, file.path(cfg$output_dir, "main", "strategy_summary.csv"))
emp2025_write_csv(strategy$strategy_returns, file.path(cfg$output_dir, "main", "strategy_returns_long.csv"))
emp2025_write_csv(quality, file.path(cfg$output_dir, "main", "alarm_quality_summary.csv"))
emp2025_write_csv(sparsity, file.path(cfg$output_dir, "main", "alarm_sparsity_summary.csv"))

audit <- emp2025_no_lookahead_audit(alarm_table)
emp2025_write_csv(audit, file.path(cfg$output_dir, "diagnostics", "no_lookahead_audit.csv"))

baseline_path <- file.path(cfg$output_dir, "diagnostics", "baseline_calibration_monitor.csv")
if (file.exists(baseline_path) && !isTRUE(cfg$force)) {
  baseline <- data.table::fread(baseline_path)
} else {
  baseline <- emp2025_run_baseline_monitor(panel, cfg)
  emp2025_write_csv(baseline, baseline_path)
}

baseline_overlap <- merge(
  alarm_table[, .(window_id, method, option_reject = !is.na(alarm_abs_idx), option_alarm_time = alarm_time)],
  baseline[, .(window_id, method, baseline_reject = reject, baseline_alarm_time = first_alarm_time)],
  by = c("window_id", "method"), all.x = TRUE
)
emp2025_write_csv(baseline_overlap, file.path(cfg$output_dir, "diagnostics", "baseline_vs_option_alarm_overlap.csv"))

placebo_types <- c("lagged_spot_only", "block_permuted_week", "delayed_option_168h")
placebo_alarm_path <- file.path(cfg$output_dir, "placebo", "placebo_alarm_table.csv")
if (file.exists(placebo_alarm_path) && !isTRUE(cfg$force)) {
  placebo_alarm <- data.table::fread(placebo_alarm_path)
} else {
  plist <- list()
  for (ptype in placebo_types) {
    message("Running placebo: ", ptype)
    Zp <- emp2025_make_placebo_z(panel, ptype, cfg)
    plist[[ptype]] <- emp2025_run_windows(panel, cfg, z_override = Zp, label = ptype)
  }
  placebo_alarm <- data.table::rbindlist(plist, fill = TRUE)
  emp2025_write_csv(placebo_alarm, placebo_alarm_path)
}
emp2025_write_csv(emp2025_summarize_placebo(placebo_alarm),
                  file.path(cfg$output_dir, "placebo", "placebo_monitoring_summary.csv"))

tc_path <- file.path(cfg$output_dir, "robustness", "transaction_cost_sensitivity.csv")
if (file.exists(tc_path) && !isTRUE(cfg$force)) {
  tc <- data.table::fread(tc_path)
} else {
  tc <- emp2025_transaction_cost_sensitivity(panel, alarm_table, cfg)
  emp2025_write_csv(tc, tc_path)
}

rap_path <- file.path(cfg$output_dir, "robustness", "random_alarm_placebo.csv")
if (file.exists(rap_path) && !isTRUE(cfg$force)) {
  rap <- data.table::fread(rap_path)
} else {
  rap <- emp2025_random_alarm_placebo(panel, alarm_table, cfg)
  emp2025_write_csv(rap, rap_path)
}

boot_path <- file.path(cfg$output_dir, "robustness", "bootstrap_strategy_uncertainty.csv")
if (file.exists(boot_path) && !isTRUE(cfg$force)) {
  boot <- data.table::fread(boot_path)
} else {
  boot <- emp2025_bootstrap_uncertainty(strategy$strategy_returns, cfg)
  emp2025_write_csv(boot, boot_path)
}

reality_path <- file.path(cfg$output_dir, "robustness", "reality_check_random_alarm.csv")
if (file.exists(reality_path) && !isTRUE(cfg$force)) {
  reality_check <- data.table::fread(reality_path)
} else {
  reality_check <- emp2025_reality_check_random_alarm(panel, alarm_table, cfg, strategy$strategy_summary)
  emp2025_write_csv(reality_check, reality_path)
}

mt_path <- file.path(cfg$output_dir, "robustness", "multiple_testing_adjustment.csv")
if (file.exists(mt_path) && !isTRUE(cfg$force)) {
  mt <- data.table::fread(mt_path)
} else {
  mt <- emp2025_multiple_testing(rap, strategy$strategy_summary, alarm_table,
                                 reality_check = reality_check)
  emp2025_write_csv(mt, mt_path)
}

metadata <- data.table::data.table(
  run_time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  sample_start = as.character(cfg$sample_start),
  sample_end = as.character(cfg$sample_end),
  panel_csv = cfg$panel_csv,
  smoke_synthetic = cfg$smoke_synthetic,
  panel_rows = nrow(panel),
  training_size = cfg$training_size,
  monitor_size = cfg$monitor_size,
  refit_every = cfg$refit_every,
  n_windows = data.table::uniqueN(alarm_table$window_id),
  n_methods = data.table::uniqueN(alarm_table$method),
  n_cv_sims = ifelse(isTRUE(cfg$smoke), min(25L, cfg$n_cv_sims), cfg$n_cv_sims),
  n_random_alarm = ifelse(isTRUE(cfg$smoke), min(50L, cfg$n_random_alarm), cfg$n_random_alarm),
  n_bootstrap = ifelse(isTRUE(cfg$smoke), min(50L, cfg$n_bootstrap), cfg$n_bootstrap)
)
emp2025_write_csv(metadata, file.path(cfg$output_dir, "run_metadata.csv"))

saveRDS(list(config = cfg, panel = panel, alarm_table = alarm_table,
             strategy = strategy, quality = quality, sparsity = sparsity,
             baseline = baseline, placebo_alarm = placebo_alarm,
             transaction_cost = tc, random_alarm_placebo = rap,
             bootstrap = boot, reality_check = reality_check, multiple_testing = mt),
        file.path(cfg$output_dir, "results_full_2025.rds"))

message("Done. Results written to: ", cfg$output_dir)

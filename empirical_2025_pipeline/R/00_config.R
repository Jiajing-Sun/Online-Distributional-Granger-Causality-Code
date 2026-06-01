# Configuration and small utilities for the April 2025--March 2026 empirical run.

emp2025_root <- function() {
  this <- tryCatch(normalizePath(dirname(sys.frame(1)$ofile), winslash = "/", mustWork = FALSE),
                   error = function(e) "")
  if (nzchar(this) && basename(this) == "R") return(dirname(this))
  wd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if (basename(wd) == "empirical_2025_pipeline") return(wd)
  cand <- file.path(wd, "empirical_2025_pipeline")
  if (dir.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  stop("Run from the JoE project root or empirical_2025_pipeline.")
}

emp2025_project_root <- function() dirname(emp2025_root())

emp2025_path <- function(...) file.path(emp2025_root(), ...)

emp2025_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

emp2025_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  named <- grep("^--[A-Za-z0-9_.-]+=", args, value = TRUE)
  for (x in named) {
    key <- sub("^--([^=]+)=.*$", "\\1", x)
    val <- sub("^--[^=]+=", "", x)
    out[[key]] <- val
  }
  out
}

emp2025_bool <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || !nzchar(x[1])) return(default)
  y <- tolower(trimws(as.character(x[1])))
  if (y %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (y %in% c("0", "false", "f", "no", "n")) return(FALSE)
  default
}

emp2025_num <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x) || !nzchar(x[1])) return(default)
  y <- suppressWarnings(as.numeric(x[1]))
  if (is.finite(y)) y else default
}

emp2025_int <- function(x, default = NA_integer_) {
  y <- emp2025_num(x, default = NA_real_)
  if (is.finite(y)) as.integer(round(y)) else default
}

emp2025_time <- function(x, default) {
  if (is.null(x) || !length(x) || !nzchar(x[1])) return(default)
  y <- suppressWarnings(as.POSIXct(x[1], tz = "UTC"))
  if (is.na(y)) stop("Cannot parse UTC time: ", x[1])
  y
}

emp2025_split_paths <- function(x) {
  if (is.null(x) || !length(x) || !nzchar(x[1])) return(character())
  y <- unlist(strsplit(x[1], ",", fixed = TRUE), use.names = FALSE)
  trimws(y[nzchar(trimws(y))])
}

emp2025_default_config <- function(args = emp2025_parse_args()) {
  smoke <- emp2025_bool(args$smoke, FALSE)
  out_dir <- args$output_dir
  if (is.null(out_dir) || !nzchar(out_dir)) {
    out_dir <- emp2025_path(if (isTRUE(smoke)) "output_apr2025_mar2026_smoke" else "output_apr2025_mar2026_ncv20000")
  }
  default_training <- if (isTRUE(smoke)) 480L else 2160L
  default_monitor <- if (isTRUE(smoke)) 96L else 336L

  list(
    deribit_dir = if (!is.null(args$deribit_dir)) args$deribit_dir else
      Sys.getenv("DERIBIT_DIR", unset = "data/raw_deribit"),
    extra_deribit_files = emp2025_split_paths(args$extra_deribit_files),
    panel_csv = if (!is.null(args$panel_csv)) args$panel_csv else "",
    smoke_synthetic = emp2025_bool(args$smoke_synthetic, smoke),
    output_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    symbol = if (!is.null(args$symbol)) args$symbol else "BTCUSDT",
    binance_interval = "1h",
    sample_start = emp2025_time(args$sample_start, as.POSIXct("2025-04-01 00:00:00", tz = "UTC")),
    sample_end = emp2025_time(args$sample_end, as.POSIXct("2026-03-31 23:00:00", tz = "UTC")),
    keep_outcomes_inside_sample = TRUE,
    training_size = emp2025_int(args$training_size, default_training),
    monitor_size = emp2025_int(args$monitor_size, default_monitor),
    refit_every = emp2025_int(args$refit_every, 24L),
    taus = c(0.05, 0.10),
    alpha = emp2025_num(args$alpha, 0.05),
    gammas = c(0, 0.15),
    cvm_weights = c("U", "Early", "Mid", "Late"),
    theta_dict = c(-1, -0.5, 0.5, 1),
    eprocess_restart_every = 24L,
    hold_periods = c(6L, 24L),
    strategy_sides = c("flat", "short"),
    transaction_cost_bps = c(0, 5, 10, 20),
    short_funding_bps_per_hour = c(0, 0.25, 0.50, 1.00),
    n_cv_sims = emp2025_int(args$n_cv_sims, 20000L),
    n_random_alarm = emp2025_int(args$n_random_alarm, 2000L),
    n_bootstrap = emp2025_int(args$n_bootstrap, 1000L),
    bootstrap_blocks = c(24L, 72L, 168L),
    seed = emp2025_int(args$seed, 20250529L),
    ncores = emp2025_int(args$ncores, max(1L, parallel::detectCores(logical = TRUE) - 2L)),
    smoke = smoke,
    force = emp2025_bool(args$force, FALSE)
  )
}

emp2025_prepare_dirs <- function(cfg) {
  dirs <- file.path(cfg$output_dir, c("panel", "main", "diagnostics", "placebo", "robustness", "cache", "logs"))
  invisible(lapply(dirs, emp2025_dir_create))
  emp2025_dir_create(file.path(cfg$output_dir, "figures"))
  cfg$output_dir
}

emp2025_write_csv <- function(x, path) {
  emp2025_dir_create(dirname(path))
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  data.table::fwrite(x, tmp)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

emp2025_source_monitor_code <- function() {
  root <- emp2025_project_root()
  mon_root <- file.path(root, "code_from_march_email", "empirical")
  if (!dir.exists(mon_root)) stop("Cannot find monitoring code at: ", mon_root)
  options(distgc.project_root = normalizePath(mon_root, winslash = "/", mustWork = TRUE))
  source(file.path(mon_root, "R", "utils.R"), local = .GlobalEnv)
  source(file.path(mon_root, "R", "fit_distributional_models.R"), local = .GlobalEnv)
  source(file.path(mon_root, "R", "monitors.R"), local = .GlobalEnv)
  source(file.path(mon_root, "R", "eprocess.R"), local = .GlobalEnv)
  invisible(mon_root)
}

emp2025_require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(pkgs)
}

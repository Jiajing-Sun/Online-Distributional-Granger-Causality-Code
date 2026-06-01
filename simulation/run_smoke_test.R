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

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (!nzchar(p)) p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
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
.this_file <- .bootstrap_get_script_path()
.this_dir <- if (nzchar(.this_file)) dirname(.this_file) else getwd()
setwd(.bootstrap_find_root(.this_dir))
rm(.bootstrap_get_script_path, .bootstrap_find_root, .this_file, .this_dir)
source(file.path("R", "utils.R"))
set_script_wd()
source_project("R", "simulation_runner.R")

dir_create("output")

nrep <- 1000  # quick smoke test only
ncores <- 1L

cfg <- data.frame(
  design = c("N1", "A1", "A4"),
  m = c(200, 200, 200),
  T = c(1, 1, 1),
  c_val = c(0, 0.25, 0.25),
  break_frac = c(NA, 0.50, 0.50),
  stringsAsFactors = FALSE
)

# Quantile smoke test
cfg_q <- cfg
cfg_q$instrument_type <- c("z", "z", "asym")
run_simulation_grid(
  config_grid = cfg_q,
  model_type = "quantile",
  tau_grid = c(0.05, 0.10, 0.50, 0.90, 0.95),
  include_eprocess = TRUE,
  eprocess_tau_grid = c(0.05, 0.10),
  eprocess_feature = "z",
  nrep = nrep,
  ncores = ncores,
  gamma_vec = c(0, 0.15),
  cvm_weights = c("U", "Late", "Early", "Mid"),
  output_dir = file.path("output", "smoke_quantile"),
  file_stub = "smoke_quantile"
)

# Expectile smoke test
cfg_e <- cfg
cfg_e$instrument_type <- c("z", "z", "asym")
run_simulation_grid(
  config_grid = cfg_e,
  model_type = "expectile",
  tau_grid = c(0.10, 0.25, 0.50, 0.75, 0.90),
  include_eprocess = FALSE,
  nrep = nrep,
  ncores = ncores,
  gamma_vec = c(0, 0.15),
  cvm_weights = c("U", "Late", "Early", "Mid"),
  output_dir = file.path("output", "smoke_expectile"),
  file_stub = "smoke_expectile"
)

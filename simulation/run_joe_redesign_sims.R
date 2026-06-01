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

.env_flag_optional <- function(name) {
  val <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(val)) return(NULL)
  val %in% c("1", "true", "yes", "y")
}

.env_num_vec <- function(name, default = NULL) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) return(default)
  parts <- unlist(strsplit(val, "[,;[:space:]]+"))
  parts <- parts[nzchar(parts)]
  out <- suppressWarnings(as.numeric(parts))
  out <- out[!is.na(out)]
  if (!length(out)) default else out
}

.this_file <- .bootstrap_get_script_path()
.this_dir <- if (nzchar(.this_file)) dirname(.this_file) else getwd()
setwd(.bootstrap_find_root(.this_dir))
rm(.bootstrap_get_script_path, .bootstrap_find_root, .this_file, .this_dir)

source(file.path("R", "utils.R"))
set_script_wd()
source_project("R", "joe_redesign_runner.R")

package <- Sys.getenv("JOE_SIM_PACKAGE", unset = "pilot")
B_cal <- .env_int("JOE_B_CAL", default = NULL)
B_eval <- .env_int("JOE_B_EVAL", default = NULL)
B_alt <- .env_int("JOE_B_ALT", default = NULL)
ncores <- .env_int("JOE_NCORES", default = NULL)
reserve <- .env_int("JOE_RESERVE_CORES", default = 3L)
include_expectile <- .env_flag("JOE_INCLUDE_EXPECTILE", default = FALSE)
include_n3 <- .env_flag_optional("JOE_INCLUDE_N3")
include_n4 <- .env_flag_optional("JOE_INCLUDE_N4")
resume <- .env_flag("JOE_RESUME", default = TRUE)
checkpoint_batch_size <- .env_int("JOE_CHECKPOINT_BATCH_SIZE", default = NULL)
dgp_variant_id <- Sys.getenv("JOE_DGP_VARIANT", unset = "baseline")
tau_grid_override <- Sys.getenv("JOE_TAU_GRID", unset = "")
if (!nzchar(tau_grid_override) || tau_grid_override == "default") tau_grid_override <- NULL
m_override <- .env_num_vec("JOE_M_GRID", default = NULL)
T_override <- .env_num_vec("JOE_T_GRID", default = NULL)
c_override <- .env_num_vec("JOE_C_GRID", default = NULL)
break_override <- .env_num_vec("JOE_BREAK_GRID", default = NULL)
output_dir <- Sys.getenv("JOE_OUTPUT_DIR", unset = "")
if (!nzchar(output_dir)) output_dir <- NULL

message("Detected logical cores: ", parallel::detectCores(logical = TRUE))
message("Requested reserve cores: ", reserve)
message("Resume enabled: ", resume)
if (!is.null(checkpoint_batch_size)) {
  message("Checkpoint batch size: ", checkpoint_batch_size)
}

run_joe_redesign_package(
  package = package,
  B_cal = B_cal,
  B_eval = B_eval,
  B_alt = B_alt,
  ncores = ncores,
  reserve_cores = reserve,
  output_dir = output_dir,
  include_expectile = include_expectile,
  include_n3 = include_n3,
  include_n4 = include_n4,
  tau_grid_override = tau_grid_override,
  dgp_variant_id = dgp_variant_id,
  m_override = m_override,
  T_override = T_override,
  c_override = c_override,
  break_override = break_override,
  resume = resume,
  checkpoint_batch_size = checkpoint_batch_size
)

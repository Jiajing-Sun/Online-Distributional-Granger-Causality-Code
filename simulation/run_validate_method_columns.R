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
source_project("R", "method_catalog.R")

required_quantile <- get_method_catalog("all", "quantile")
required_expectile <- get_method_catalog("all", "expectile")

check_summary <- function(path, required_methods) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  needed <- as.vector(outer(required_methods, c("_rej_rate", "_fa_rate", "_det_rate", "_avg_delay", "_med_delay", "_arl0_mean"), paste0))
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop(sprintf("Missing %s required summary columns in %s. First few: %s",
                 length(missing), path, paste(head(missing, 10), collapse = ", ")), call. = FALSE)
  }
  message("OK: ", path)
}

check_summary(project_path("output", "table1_null_quantile", "table1_null_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "table1_null_expectile", "table1_null_expectile_summary.csv"), required_expectile)
check_summary(project_path("output", "table2_abrupt_quantile", "table2_abrupt_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "table2_abrupt_expectile", "table2_abrupt_expectile_summary.csv"), required_expectile)
check_summary(project_path("output", "table3_a3_quantile", "table3_a3_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "table3_a3_expectile", "table3_a3_expectile_summary.csv"), required_expectile)
check_summary(project_path("output", "table3_a4_quantile", "table3_a4_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "table3_a4_expectile", "table3_a4_expectile_summary.csv"), required_expectile)
check_summary(project_path("output", "table4_gradual_quantile", "table4_gradual_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "table4_gradual_expectile", "table4_gradual_expectile_summary.csv"), required_expectile)
check_summary(project_path("output", "appendix_c1_quantile", "appendix_c1_quantile_summary.csv"), required_quantile)
check_summary(project_path("output", "appendix_c1_expectile", "appendix_c1_expectile_summary.csv"), required_expectile)
message("All required method columns are present.")

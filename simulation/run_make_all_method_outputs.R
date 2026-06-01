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
source_project("make_simulation_summary_figures.R")
source_project("make_full_method_summary_tables.R")

make_all_simulation_summary_figures(results_root = get_project_root(), method_set = "paper")
make_all_simulation_summary_figures(results_root = get_project_root(), method_set = "all")
make_full_method_summary_tables(results_root = get_project_root())

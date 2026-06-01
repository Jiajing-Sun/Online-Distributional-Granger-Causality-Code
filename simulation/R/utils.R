# ==============================================================
# utils.R -- small helpers for the online distributional-GC simulation package
# ==============================================================

is_windows <- function() {
  tolower(Sys.info()[["sysname"]]) == "windows"
}

get_rstudio_script_path <- function() {
  if (!requireNamespace("rstudioapi", quietly = TRUE)) return("")
  if (!rstudioapi::isAvailable()) return("")

  p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (!nzchar(p)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
  }
  if (!nzchar(p)) return("")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

detect_current_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file.arg <- grep("^--file=", args, value = TRUE)
  if (length(file.arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file.arg[1]), winslash = "/", mustWork = FALSE))
  }

  frames <- sys.frames()
  ofiles <- vapply(frames, function(fr) {
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      val <- get("ofile", envir = fr, inherits = FALSE)
      if (length(val) && !is.null(val) && nzchar(val[1])) return(as.character(val[1]))
    }
    NA_character_
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0L) {
    return(normalizePath(ofiles[length(ofiles)], winslash = "/", mustWork = FALSE))
  }

  p <- get_rstudio_script_path()
  if (nzchar(p)) return(p)
  ""
}

find_project_root <- function(start_dir = getwd(), max_up = 12L) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  if (file.exists(cur) && !dir.exists(cur)) {
    cur <- dirname(cur)
  }

  for (ii in seq_len(max_up + 1L)) {
    has_utils <- file.exists(file.path(cur, "R", "utils.R"))
    has_runner <- any(file.exists(file.path(cur, c(
      "run_table1_null.R", "run_table2_abrupt.R", "run_smoke_test.R",
      "run_all_paper_sims.R", "README.md"
    ))))
    if (has_utils && has_runner) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  normalizePath(start_dir, winslash = "/", mustWork = FALSE)
}

get_project_root <- function(refresh = FALSE) {
  root <- getOption("distgc.project_root", default = NULL)
  if (!refresh && !is.null(root) && dir.exists(root) && file.exists(file.path(root, "R", "utils.R"))) {
    return(normalizePath(root, winslash = "/", mustWork = FALSE))
  }

  script.path <- detect_current_script_path()
  start_dir <- if (nzchar(script.path)) dirname(script.path) else getwd()
  root <- find_project_root(start_dir)
  options(distgc.project_root = root)
  normalizePath(root, winslash = "/", mustWork = FALSE)
}

set_script_wd <- function() {
  root <- get_project_root(refresh = TRUE)
  setwd(root)
  options(distgc.project_root = root)
  invisible(root)
}

project_path <- function(...) {
  file.path(get_project_root(), ...)
}

source_project <- function(..., local = parent.frame(), chdir = FALSE) {
  source(project_path(...), local = local, chdir = chdir)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_solve <- function(A, ridge = 1e-10) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  q <- nrow(A)
  A2 <- A + ridge * diag(q)

  out <- tryCatch({
    R <- chol(A2)
    chol2inv(R)
  }, error = function(e) {
    solve(A2)
  })
  out
}

safe_qr_solve <- function(A, b, ridge = 1e-8) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  q <- nrow(A)
  A2 <- A + ridge * diag(q)
  tryCatch(qr.solve(A2, b), error = function(e) solve(A2, b))
}

make_cluster <- function(ncores) {
  if (is.null(ncores) || is.na(ncores) || ncores <= 1L) return(NULL)
  parallel::makeCluster(as.integer(ncores), type = "PSOCK")
}

close_parallel_connections <- function() {
  cons <- try(showConnections(all = TRUE), silent = TRUE)
  if (inherits(cons, "try-error") || is.null(cons) || nrow(cons) == 0L) return(invisible(NULL))

  desc <- cons[, "description"]
  ids <- suppressWarnings(as.integer(rownames(cons)))
  keep <- grepl("localhost", desc, fixed = TRUE) | grepl("sock", desc, ignore.case = TRUE)
  ids <- ids[keep & !is.na(ids)]
  for (id in ids) {
    try(close(getConnection(id)), silent = TRUE)
  }
  invisible(NULL)
}

stop_cluster <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
  close_parallel_connections()
  invisible(NULL)
}

tic <- function() {
  assign(".tic_time__", proc.time()[3], envir = .GlobalEnv)
}

toc <- function(msg = "Elapsed") {
  t0 <- get(".tic_time__", envir = .GlobalEnv)
  dt <- proc.time()[3] - t0
  message(sprintf("%s: %.2f sec", msg, dt))
  invisible(dt)
}

log_sum_exp <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(-Inf)
  mx <- max(x)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}

normalize_weights <- function(w) {
  w <- as.numeric(w)
  if (any(w < 0)) stop("Weights must be nonnegative.")
  s <- sum(w)
  if (s <= 0) stop("Weights must sum to a positive number.")
  w / s
}

mat_cumsum <- function(X) {
  X <- as.matrix(X)
  out <- apply(X, 2, cumsum)
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L)
  }
  colnames(out) <- colnames(X)
  out
}

trimmed_mean <- function(x, trim = 0.05) {
  stats::mean(x, trim = trim, na.rm = TRUE)
}

require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed. Please run install.packages('%s').", pkg, pkg))
  }
}

get_default_ncores <- function(reserve = 1L) {
  nc <- parallel::detectCores(logical = TRUE)
  nc <- if (is.na(nc)) 1L else nc
  max(1L, as.integer(nc - reserve))
}

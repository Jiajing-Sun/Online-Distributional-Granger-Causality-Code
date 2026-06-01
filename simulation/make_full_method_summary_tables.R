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
# ============================================================
# Create pooled CSV and LaTeX tables for all monitoring methods.
# ============================================================

source_project("R", "method_catalog.R")

resolve_results_root <- function(results_root = ".") {
  if (missing(results_root) || is.null(results_root) || identical(results_root, ".")) return(get_project_root())
  if (!grepl("^(/|[A-Za-z]:[/\\])", results_root)) return(project_path(results_root))
  results_root
}

read_summary <- function(results_root, design_name) {
  results_root <- resolve_results_root(results_root)
  f <- file.path(results_root, "output", design_name, paste0(design_name, "_summary.csv"))
  if (!file.exists(f)) stop("Missing summary file: ", f, call. = FALSE)
  read.csv(f, stringsAsFactors = FALSE)
}

pool_metric <- function(df, methods, suffix) {
  out <- sapply(methods, function(m) {
    col <- paste0(m, suffix)
    if (!col %in% names(df)) return(NA_real_)
    mean(df[[col]], na.rm = TRUE)
  })
  as.numeric(out)
}

align_metric <- function(target_methods, source_methods, values) {
  out <- rep(NA_real_, length(target_methods))
  names(out) <- target_methods
  if (length(source_methods) == 0L) return(out)
  idx <- match(source_methods, target_methods)
  keep <- which(!is.na(idx))
  if (length(keep) > 0L) out[idx[keep]] <- values[keep]
  out
}

fmt_pct <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", 100 * x))
fmt_num <- function(x) ifelse(is.na(x), "--", sprintf("%.1f", x))

write_table_tex <- function(df, path, caption, label) {
  cols <- names(df)
  aligns <- paste0("l", paste(rep("r", length(cols) - 1L), collapse = ""))
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", aligns, "}"),
    "\\toprule",
    paste(cols, collapse = " & "), "\\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(df))) {
    lines <- c(lines, paste(df[i, ], collapse = " & "), "\\\\")
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, con = path)
}

make_det_delay_table <- function(results_root, qfile, efile, stem, caption, label, meth_q, meth_e, outdir) {
  qd <- read_summary(results_root, qfile)
  ed <- read_summary(results_root, efile)
  e_det <- align_metric(meth_q, meth_e, pool_metric(ed, meth_e, "_det_rate"))
  e_del <- align_metric(meth_q, meth_e, pool_metric(ed, meth_e, "_avg_delay"))
  df <- data.frame(
    Method = pretty_method_name(meth_q),
    Q_Det = fmt_pct(pool_metric(qd, meth_q, "_det_rate")),
    Q_Delay = fmt_num(pool_metric(qd, meth_q, "_avg_delay")),
    E_Det = fmt_pct(e_det),
    E_Delay = fmt_num(e_del),
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(outdir, paste0(stem, ".csv")), row.names = FALSE)
  write_table_tex(df, file.path(outdir, paste0(stem, ".tex")), caption, label)
}

make_full_method_summary_tables <- function(results_root = ".") {
  results_root <- resolve_results_root(results_root)
  outdir <- file.path(results_root, "summary_all_methods")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  meth_q <- get_method_catalog("all", "quantile")
  meth_e <- get_method_catalog("all", "expectile")

  q1 <- read_summary(results_root, "table1_null_quantile")
  e1 <- read_summary(results_root, "table1_null_expectile")
  size_df <- data.frame(
    Method = pretty_method_name(meth_q),
    Quantile = fmt_pct(pool_metric(q1, meth_q, "_rej_rate")),
    Expectile = fmt_pct(align_metric(meth_q, meth_e, pool_metric(e1, meth_e, "_rej_rate"))),
    stringsAsFactors = FALSE
  )
  write.csv(size_df, file.path(outdir, "table_size_all_methods.csv"), row.names = FALSE)
  write_table_tex(size_df, file.path(outdir, "table_size_all_methods.tex"),
                  "Pooled null rejection frequencies for all monitoring methods.",
                  "tab:allmethods_size")

  make_det_delay_table(results_root, "table2_abrupt_quantile", "table2_abrupt_expectile", "table_abrupt_all_methods",
                       "Pooled abrupt-design detection rates and conditional delays for all monitoring methods.",
                       "tab:allmethods_abrupt", meth_q, meth_e, outdir)
  make_det_delay_table(results_root, "table3_a3_quantile", "table3_a3_expectile", "table_a3_all_methods",
                       "Pooled A3 scale-design detection rates and conditional delays for all monitoring methods.",
                       "tab:allmethods_a3", meth_q, meth_e, outdir)
  make_det_delay_table(results_root, "table3_a4_quantile", "table3_a4_expectile", "table_a4_all_methods",
                       "Pooled A4 downside-tail design detection rates and conditional delays for all monitoring methods.",
                       "tab:allmethods_a4", meth_q, meth_e, outdir)
  make_det_delay_table(results_root, "appendix_c1_quantile", "appendix_c1_expectile", "table_c1_all_methods",
                       "Pooled contaminated-training design detection rates and conditional delays for all monitoring methods.",
                       "tab:allmethods_c1", meth_q, meth_e, outdir)
  make_det_delay_table(results_root, "table4_gradual_quantile", "table4_gradual_expectile", "table_gradual_all_methods",
                       "Pooled gradual-design detection rates and conditional delays for all monitoring methods.",
                       "tab:allmethods_gradual", meth_q, meth_e, outdir)

  invisible(outdir)
}

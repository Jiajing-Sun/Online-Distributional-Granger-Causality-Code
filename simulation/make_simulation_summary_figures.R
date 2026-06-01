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
# Simulation-summary graphics for the online distributional
# causality paper. Supports compact paper figures and
# comprehensive all-method figures.
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
  if (!file.exists(f)) {
    stop("Missing summary file: ", f, call. = FALSE)
  }
  read.csv(f, stringsAsFactors = FALSE)
}

avg_metric <- function(df, suffix, methods) {
  vals <- sapply(methods, function(m) {
    col <- paste0(m, suffix)
    if (!col %in% names(df)) return(NA_real_)
    mean(df[[col]], na.rm = TRUE)
  })
  names(vals) <- pretty_method_name(methods)
  vals
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

save_plot_both <- function(basefile, width = 9, height = 5, expr_fun) {
  pdf(paste0(basefile, ".pdf"), width = width, height = height)
  expr_fun()
  dev.off()

  png(paste0(basefile, ".png"), width = width, height = height, units = "in", res = 300)
  expr_fun()
  dev.off()
}

bar_pattern <- function(n) {
  list(
    density = rep(c(10, 18, 26, 34, 42, 50, 58), length.out = n),
    angle = rep(c(45, -45, 0, 90, 135), length.out = n)
  )
}

method_colors <- function(labels) {
  cols <- vapply(labels, function(x) {
    if (grepl("E-process", x)) return("#6B7280")
    if (grepl("^SSMS", x)) return("#2B6CB0")
    if (grepl("^RSMS", x)) return("#D97706")
    if (grepl("^HAC", x)) return("#2F855A")
    "#4B5563"
  }, character(1))
  unname(cols)
}

plot_bar_two_panel <- function(vals_left, vals_right, main_left, main_right,
                               ylab = "Rate", target = NULL, ylim = NULL) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 2), mar = c(10, 4, 3, 1) + 0.1)
  if (is.null(ylim)) {
    ymax <- max(c(vals_left, vals_right), na.rm = TRUE)
    ylim <- c(0, 1.10 * ymax)
  }
  cols_left <- method_colors(names(vals_left))
  cols_right <- method_colors(names(vals_right))
  pat_left <- bar_pattern(length(vals_left))
  pat_right <- bar_pattern(length(vals_right))

  bp <- barplot(vals_left, las = 2, col = cols_left, border = "black",
                density = pat_left$density, angle = pat_left$angle,
                ylim = ylim, ylab = ylab, main = main_left)
  if (!is.null(target)) abline(h = target, lty = 2)
  lbl_left <- if (grepl("delay", tolower(ylab))) sprintf("%.0f", vals_left) else sprintf("%.1f", 100 * vals_left)
  text(bp, vals_left, labels = lbl_left, pos = 3, cex = 0.65)

  bp <- barplot(vals_right, las = 2, col = cols_right, border = "black",
                density = pat_right$density, angle = pat_right$angle,
                ylim = ylim, ylab = ylab, main = main_right)
  if (!is.null(target)) abline(h = target, lty = 2)
  lbl_right <- if (grepl("delay", tolower(ylab))) sprintf("%.0f", vals_right) else sprintf("%.1f", 100 * vals_right)
  text(bp, vals_right, labels = lbl_right, pos = 3, cex = 0.65)
}

plot_weight_tradeoff <- function(det_q, det_e, del_q, del_e, outfile_base) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  save_plot_both(outfile_base, width = 10, height = 7, expr_fun = function() {
    par(mfrow = c(2, 2), mar = c(5, 5, 3, 1) + 0.1)
    line_cols <- c("#2B6CB0", "#D97706", "#2F855A")
    line_lty <- c(1, 2, 3)
    line_pch <- c(16, 17, 15)
    matplot(det_q, type = "b", pch = line_pch, lty = line_lty, col = line_cols,
            xaxt = "n", xlab = "CvM weight", ylab = "Detection rate", main = "A5 quantile: detection")
    axis(1, at = seq_len(nrow(det_q)), labels = rownames(det_q))
    legend("topleft", legend = colnames(det_q), lty = line_lty, pch = line_pch, col = line_cols, bty = "n")

    matplot(del_q, type = "b", pch = line_pch, lty = line_lty, col = line_cols,
            xaxt = "n", xlab = "CvM weight", ylab = "Average delay", main = "A5 quantile: delay")
    axis(1, at = seq_len(nrow(del_q)), labels = rownames(del_q))
    legend("topleft", legend = colnames(del_q), lty = line_lty, pch = line_pch, col = line_cols, bty = "n")

    matplot(det_e, type = "b", pch = line_pch, lty = line_lty, col = line_cols,
            xaxt = "n", xlab = "CvM weight", ylab = "Detection rate", main = "A5 expectile: detection")
    axis(1, at = seq_len(nrow(det_e)), labels = rownames(det_e))
    legend("topleft", legend = colnames(det_e), lty = line_lty, pch = line_pch, col = line_cols, bty = "n")

    matplot(del_e, type = "b", pch = line_pch, lty = line_lty, col = line_cols,
            xaxt = "n", xlab = "CvM weight", ylab = "Average delay", main = "A5 expectile: delay")
    axis(1, at = seq_len(nrow(del_e)), labels = rownames(del_e))
    legend("topleft", legend = colnames(del_e), lty = line_lty, pch = line_pch, col = line_cols, bty = "n")
  })
}

make_all_simulation_summary_figures <- function(results_root = ".", method_set = c("paper", "all")) {
  method_set <- match.arg(method_set)
  outdir <- file.path(results_root, if (method_set == "paper") "graphics_summary" else "graphics_summary_allmethods")
  ensure_dir(outdir)
  title_suffix <- if (method_set == "all") " (all methods)" else ""

  methods_q <- get_method_catalog(which = method_set, model_type = "quantile")
  methods_e <- get_method_catalog(which = method_set, model_type = "expectile")

  # 1. Null size
  q1 <- read_summary(results_root, "table1_null_quantile")
  e1 <- read_summary(results_root, "table1_null_expectile")
  size_q <- avg_metric(q1, "_rej_rate", methods_q)
  size_e <- avg_metric(e1, "_rej_rate", methods_e)
  save_plot_both(file.path(outdir, "fig_size_overview"), width = ifelse(method_set == "all", 14, 9), expr_fun = function() {
    plot_bar_two_panel(size_q, size_e,
                       main_left = paste0("Quantile size (N1-N2)", title_suffix),
                       main_right = paste0("Expectile size (N1-N2)", title_suffix),
                       ylab = "Average rejection rate",
                       target = 0.05,
                       ylim = c(0, max(c(size_q, size_e), na.rm = TRUE) * 1.15))
  })

  # 2. Abrupt designs
  q2 <- read_summary(results_root, "table2_abrupt_quantile")
  e2 <- read_summary(results_root, "table2_abrupt_expectile")
  det_q2 <- avg_metric(q2, "_det_rate", methods_q)
  det_e2 <- avg_metric(e2, "_det_rate", methods_e)
  del_q2 <- avg_metric(q2, "_avg_delay", methods_q)
  del_e2 <- avg_metric(e2, "_avg_delay", methods_e)
  save_plot_both(file.path(outdir, "fig_abrupt_detection"), width = ifelse(method_set == "all", 14, 9), expr_fun = function() {
    plot_bar_two_panel(det_q2, det_e2,
                       main_left = paste0("Abrupt A1-A2: quantile", title_suffix),
                       main_right = paste0("Abrupt A1-A2: expectile", title_suffix),
                       ylab = "Average detection rate",
                       ylim = c(0, max(c(det_q2, det_e2), na.rm = TRUE) * 1.15))
  })
  save_plot_both(file.path(outdir, "fig_abrupt_delay"), width = ifelse(method_set == "all", 14, 9), expr_fun = function() {
    plot_bar_two_panel(del_q2, del_e2,
                       main_left = paste0("Abrupt A1-A2: quantile", title_suffix),
                       main_right = paste0("Abrupt A1-A2: expectile", title_suffix),
                       ylab = "Average delay",
                       ylim = c(0, max(c(del_q2, del_e2), na.rm = TRUE) * 1.15))
  })

  # 3. Tail designs A3 and A4 detection
  a3q <- avg_metric(read_summary(results_root, "table3_a3_quantile"), "_det_rate", methods_q)
  a3e <- avg_metric(read_summary(results_root, "table3_a3_expectile"), "_det_rate", methods_e)
  a4q <- avg_metric(read_summary(results_root, "table3_a4_quantile"), "_det_rate", methods_q)
  a4e <- avg_metric(read_summary(results_root, "table3_a4_expectile"), "_det_rate", methods_e)
  save_plot_both(file.path(outdir, "fig_tail_designs_detection"), width = ifelse(method_set == "all", 16, 11), height = 8, expr_fun = function() {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar), add = TRUE)
    par(mfrow = c(2, 2), mar = c(10, 4, 3, 1) + 0.1)
    ymax <- max(c(a3q, a3e, a4q, a4e), na.rm = TRUE) * 1.15
    for (obj in list(list(a3q, "A3 quantile"), list(a3e, "A3 expectile"), list(a4q, "A4 quantile"), list(a4e, "A4 expectile"))) {
      vals <- obj[[1]]; ttl <- obj[[2]]
      pat <- bar_pattern(length(vals))
      bp <- barplot(vals, las = 2, col = method_colors(names(vals)),
                    border = "black", density = pat$density, angle = pat$angle,
                    ylim = c(0, ymax), ylab = "Average detection rate", main = ttl)
      text(bp, vals, labels = sprintf("%.1f", 100 * vals), pos = 3, cex = 0.6)
    }
  })

  # 4. Gradual design: weight tradeoff (kept focused on CvM branch)
  q4 <- read_summary(results_root, "table4_gradual_quantile")
  e4 <- read_summary(results_root, "table4_gradual_expectile")
  w_order <- c("EARLY", "MID", "U", "LATE")
  det_q <- cbind(
    SSMS = sapply(w_order, function(w) mean(q4[[paste0("SSMS_CvM_", w, "_det_rate")]], na.rm = TRUE)),
    RSMS = sapply(w_order, function(w) mean(q4[[paste0("RSMS_CvM_", w, "_det_rate")]], na.rm = TRUE)),
    HAC  = sapply(w_order, function(w) mean(q4[[paste0("HAC_CvM_", w, "_det_rate")]], na.rm = TRUE))
  )
  del_q <- cbind(
    SSMS = sapply(w_order, function(w) mean(q4[[paste0("SSMS_CvM_", w, "_avg_delay")]], na.rm = TRUE)),
    RSMS = sapply(w_order, function(w) mean(q4[[paste0("RSMS_CvM_", w, "_avg_delay")]], na.rm = TRUE)),
    HAC  = sapply(w_order, function(w) mean(q4[[paste0("HAC_CvM_", w, "_avg_delay")]], na.rm = TRUE))
  )
  det_e <- cbind(
    SSMS = sapply(w_order, function(w) mean(e4[[paste0("SSMS_CvM_", w, "_det_rate")]], na.rm = TRUE)),
    RSMS = sapply(w_order, function(w) mean(e4[[paste0("RSMS_CvM_", w, "_det_rate")]], na.rm = TRUE)),
    HAC  = sapply(w_order, function(w) mean(e4[[paste0("HAC_CvM_", w, "_det_rate")]], na.rm = TRUE))
  )
  del_e <- cbind(
    SSMS = sapply(w_order, function(w) mean(e4[[paste0("SSMS_CvM_", w, "_avg_delay")]], na.rm = TRUE)),
    RSMS = sapply(w_order, function(w) mean(e4[[paste0("RSMS_CvM_", w, "_avg_delay")]], na.rm = TRUE)),
    HAC  = sapply(w_order, function(w) mean(e4[[paste0("HAC_CvM_", w, "_avg_delay")]], na.rm = TRUE))
  )
  rownames(det_q) <- rownames(del_q) <- rownames(det_e) <- rownames(del_e) <- c("Early", "Mid", "Uniform", "Late")
  plot_weight_tradeoff(det_q, det_e, del_q, del_e, file.path(outdir, "fig_gradual_weight_tradeoff"))

  # 5. Contamination design
  c1q <- read_summary(results_root, "appendix_c1_quantile")
  c1e <- read_summary(results_root, "appendix_c1_expectile")
  det_c1q <- avg_metric(c1q, "_det_rate", methods_q)
  det_c1e <- avg_metric(c1e, "_det_rate", methods_e)
  del_c1q <- avg_metric(c1q, "_avg_delay", methods_q)
  del_c1e <- avg_metric(c1e, "_avg_delay", methods_e)
  save_plot_both(file.path(outdir, "fig_contamination_detection"), width = ifelse(method_set == "all", 14, 9), expr_fun = function() {
    plot_bar_two_panel(det_c1q, det_c1e,
                       main_left = paste0("C1 contamination: quantile", title_suffix),
                       main_right = paste0("C1 contamination: expectile", title_suffix),
                       ylab = "Average detection rate",
                       ylim = c(0, max(c(det_c1q, det_c1e), na.rm = TRUE) * 1.15))
  })
  save_plot_both(file.path(outdir, "fig_contamination_delay"), width = ifelse(method_set == "all", 14, 9), expr_fun = function() {
    plot_bar_two_panel(del_c1q, del_c1e,
                       main_left = paste0("C1 contamination: quantile", title_suffix),
                       main_right = paste0("C1 contamination: expectile", title_suffix),
                       ylab = "Average delay",
                       ylim = c(0, max(c(del_c1q, del_c1e), na.rm = TRUE) * 1.15))
  })

  message("Figures written to: ", outdir)
  invisible(outdir)
}

args <- commandArgs(trailingOnly = TRUE)
if (!interactive()) {
  root <- if (length(args) >= 1L) args[[1L]] else "."
  method_set <- if (length(args) >= 2L) args[[2L]] else "paper"
  make_all_simulation_summary_figures(results_root = root, method_set = method_set)
}

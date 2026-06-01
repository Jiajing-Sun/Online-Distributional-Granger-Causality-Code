# Rolling-window monitoring for the 2025 BTC--Deribit empirical panel.

source(file.path(dirname(sys.frame(1)$ofile), "00_config.R"))
emp2025_require_packages(c("data.table", "quantreg"))
emp2025_source_monitor_code()

emp2025_scale_by_training <- function(M, train_idx) {
  M <- as.matrix(M)
  mu <- colMeans(M[train_idx, , drop = FALSE], na.rm = TRUE)
  sdv <- apply(M[train_idx, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv < 1e-10] <- 1
  out <- sweep(sweep(M, 2, mu, "-"), 2, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

emp2025_pc1_by_training <- function(Z, train_idx) {
  Zs <- emp2025_scale_by_training(Z, train_idx)
  tr <- Zs[train_idx, , drop = FALSE]
  sv <- svd(tr, nu = 0, nv = 1)
  v <- sv$v[, 1]
  as.numeric(Zs %*% v)
}

emp2025_build_score_matrix <- function(y, X, Z, m, taus, baseline_only = FALSE) {
  X <- as.matrix(cbind(`(Intercept)` = 1, X))
  y <- as.numeric(y)
  if (baseline_only) {
    H <- matrix(1, nrow = length(y), ncol = 1)
    colnames(H) <- "baseline_hit"
  } else {
    H <- as.matrix(Z)
  }
  fit <- fit_frozen_models(X_train = X[1:m, , drop = FALSE], y_train = y[1:m],
                           tau_grid = taus, model_type = "quantile",
                           prefer_quantreg = TRUE)
  q <- ncol(H) * length(taus)
  psi <- matrix(NA_real_, nrow = length(y), ncol = q)
  cn <- character(q)
  pos <- 1L
  for (j in seq_along(taus)) {
    tau <- taus[j]
    u <- y - as.vector(X %*% fit$coefficients[, j])
    sc <- tau - as.numeric(u <= 0)
    block <- H * as.numeric(sc)
    psi[, pos:(pos + ncol(H) - 1L)] <- block
    cn[pos:(pos + ncol(H) - 1L)] <- paste0(colnames(H), "_tau", format(tau, trim = TRUE))
    pos <- pos + ncol(H)
  }
  colnames(psi) <- cn
  list(psi = psi, coefficients = fit$coefficients, converged = fit$converged, X = X)
}

emp2025_method_grid <- function() {
  ks_methods <- c("SSMS_KS_g0", "SSMS_KS_g015", "RSMS_KS_g0",
                  "RSMS_KS_g015", "HAC_KS_g0", "HAC_KS_g015")
  cvm_methods <- as.vector(outer(c("SSMS", "RSMS", "HAC"),
                                 c("U", "Early", "Mid", "Late"),
                                 paste, sep = "_CvM_"))
  e_methods <- c("EProc_Mix", "EProc_MultiStart", "EProc_Adaptive")
  data.table::rbindlist(list(
    data.table::data.table(method = ks_methods,
                           family = sub("_.*", "", ks_methods), stat = "KS",
                           gamma = c(0, 0.15, 0, 0.15, 0, 0.15), weight = NA_character_),
    data.table::data.table(method = cvm_methods,
                           family = sub("_.*", "", cvm_methods), stat = "CvM",
                           gamma = NA_real_,
                           weight = rep(c("U", "Early", "Mid", "Late"), each = 3)),
    data.table::data.table(method = e_methods,
                           family = "EProc", stat = "E", gamma = NA_real_, weight = NA_character_)
  ), fill = TRUE)
}

emp2025_stat_max <- function(method, psi, m, T, cfg) {
  if (grepl("^SSMS_KS", method)) return(ssms_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), Inf)$max_stat)
  if (grepl("^RSMS_KS", method)) return(rsms_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), Inf)$max_stat)
  if (grepl("^HAC_KS", method)) return(hac_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), Inf)$max_stat)
  if (grepl("^SSMS_CvM", method)) return(ssms_cvm_stop(psi, m, T, sub("^SSMS_CvM_", "", method), Inf)$max_stat)
  if (grepl("^RSMS_CvM", method)) return(rsms_cvm_stop(psi, m, T, sub("^RSMS_CvM_", "", method), Inf)$max_stat)
  if (grepl("^HAC_CvM", method)) return(hac_cvm_stop(psi, m, T, sub("^HAC_CvM_", "", method), Inf)$max_stat)
  stop("Unknown method: ", method)
}

emp2025_calibrate_monitor_cvs <- function(q, m, monitor_size, methods, cfg) {
  cv_path <- file.path(cfg$output_dir, "cache", sprintf("critical_values_q%d_m%d_h%d.csv", q, m, monitor_size))
  B <- if (isTRUE(cfg$smoke)) min(25L, cfg$n_cv_sims) else cfg$n_cv_sims
  if (file.exists(cv_path) && !isTRUE(cfg$force)) {
    cv <- data.table::fread(cv_path)
    ok_methods <- all(methods %in% cv$method)
    ok_alpha <- "alpha" %in% names(cv) && all(abs(cv[method %in% methods, alpha] - cfg$alpha) < 1e-12)
    ok_sims <- "n_sims" %in% names(cv) && all(cv[method %in% methods, n_sims] >= B)
    if (ok_methods && ok_alpha && ok_sims) return(cv[method %in% methods])
  }
  set.seed(cfg$seed + q + monitor_size)
  T <- monitor_size / m
  if (!isTRUE(cfg$smoke) && B < 20000L) {
    message("Critical-value warning: n_cv_sims=", B,
            " is suitable for pilot runs but low for final empirical reporting.")
  }
  if (T < 1) {
    message("Critical-value note: empirical T_ratio=", signif(T, 4),
            " is below the manuscript's standard tabulated grid; using cached finite-sample iid-Gaussian MC.")
  }
  out <- vector("list", length(methods))
  for (ii in seq_along(methods)) {
    mm <- methods[ii]
    vals <- numeric(B)
    for (b in seq_len(B)) {
      psi <- matrix(stats::rnorm((m + monitor_size) * q), ncol = q)
      vals[b] <- emp2025_stat_max(mm, psi, m, T, cfg)
    }
    out[[ii]] <- data.table::data.table(
      method = mm, q = q, m = m, monitor_size = monitor_size, T_ratio = T,
      alpha = cfg$alpha, critical_value = as.numeric(stats::quantile(vals, 1 - cfg$alpha, names = FALSE)),
      n_sims = B,
      cv_source = "finite_sample_iid_gaussian_mc",
      cv_note = "Empirical T may be outside manuscript tabulated Brownian-CV grid; cache and validate before final reporting."
    )
  }
  cv <- data.table::rbindlist(out)
  old <- if (file.exists(cv_path)) data.table::fread(cv_path) else data.table::data.table()
  if (nrow(old)) {
    old <- old[!(old$method %in% methods & old$q == q & old$m == m & old$monitor_size == monitor_size)]
  }
  emp2025_write_csv(data.table::rbindlist(list(old, cv), fill = TRUE), cv_path)
  cv
}

emp2025_run_asymptotic_method <- function(method, psi, m, monitor_size, cv) {
  T <- monitor_size / m
  if (grepl("^SSMS_KS", method)) return(ssms_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), cv))
  if (grepl("^RSMS_KS", method)) return(rsms_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), cv))
  if (grepl("^HAC_KS", method)) return(hac_ks_stop(psi, m, T, ifelse(grepl("015", method), 0.15, 0), cv))
  if (grepl("^SSMS_CvM", method)) return(ssms_cvm_stop(psi, m, T, sub("^SSMS_CvM_", "", method), cv))
  if (grepl("^RSMS_CvM", method)) return(rsms_cvm_stop(psi, m, T, sub("^RSMS_CvM_", "", method), cv))
  if (grepl("^HAC_CvM", method)) return(hac_cvm_stop(psi, m, T, sub("^HAC_CvM_", "", method), cv))
  stop("Unknown asymptotic method: ", method)
}

emp2025_run_eprocess_method <- function(method, y, X, z_index, m, monitor_size, taus, cfg, coeff_mat) {
  X <- as.matrix(cbind(`(Intercept)` = 1, X))
  restart_grid <- 1L
  feature_type <- "z"
  if (method == "EProc_MultiStart") restart_grid <- seq.int(1L, monitor_size, by = cfg$eprocess_restart_every)
  if (method == "EProc_Adaptive") feature_type <- "zminus"
  run_quantile_eprocess(
    y = y, X = X, zlag = z_index, m = m, n_monitor = monitor_size,
    alpha = cfg$alpha, tau_grid = taus, coeff_mat = coeff_mat,
    feature_type = feature_type, theta_dict = cfg$theta_dict,
    restart_grid = restart_grid, return_path = TRUE
  )
}

emp2025_make_windows <- function(n, cfg) {
  empty <- data.table::data.table(
    window_id = integer(), train_start_idx = integer(), train_end_idx = integer(),
    monitor_start_idx = integer(), monitor_end_idx = integer()
  )
  needed <- cfg$training_size + cfg$monitor_size
  if (!is.finite(n) || n < needed) return(empty)
  starts <- seq.int(1L, n - needed + 1L, by = max(1L, cfg$refit_every))
  if (isTRUE(cfg$smoke)) starts <- head(starts, 3L)
  data.table::data.table(
    window_id = seq_along(starts),
    train_start_idx = starts,
    train_end_idx = starts + cfg$training_size - 1L,
    monitor_start_idx = starts + cfg$training_size,
    monitor_end_idx = starts + cfg$training_size + cfg$monitor_size - 1L
  )
}

emp2025_run_windows <- function(panel, cfg, z_override = NULL, label = "option") {
  x_cols <- c("ret_1", "ret_mom", "rv_ann_24", "range_1", "log_quote_volume")
  z_cols <- if (is.null(z_override)) c("iv_rv_spread", "skew_proxy", "term_slope", "put_call_imbalance", "activity_log") else names(z_override)
  Zfull <- if (is.null(z_override)) as.matrix(panel[, ..z_cols]) else as.matrix(z_override)
  Xfull <- as.matrix(panel[, ..x_cols])
  yfull <- panel$y
  methods <- emp2025_method_grid()
  asym_methods <- methods[stat != "E", method]
  wins <- emp2025_make_windows(nrow(panel), cfg)
  if (!nrow(wins)) {
    stop("No rolling windows available. Panel rows=", nrow(panel),
         ", required at least training_size + monitor_size = ",
         cfg$training_size + cfg$monitor_size, ".")
  }

  cv_main <- emp2025_calibrate_monitor_cvs(q = length(cfg$taus) * ncol(Zfull),
                                           m = cfg$training_size,
                                           monitor_size = cfg$monitor_size,
                                           methods = asym_methods, cfg = cfg)
  rows <- vector("list", nrow(wins) * nrow(methods))
  rr <- 1L
  for (ww in seq_len(nrow(wins))) {
    w <- wins[ww]
    idx <- w$train_start_idx:w$monitor_end_idx
    train_idx_local <- seq_len(cfg$training_size)
    X <- Xfull[idx, , drop = FALSE]
    Z <- emp2025_scale_by_training(Zfull[idx, , drop = FALSE], train_idx_local)
    y <- yfull[idx]
    sc <- emp2025_build_score_matrix(y, X, Z, cfg$training_size, cfg$taus)
    z_index <- emp2025_pc1_by_training(Z, train_idx_local)

    for (ii in seq_len(nrow(methods))) {
      mm <- methods$method[ii]
      if (methods$stat[ii] == "E") {
        fit <- emp2025_run_eprocess_method(mm, y, X, z_index, cfg$training_size,
                                           cfg$monitor_size, cfg$taus, cfg, sc$coefficients)
        stop_k <- fit$stop_k
        max_path <- max(fit$log_path, na.rm = TRUE)
        threshold <- fit$log_threshold
        reject <- isTRUE(fit$reject)
      } else {
        threshold <- cv_main[method == mm, critical_value][1]
        fit <- emp2025_run_asymptotic_method(mm, sc$psi, cfg$training_size,
                                             cfg$monitor_size, threshold)
        stop_k <- fit$stop_k
        max_path <- fit$max_stat
        reject <- isTRUE(fit$reject)
      }
      alarm_abs_idx <- if (reject) w$monitor_start_idx + stop_k - 1L else NA_integer_
      trade_start_idx <- if (reject) alarm_abs_idx + 1L else NA_integer_
      rows[[rr]] <- data.table::data.table(
        sample_label = label,
        method = mm, family = methods$family[ii], stat = methods$stat[ii],
        gamma = methods$gamma[ii], weight = methods$weight[ii],
        threshold = threshold, threshold_scale = ifelse(methods$stat[ii] == "E", "log_e", "statistic"),
        alarm_k = ifelse(reject, stop_k, NA_integer_), max_path = max_path,
        max_path_over_threshold = max_path / threshold,
        window_id = w$window_id,
        train_start_idx = w$train_start_idx, train_end_idx = w$train_end_idx,
        monitor_start_idx = w$monitor_start_idx, monitor_end_idx = w$monitor_end_idx,
        alarm_abs_idx = alarm_abs_idx, trade_start_idx = trade_start_idx,
        train_start_time = panel$hour_end[w$train_start_idx],
        train_end_time = panel$hour_end[w$train_end_idx],
        monitor_start_time = panel$hour_end[w$monitor_start_idx],
        monitor_end_time = panel$hour_end[w$monitor_end_idx],
        alarm_time = if (reject) panel$hour_end[alarm_abs_idx] else as.POSIXct(NA, tz = "UTC"),
        trade_start_time = if (reject && trade_start_idx <= nrow(panel)) panel$hour_end[trade_start_idx] else as.POSIXct(NA, tz = "UTC"),
        q_dim = ncol(sc$psi), T_ratio = cfg$monitor_size / cfg$training_size
      )
      rr <- rr + 1L
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

emp2025_run_baseline_monitor <- function(panel, cfg) {
  x_cols <- c("ret_1", "ret_mom", "rv_ann_24", "range_1", "log_quote_volume")
  Xfull <- as.matrix(panel[, ..x_cols])
  yfull <- panel$y
  methods <- emp2025_method_grid()[stat != "E"]
  wins <- emp2025_make_windows(nrow(panel), cfg)
  if (!nrow(wins)) {
    stop("No baseline-monitor windows available. Panel rows=", nrow(panel),
         ", required at least training_size + monitor_size = ",
         cfg$training_size + cfg$monitor_size, ".")
  }
  cv_base <- emp2025_calibrate_monitor_cvs(q = length(cfg$taus), m = cfg$training_size,
                                           monitor_size = cfg$monitor_size,
                                           methods = methods$method, cfg = cfg)
  rows <- list()
  rr <- 1L
  for (ww in seq_len(nrow(wins))) {
    w <- wins[ww]
    idx <- w$train_start_idx:w$monitor_end_idx
    sc <- emp2025_build_score_matrix(yfull[idx], Xfull[idx, , drop = FALSE], NULL,
                                     cfg$training_size, cfg$taus, baseline_only = TRUE)
    for (mm in methods$method) {
      threshold <- cv_base[method == mm, critical_value][1]
      fit <- emp2025_run_asymptotic_method(mm, sc$psi, cfg$training_size, cfg$monitor_size, threshold)
      reject <- isTRUE(fit$reject)
      alarm_abs_idx <- if (reject) w$monitor_start_idx + fit$stop_k - 1L else NA_integer_
      rows[[rr]] <- data.table::data.table(
        window_id = w$window_id, method = mm, statistic = methods[method == mm, stat][1],
        stat_cv_ratio = fit$max_stat / threshold, reject = reject,
        first_alarm_time = if (reject) panel$hour_end[alarm_abs_idx] else as.POSIXct(NA, tz = "UTC"),
        threshold = threshold, max_stat = fit$max_stat
      )
      rr <- rr + 1L
    }
  }
  data.table::rbindlist(rows)
}

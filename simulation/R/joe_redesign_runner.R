# ==============================================================
# joe_redesign_runner.R -- JoE-oriented size-adjusted simulation runner
# ==============================================================

source_project("R", "utils.R")
source_project("R", "critical_values.R")
source_project("R", "weights.R")
source_project("R", "monitors.R")
source_project("R", "fit_distributional_models.R")
source_project("R", "dgp_online_distgc.R")
source_project("R", "eprocess.R")

reserve_cores_n <- function(reserve = 3L, override = NULL) {
  if (!is.null(override) && is.finite(as.numeric(override))) {
    return(max(1L, as.integer(override)))
  }
  get_default_ncores(reserve = reserve)
}

make_tau_grid <- function(model_type = c("quantile", "expectile")) {
  model_type <- match.arg(model_type)
  if (model_type == "quantile") {
    return(c(0.05, 0.10, 0.50, 0.90, 0.95))
  }
  c(0.10, 0.25, 0.50, 0.75, 0.90)
}

tau_grid_id <- function(model_type = c("quantile", "expectile")) {
  model_type <- match.arg(model_type)
  if (model_type == "quantile") "Q_005_010_050_090_095" else "E_010_025_050_075_090"
}

make_tau_grid_override <- function(model_type = c("quantile", "expectile"),
                                   tau_grid_override = NULL) {
  model_type <- match.arg(model_type)
  if (is.null(tau_grid_override) || !nzchar(tau_grid_override) || tau_grid_override == "default") {
    return(make_tau_grid(model_type))
  }
  if (model_type != "quantile") return(make_tau_grid(model_type))
  switch(
    tau_grid_override,
    Q5_tail = c(0.05, 0.10, 0.50, 0.90, 0.95),
    Q3_midtail = c(0.10, 0.50, 0.90),
    Q3_tail = c(0.05, 0.50, 0.95),
    stop("Unknown tau_grid_override: ", tau_grid_override)
  )
}

tau_grid_id_override <- function(model_type = c("quantile", "expectile"),
                                 tau_grid_override = NULL) {
  model_type <- match.arg(model_type)
  if (is.null(tau_grid_override) || !nzchar(tau_grid_override) || tau_grid_override == "default") {
    return(tau_grid_id(model_type))
  }
  if (model_type != "quantile") return(tau_grid_id(model_type))
  switch(
    tau_grid_override,
    Q5_tail = "Q_005_010_050_090_095",
    Q3_midtail = "Q_010_050_090",
    Q3_tail = "Q_005_050_095",
    stop("Unknown tau_grid_override: ", tau_grid_override)
  )
}

dgp_variant_params <- function(dgp_variant_id = "baseline") {
  variants <- data.frame(
    dgp_variant_id = c("baseline", "baseline_ar05", "low_ar03_mild_garch",
                       "low_ar02_milder_garch", "iid_mild_garch"),
    phi_y = c(0.50, 0.50, 0.30, 0.20, 0.00),
    phi_z = c(0.50, 0.50, 0.30, 0.20, 0.00),
    n3_df = c(5, 5, 7, 8, 8),
    n4_phi_z = c(0.80, 0.80, 0.65, 0.60, 0.50),
    garch_alpha = c(0.05, 0.05, 0.04, 0.03, 0.03),
    garch_beta = c(0.85, 0.85, 0.80, 0.75, 0.70),
    stringsAsFactors = FALSE
  )
  variants$garch_omega <- 1 - variants$garch_alpha - variants$garch_beta
  idx <- match(dgp_variant_id, variants$dgp_variant_id)
  if (is.na(idx)) stop("Unknown dgp_variant_id: ", dgp_variant_id)
  variants[idx, , drop = FALSE]
}

instrument_description <- function(instrument_id) {
  switch(
    instrument_id,
    z = "H = Z",
    asym = "H = (Z, Zminus)",
    scale = "H = (Z, Z^2 - Phase-I mean)",
    rich = "H = (Z, Zminus, Zminus^2 centered, Z^2 centered)",
    instrument_id
  )
}

make_main_method_grid <- function(include_eprocess = TRUE) {
  out <- data.frame(
    method_id = c(
      "SSMS_KS_g0",
      "RSMS_KS_g0_LDL",
      "HAC_KS_g0_Bartlett_m13",
      "SSMS_CvM_Late",
      "RSMS_CvM_Late_LDL",
      "HAC_CvM_Late_Bartlett_m13"
    ),
    family = "asymptotic",
    standardizer = c("SSMS", "RSMS", "HAC", "SSMS", "RSMS", "HAC"),
    statistic = c("KS", "KS", "KS", "CvM", "CvM", "CvM"),
    gamma = c(0, 0, 0, 0, 0, 0),
    weight = c(NA, NA, NA, "Late", "Late", "Late"),
    hac_kernel = c(NA, NA, "bartlett", NA, NA, "bartlett"),
    hac_bw_rule = c(NA, NA, "m13", NA, NA, "m13"),
    hac_bw_mult = c(NA, NA, 1, NA, NA, 1),
    rsms_prewhitening = c(NA, "LDL_lag0", NA, NA, "LDL_lag0", NA),
    rsms_ridge = c(NA, 1e-8, NA, NA, 1e-8, NA),
    main_text = TRUE,
    stringsAsFactors = FALSE
  )

  if (include_eprocess) {
    ep <- data.frame(
      method_id = c("EPROC_MIX_default", "EPROC_BANK_dyadic"),
      family = "eprocess",
      standardizer = "EPROC",
      statistic = "logE",
      gamma = NA_real_,
      weight = NA_character_,
      hac_kernel = NA_character_,
      hac_bw_rule = NA_character_,
      hac_bw_mult = NA_real_,
      rsms_prewhitening = NA_character_,
      rsms_ridge = NA_real_,
      main_text = TRUE,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, ep)
  }

  out
}

build_design_grid <- function(package = c("pilot", "minimum", "full"),
                              include_expectile = FALSE,
                              tau_grid_override = NULL,
                              dgp_variant_id = "baseline",
                              include_n3 = NULL,
                              include_n4 = NULL,
                              m_override = NULL,
                              T_override = NULL,
                              c_override = NULL,
                              break_override = NULL) {
  package <- match.arg(package)
  if (is.null(include_n3)) include_n3 <- package %in% c("pilot", "full")
  if (is.null(include_n4)) include_n4 <- package %in% c("pilot", "full")
  dgp_params <- dgp_variant_params(dgp_variant_id)

  m_vec <- if (package == "pilot") 200 else c(200, 500)
  T_vec <- if (package == "pilot") 1 else c(1, 2, 5)
  c_vec <- if (package == "pilot") 0.25 else c(0.10, 0.25, 0.50)
  break_vec <- if (package == "pilot") 0.50 else c(0.25, 0.50, 0.75)
  if (!is.null(m_override)) m_vec <- as.integer(m_override)
  if (!is.null(T_override)) T_vec <- as.numeric(T_override)
  if (!is.null(c_override)) c_vec <- as.numeric(c_override)
  if (!is.null(break_override)) break_vec <- as.numeric(break_override)
  model_types <- if (include_expectile) c("quantile", "expectile") else "quantile"

  null_designs <- c("N1", "N2")
  if (include_n3) null_designs <- c(null_designs, "N3")
  if (include_n4) null_designs <- c(null_designs, "N4")

  null_rows <- list()
  ii <- 1L
  for (model_type in model_types) {
    tau_values <- make_tau_grid_override(model_type, tau_grid_override)
    tau_id <- tau_grid_id_override(model_type, tau_grid_override)
    for (nd in null_designs) {
      inst_vec <- if (nd == "N4") "scale" else c("z", "asym", "scale")
      for (inst in inst_vec) {
        for (m in m_vec) for (T0 in T_vec) {
          null_rows[[ii]] <- data.frame(
            design_id = nd,
            design_family = "null",
            null_match_id = nd,
            dgp_variant_id = dgp_params$dgp_variant_id,
            phi_y = dgp_params$phi_y,
            phi_z = dgp_params$phi_z,
            n3_df = dgp_params$n3_df,
            n4_phi_z = dgp_params$n4_phi_z,
            garch_omega = dgp_params$garch_omega,
            garch_alpha = dgp_params$garch_alpha,
            garch_beta = dgp_params$garch_beta,
            model_type = model_type,
            tau_grid_id = tau_id,
            tau_grid = paste(tau_values, collapse = ";"),
            tau_weight_scheme = "equal",
            instrument_id = inst,
            instrument_description = instrument_description(inst),
            q_expected = length(tau_values) * ifelse(inst == "z", 1L, ifelse(inst %in% c("asym", "scale"), 2L, 4L)),
            m = m,
            T = T0,
            n_monitor = as.integer(round(m * T0)),
            c_val = 0,
            break_frac = NA_real_,
            k_star = NA_integer_,
            include_eprocess = model_type == "quantile",
            eprocess_tau_grid_id = ifelse(model_type == "quantile", "Q_005_010", NA_character_),
            eprocess_feature_id = ifelse(inst == "asym", "zminus", "z"),
            notes = "null calibration/evaluation",
            stringsAsFactors = FALSE
          )
          ii <- ii + 1L
        }
      }
    }
  }
  null_df <- do.call(rbind, null_rows)

  alt_map <- data.frame(
    design_id = c("A1", "A2", "A3A", "A3B", "A4", "A5", "C1"),
    null_match_id = c("N1", "N2", "N1", "N1", "N1", "N1", "N1"),
    instrument_id = c("z", "z", "z", "scale", "asym", "z", "z"),
    notes = c(
      "abrupt location",
      "abrupt location with GARCH",
      "under-instrumented scale negative control",
      "enriched scale alternative",
      "downside-tail asymmetric instrument",
      "gradual onset",
      "training contamination robustness"
    ),
    stringsAsFactors = FALSE
  )

  alt_rows <- list()
  ii <- 1L
  for (model_type in model_types) {
    tau_values <- make_tau_grid_override(model_type, tau_grid_override)
    tau_id <- tau_grid_id_override(model_type, tau_grid_override)
    for (rr in seq_len(nrow(alt_map))) {
      for (m in m_vec) for (T0 in T_vec) for (cc in c_vec) for (bb in break_vec) {
        if (alt_map$design_id[rr] == "C1" && bb != break_vec[1]) next
        inst <- alt_map$instrument_id[rr]
        alt_rows[[ii]] <- data.frame(
          design_id = alt_map$design_id[rr],
          design_family = "alternative",
          null_match_id = alt_map$null_match_id[rr],
          dgp_variant_id = dgp_params$dgp_variant_id,
          phi_y = dgp_params$phi_y,
          phi_z = dgp_params$phi_z,
          n3_df = dgp_params$n3_df,
          n4_phi_z = dgp_params$n4_phi_z,
          garch_omega = dgp_params$garch_omega,
          garch_alpha = dgp_params$garch_alpha,
          garch_beta = dgp_params$garch_beta,
          model_type = model_type,
          tau_grid_id = tau_id,
          tau_grid = paste(tau_values, collapse = ";"),
          tau_weight_scheme = "equal",
          instrument_id = inst,
          instrument_description = instrument_description(inst),
          q_expected = length(tau_values) * ifelse(inst == "z", 1L, ifelse(inst %in% c("asym", "scale"), 2L, 4L)),
          m = m,
          T = T0,
          n_monitor = as.integer(round(m * T0)),
          c_val = cc,
          break_frac = bb,
          k_star = ifelse(alt_map$design_id[rr] == "C1", 0L, as.integer(floor(bb * m * T0))),
          include_eprocess = model_type == "quantile",
          eprocess_tau_grid_id = ifelse(model_type == "quantile", "Q_005_010", NA_character_),
          eprocess_feature_id = ifelse(inst == "asym", "zminus", "z"),
          notes = alt_map$notes[rr],
          stringsAsFactors = FALSE
        )
        ii <- ii + 1L
      }
    }
  }
  alt_df <- do.call(rbind, alt_rows)

  out <- rbind(null_df, alt_df)
  out$run_id <- seq_len(nrow(out))
  out[, c("run_id", setdiff(names(out), "run_id"))]
}

get_brownian_cv_for_method <- function(cv, method_row, T, q, alpha = 0.05) {
  if (identical(method_row$family, "eprocess")) return(-log(alpha))
  get_critical_value(
    cv = cv,
    stat = method_row$standardizer,
    type = method_row$statistic,
    T = T,
    q = q,
    alpha = alpha,
    gamma = method_row$gamma,
    weight = ifelse(is.na(method_row$weight), "U", method_row$weight)
  )
}

compute_all_stat_paths <- function(psi, m, T, method_grid,
                                   cv,
                                   alpha = 0.05) {
  method_grid <- method_grid[method_grid$family == "asymptotic", , drop = FALSE]
  out <- vector("list", nrow(method_grid))

  for (ii in seq_len(nrow(method_grid))) {
    meth <- method_grid[ii, , drop = FALSE]
    q <- ncol(psi)
    crit <- get_brownian_cv_for_method(cv, meth, T = T, q = q, alpha = alpha)

    delta_off <- NA_real_
    lrv_offdiag_maxcorr <- NA_real_
    sigma0_condition <- NA_real_
    min_ldl_diag <- NA_real_
    hac_bandwidth <- NA_integer_

    if (meth$standardizer == "SSMS") {
      seq <- compute_ssms_sequence(psi = psi, m = m, T = T, gamma = meth$gamma)
      if (meth$statistic == "KS") {
        stat_path <- seq$M
      } else {
        k_vec <- seq_along(seq$M)
        w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = meth$weight)
        stat_path <- cumsum(w_vec * seq$M) / m
      }
    } else if (meth$standardizer == "RSMS") {
      seq <- compute_rsms_sequence(psi = psi, m = m, T = T, gamma = meth$gamma,
                                   ridge = ifelse(is.na(meth$rsms_ridge), 1e-8, meth$rsms_ridge))
      delta_off <- seq$delta_off
      lrv_offdiag_maxcorr <- seq$lrv_offdiag_maxcorr
      sigma0_condition <- seq$sigma0_condition
      min_ldl_diag <- seq$min_ldl_diag
      hac_bandwidth <- seq$rsms_hac_bandwidth
      if (meth$statistic == "KS") {
        stat_path <- seq$M
      } else {
        k_vec <- seq_along(seq$M)
        w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = meth$weight)
        stat_path <- cumsum(w_vec * seq$M) / m
      }
    } else if (meth$standardizer == "HAC") {
      seq <- compute_hac_sequence(
        psi = psi,
        m = m,
        T = T,
        gamma = meth$gamma,
        kernel = meth$hac_kernel,
        bw_rule = meth$hac_bw_rule,
        bw_const = meth$hac_bw_mult
      )
      hac_bandwidth <- seq$h
      if (meth$statistic == "KS") {
        stat_path <- seq$M
      } else {
        k_vec <- seq_along(seq$M)
        w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = meth$weight)
        stat_path <- cumsum(w_vec * seq$M) / m
      }
    } else {
      stop("Unknown standardizer: ", meth$standardizer)
    }

    out[[ii]] <- list(
      method = meth,
      stat_path = as.numeric(stat_path),
      max_stat = max(stat_path, na.rm = TRUE),
      brownian_cv = crit,
      delta_off = delta_off,
      lrv_offdiag_maxcorr = lrv_offdiag_maxcorr,
      sigma0_condition = sigma0_condition,
      min_ldl_diag = min_ldl_diag,
      hac_bandwidth = hac_bandwidth
    )
  }

  out
}

apply_threshold_to_path <- function(stat_path, critical_value, n_monitor,
                                    k_star = NA_integer_, design_id = NULL) {
  stat_path <- as.numeric(stat_path)
  stop_time <- which(stat_path > critical_value)[1]
  if (is.na(stop_time)) stop_time <- n_monitor + 1L

  overall <- stop_time <= n_monitor
  if (is.na(k_star)) {
    return(list(
      stop_time = stop_time,
      overall_reject = overall,
      prebreak_false_alarm = overall,
      postbreak_detect = NA,
      delay = NA_real_,
      relative_delay = NA_real_,
      arl0 = min(stop_time, n_monitor)
    ))
  }

  if (!is.null(design_id) && toupper(design_id) == "C1") {
    delay <- if (overall) stop_time else NA_real_
    return(list(
      stop_time = stop_time,
      overall_reject = overall,
      prebreak_false_alarm = NA,
      postbreak_detect = overall,
      delay = delay,
      relative_delay = if (overall) delay / n_monitor else NA_real_,
      arl0 = NA_real_
    ))
  }

  fa <- stop_time <= k_star
  det <- stop_time > k_star && stop_time <= n_monitor
  delay <- if (det) stop_time - k_star else NA_real_
  denom <- max(1, n_monitor - k_star)

  list(
    stop_time = stop_time,
    overall_reject = overall,
    prebreak_false_alarm = fa,
    postbreak_detect = det,
    delay = delay,
    relative_delay = if (det) delay / denom else NA_real_,
    arl0 = NA_real_
  )
}

cv_match_cols <- function() {
  c("design_id", "model_type", "tau_grid_id", "instrument_id",
    "m", "T", "q", "method_id")
}

make_group_key <- function(df, keys) {
  pieces <- lapply(keys, function(kk) {
    x <- df[[kk]]
    ifelse(is.na(x), "NA", as.character(x))
  })
  do.call(paste, c(pieces, sep = "|"))
}

conf_num <- function(conf, name, default) {
  if (!name %in% names(conf)) return(default)
  val <- suppressWarnings(as.numeric(conf[[name]][1]))
  if (!is.finite(val)) default else val
}

conf_chr <- function(conf, name, default) {
  if (!name %in% names(conf)) return(default)
  val <- as.character(conf[[name]][1])
  if (length(val) == 0L || is.na(val) || !nzchar(val)) default else val
}

lookup_empirical_cv <- function(cv_map, conf, q, method_id) {
  if (is.null(cv_map) || nrow(cv_map) == 0L) return(NA_real_)
  sub <- cv_map[
    cv_map$design_id == conf$null_match_id &
      cv_map$model_type == conf$model_type &
      cv_map$tau_grid_id == conf$tau_grid_id &
      cv_map$instrument_id == conf$instrument_id &
      cv_map$m == conf$m &
      cv_map$T == conf$T &
      cv_map$q == q &
      cv_map$method_id == method_id,
    , drop = FALSE
  ]
  if ("dgp_variant_id" %in% names(cv_map)) {
    sub <- sub[sub$dgp_variant_id == conf_chr(conf, "dgp_variant_id", "baseline"), , drop = FALSE]
  }
  if (nrow(sub) != 1L) return(NA_real_)
  sub$empirical_null_cv[1]
}

run_one_joe_replication <- function(conf,
                                    rep_id,
                                    method_grid,
                                    cv,
                                    cv_map = NULL,
                                    alpha = 0.05,
                                    seed_base = 20260528,
                                    phase = c("calibration", "evaluation", "alternative"),
                                    prefer_quantreg = TRUE) {
  phase <- match.arg(phase)
  seed <- as.integer(seed_base + 100000L * conf$run_id + rep_id)
  set.seed(seed)

  tau_grid <- as.numeric(strsplit(conf$tau_grid, ";", fixed = TRUE)[[1]])

  sim <- simulate_online_distgc(
    design = conf$design_id,
    m = conf$m,
    T = conf$T,
    c_val = conf$c_val,
    break_frac = ifelse(is.na(conf$break_frac), 0.50, conf$break_frac),
    phi_y = conf_num(conf, "phi_y", 0.5),
    phi_z = conf_num(conf, "phi_z", 0.5),
    n3_df = conf_num(conf, "n3_df", 5),
    n4_phi_z = conf_num(conf, "n4_phi_z", 0.8),
    garch_omega = conf_num(conf, "garch_omega", 0.10),
    garch_alpha = conf_num(conf, "garch_alpha", 0.05),
    garch_beta = conf_num(conf, "garch_beta", 0.85)
  )

  score_obj <- score_matrix_from_frozen(
    y = sim$y,
    X = sim$X,
    zlag = sim$zlag,
    m = conf$m,
    tau_grid = tau_grid,
    model_type = conf$model_type,
    instrument_type = conf$instrument_id,
    tau_weight_scheme = conf$tau_weight_scheme,
    prefer_quantreg = prefer_quantreg
  )

  q <- ncol(score_obj$psi)
  paths <- compute_all_stat_paths(
    psi = score_obj$psi,
    m = conf$m,
    T = conf$T,
    method_grid = method_grid,
    cv = cv,
    alpha = alpha
  )

  rows <- list()
  rr <- 1L
  n_monitor <- sim$n_monitor
  k_star <- if (is.na(sim$k_star)) NA_integer_ else sim$k_star

  for (pp in paths) {
    meth <- pp$method
    brown_metrics <- apply_threshold_to_path(
      pp$stat_path,
      critical_value = pp$brownian_cv,
      n_monitor = n_monitor,
      k_star = k_star,
      design_id = conf$design_id
    )

    rows[[rr]] <- make_replication_row(
      conf = conf, rep_id = rep_id, seed = seed, q = q, phase = phase,
      method = meth, threshold_type = "brownian", critical_value = pp$brownian_cv,
      max_stat = pp$max_stat, metrics = brown_metrics,
      diagnostics = pp, all_fits_converged = all(score_obj$converged)
    )
    rr <- rr + 1L

    emp_cv <- lookup_empirical_cv(cv_map, conf = conf, q = q, method_id = meth$method_id)
    if (is.finite(emp_cv) && phase != "calibration") {
      emp_metrics <- apply_threshold_to_path(
        pp$stat_path,
        critical_value = emp_cv,
        n_monitor = n_monitor,
        k_star = k_star,
        design_id = conf$design_id
      )
      rows[[rr]] <- make_replication_row(
        conf = conf, rep_id = rep_id, seed = seed, q = q, phase = phase,
        method = meth, threshold_type = "empirical", critical_value = emp_cv,
        max_stat = pp$max_stat, metrics = emp_metrics,
        diagnostics = pp, all_fits_converged = all(score_obj$converged)
      )
      rr <- rr + 1L
    }
  }

  if (isTRUE(conf$include_eprocess) && conf$model_type == "quantile") {
    e_tau <- c(0.05, 0.10)
    feature <- conf$eprocess_feature_id

    ep_single <- run_quantile_eprocess(
      y = sim$y, X = sim$X, zlag = sim$zlag,
      m = conf$m, n_monitor = n_monitor, alpha = alpha,
      tau_grid = e_tau, feature_type = feature,
      restart_grid = 1L, prefer_quantreg = prefer_quantreg,
      return_path = TRUE
    )

    ep_bank <- run_quantile_eprocess(
      y = sim$y, X = sim$X, zlag = sim$zlag,
      m = conf$m, n_monitor = n_monitor, alpha = alpha,
      tau_grid = e_tau, feature_type = feature,
      restart_grid = make_restart_grid(n_monitor, type = "dyadic"),
      prefer_quantreg = prefer_quantreg,
      return_path = TRUE
    )

    for (ep in list(
      list(id = "EPROC_MIX_default", obj = ep_single),
      list(id = "EPROC_BANK_dyadic", obj = ep_bank)
    )) {
      meth <- method_grid[method_grid$method_id == ep$id, , drop = FALSE]
      if (nrow(meth) != 1L) next
      stat_path <- ep$obj$log_path
      crit <- -log(alpha)
      metrics <- apply_threshold_to_path(
        stat_path = stat_path,
        critical_value = crit,
        n_monitor = n_monitor,
        k_star = k_star,
        design_id = conf$design_id
      )
      rows[[rr]] <- make_replication_row(
        conf = conf, rep_id = rep_id, seed = seed, q = q, phase = phase,
        method = meth, threshold_type = "eprocess_fixed", critical_value = crit,
        max_stat = max(stat_path, na.rm = TRUE), metrics = metrics,
        diagnostics = list(delta_off = NA_real_, lrv_offdiag_maxcorr = NA_real_,
                           sigma0_condition = NA_real_, min_ldl_diag = NA_real_,
                           hac_bandwidth = NA_integer_),
        all_fits_converged = all(score_obj$converged)
      )
      rr <- rr + 1L
    }
  }

  do.call(rbind, rows)
}

make_replication_row <- function(conf, rep_id, seed, q, phase, method,
                                 threshold_type, critical_value, max_stat,
                                 metrics, diagnostics, all_fits_converged) {
  data.frame(
    phase = phase,
    rep_id = rep_id,
    seed = seed,
    run_id = conf$run_id,
    design_id = conf$design_id,
    design_family = conf$design_family,
    null_match_id = conf$null_match_id,
    dgp_variant_id = conf_chr(conf, "dgp_variant_id", "baseline"),
    phi_y = conf_num(conf, "phi_y", 0.5),
    phi_z = conf_num(conf, "phi_z", 0.5),
    n3_df = conf_num(conf, "n3_df", 5),
    n4_phi_z = conf_num(conf, "n4_phi_z", 0.8),
    garch_omega = conf_num(conf, "garch_omega", 0.10),
    garch_alpha = conf_num(conf, "garch_alpha", 0.05),
    garch_beta = conf_num(conf, "garch_beta", 0.85),
    model_type = conf$model_type,
    tau_grid_id = conf$tau_grid_id,
    instrument_id = conf$instrument_id,
    m = conf$m,
    T = conf$T,
    n_monitor = conf$n_monitor,
    q = q,
    c_val = conf$c_val,
    break_frac = conf$break_frac,
    k_star = conf$k_star,
    method_id = method$method_id,
    family = method$family,
    standardizer = method$standardizer,
    statistic = method$statistic,
    gamma = method$gamma,
    weight = method$weight,
    hac_kernel = method$hac_kernel,
    hac_bw_rule = method$hac_bw_rule,
    hac_bw_mult = method$hac_bw_mult,
    hac_bandwidth = diagnostics$hac_bandwidth,
    rsms_prewhitening = method$rsms_prewhitening,
    rsms_ridge = method$rsms_ridge,
    threshold_type = threshold_type,
    critical_value = critical_value,
    max_stat = max_stat,
    stop_time = metrics$stop_time,
    overall_reject = metrics$overall_reject,
    prebreak_false_alarm = metrics$prebreak_false_alarm,
    postbreak_detect = metrics$postbreak_detect,
    delay = metrics$delay,
    relative_delay = metrics$relative_delay,
    arl0 = metrics$arl0,
    delta_off = diagnostics$delta_off,
    lrv_offdiag_maxcorr = diagnostics$lrv_offdiag_maxcorr,
    sigma0_condition = diagnostics$sigma0_condition,
    min_ldl_diag = diagnostics$min_ldl_diag,
    all_fits_converged = all_fits_converged,
    stringsAsFactors = FALSE
  )
}

read_csv_if_exists <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}

append_csv <- function(x, path) {
  dir_create(dirname(path))
  has_file <- file.exists(path) && file.info(path)$size > 0
  utils::write.table(
    x,
    file = path,
    sep = ",",
    row.names = FALSE,
    col.names = !has_file,
    append = has_file,
    quote = TRUE,
    na = ""
  )
  invisible(path)
}

job_keys_from_rows <- function(df) {
  if (is.null(df) || nrow(df) == 0L || !all(c("run_id", "rep_id") %in% names(df))) {
    return(character(0))
  }
  paste(df$run_id, df$rep_id, sep = "|")
}

run_replication_grid <- function(grid, B, method_grid, cv, cv_map = NULL,
                                 alpha = 0.05, seed_base = 20260528,
                                 phase = c("calibration", "evaluation", "alternative"),
                                 ncores = 1L, prefer_quantreg = TRUE,
                                 output_path = NULL,
                                 resume = TRUE,
                                 checkpoint_batch_size = NULL) {
  phase <- match.arg(phase)
  grid <- as.data.frame(grid, stringsAsFactors = FALSE)
  jobs <- expand.grid(row_id = seq_len(nrow(grid)), rep_id = seq_len(B))
  jobs$run_id <- grid$run_id[jobs$row_id]
  jobs$job_key <- paste(jobs$run_id, jobs$rep_id, sep = "|")

  if (!is.null(output_path) && !isTRUE(resume) && file.exists(output_path)) {
    unlink(output_path)
  }

  existing <- if (!is.null(output_path) && isTRUE(resume)) read_csv_if_exists(output_path) else NULL
  done_keys <- job_keys_from_rows(existing)
  if (length(done_keys) > 0L) {
    jobs <- jobs[!jobs$job_key %in% done_keys, , drop = FALSE]
    message(sprintf("Resume enabled for %s: found %s completed jobs, %s remaining.",
                    basename(output_path), length(unique(done_keys)), nrow(jobs)))
  }

  if (nrow(jobs) == 0L) {
    if (!is.null(existing)) return(existing)
    return(data.frame())
  }

  if (is.null(checkpoint_batch_size)) {
    checkpoint_batch_size <- max(1L, as.integer(max(1L, ncores) * 4L))
  }
  checkpoint_batch_size <- max(1L, as.integer(checkpoint_batch_size))

  project_root <- get_project_root()
  cl <- make_cluster(ncores)
  on.exit(stop_cluster(cl), add = TRUE)

  if (!is.null(cl)) {
    parallel::clusterExport(
      cl,
      varlist = c("project_root"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      setwd(project_root)
      options(distgc.project_root = project_root)
      source(file.path(project_root, "R", "utils.R"))
      source(file.path(project_root, "R", "critical_values.R"))
      source(file.path(project_root, "R", "weights.R"))
      source(file.path(project_root, "R", "monitors.R"))
      source(file.path(project_root, "R", "fit_distributional_models.R"))
      source(file.path(project_root, "R", "dgp_online_distgc.R"))
      source(file.path(project_root, "R", "eprocess.R"))
      source(file.path(project_root, "R", "joe_redesign_runner.R"))
      NULL
    })
  }

  batches <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / checkpoint_batch_size))
  collected <- list()
  cc <- 1L

  for (bb in seq_along(batches)) {
    batch_jobs <- jobs[batches[[bb]], , drop = FALSE]
    message(sprintf("  %s batch %s/%s: %s jobs", phase, bb, length(batches), nrow(batch_jobs)))

    if (!is.null(cl)) {
      parallel::clusterExport(
        cl,
        varlist = c("grid", "batch_jobs", "method_grid", "cv", "cv_map", "alpha",
                    "seed_base", "phase", "prefer_quantreg"),
        envir = environment()
      )
      out_batch <- parallel::parLapply(cl, seq_len(nrow(batch_jobs)), function(jj) {
        conf <- grid[batch_jobs$row_id[jj], , drop = FALSE]
        run_one_joe_replication(
          conf = conf,
          rep_id = batch_jobs$rep_id[jj],
          method_grid = method_grid,
          cv = cv,
          cv_map = cv_map,
          alpha = alpha,
          seed_base = seed_base,
          phase = phase,
          prefer_quantreg = prefer_quantreg
        )
      })
    } else {
      out_batch <- lapply(seq_len(nrow(batch_jobs)), function(jj) {
        conf <- grid[batch_jobs$row_id[jj], , drop = FALSE]
        run_one_joe_replication(
          conf = conf,
          rep_id = batch_jobs$rep_id[jj],
          method_grid = method_grid,
          cv = cv,
          cv_map = cv_map,
          alpha = alpha,
          seed_base = seed_base,
          phase = phase,
          prefer_quantreg = prefer_quantreg
        )
      })
    }

    batch_df <- do.call(rbind, out_batch)
    if (!is.null(output_path)) {
      append_csv(batch_df, output_path)
    } else {
      collected[[cc]] <- batch_df
      cc <- cc + 1L
    }
  }

  if (!is.null(output_path)) {
    return(read.csv(output_path, stringsAsFactors = FALSE))
  }
  do.call(rbind, collected)
}

calibrate_null_cv <- function(null_stats_long, alpha = 0.05) {
  df <- null_stats_long[
    null_stats_long$family == "asymptotic" &
      null_stats_long$threshold_type == "brownian",
    , drop = FALSE
  ]

  keys <- c("design_id", "model_type", "tau_grid_id", "instrument_id", "m", "T",
            "q", "dgp_variant_id", "method_id", "standardizer", "statistic", "gamma", "weight",
            "hac_kernel", "hac_bw_rule", "hac_bw_mult", "rsms_prewhitening", "rsms_ridge")
  key <- make_group_key(df, keys)
  groups <- split(seq_len(nrow(df)), key)

  out <- lapply(groups, function(idx) {
    sub <- df[idx, , drop = FALSE]
    cv_emp <- as.numeric(stats::quantile(sub$max_stat, probs = 1 - alpha, names = FALSE, type = 8, na.rm = TRUE))
    data.frame(
      null_match_id = sub$design_id[1],
      design_id = sub$design_id[1],
      dgp_variant_id = sub$dgp_variant_id[1],
      phi_y = sub$phi_y[1],
      phi_z = sub$phi_z[1],
      n3_df = sub$n3_df[1],
      n4_phi_z = sub$n4_phi_z[1],
      garch_omega = sub$garch_omega[1],
      garch_alpha = sub$garch_alpha[1],
      garch_beta = sub$garch_beta[1],
      model_type = sub$model_type[1],
      tau_grid_id = sub$tau_grid_id[1],
      instrument_id = sub$instrument_id[1],
      m = sub$m[1],
      T = sub$T[1],
      q = sub$q[1],
      method_id = sub$method_id[1],
      standardizer = sub$standardizer[1],
      statistic = sub$statistic[1],
      gamma = sub$gamma[1],
      weight = sub$weight[1],
      hac_kernel = sub$hac_kernel[1],
      hac_bw_rule = sub$hac_bw_rule[1],
      hac_bw_mult = sub$hac_bw_mult[1],
      hac_bandwidth = suppressWarnings(stats::median(sub$hac_bandwidth, na.rm = TRUE)),
      rsms_prewhitening = sub$rsms_prewhitening[1],
      rsms_ridge = sub$rsms_ridge[1],
      nominal_alpha = alpha,
      brownian_cv = sub$critical_value[1],
      empirical_null_cv = cv_emp,
      empirical_cv_se_bootstrap = NA_real_,
      B_cal = nrow(sub),
      calibration_max_stat_mean = mean(sub$max_stat, na.rm = TRUE),
      calibration_max_stat_sd = stats::sd(sub$max_stat, na.rm = TRUE),
      calibration_size_brownian = mean(sub$overall_reject, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

summarise_size <- function(eval_long) {
  df <- eval_long[eval_long$design_family == "null", , drop = FALSE]
  keys <- c("design_id", "dgp_variant_id", "model_type", "tau_grid_id", "instrument_id", "m", "T",
            "q", "method_id", "threshold_type")
  key <- make_group_key(df, keys)
  groups <- split(seq_len(nrow(df)), key)

  out <- lapply(groups, function(idx) {
    sub <- df[idx, , drop = FALSE]
    p <- mean(sub$overall_reject, na.rm = TRUE)
    B <- sum(!is.na(sub$overall_reject))
    se <- sqrt(p * (1 - p) / max(1, B))
    data.frame(
      design_id = sub$design_id[1],
      dgp_variant_id = sub$dgp_variant_id[1],
      phi_y = sub$phi_y[1],
      phi_z = sub$phi_z[1],
      n3_df = sub$n3_df[1],
      n4_phi_z = sub$n4_phi_z[1],
      garch_omega = sub$garch_omega[1],
      garch_alpha = sub$garch_alpha[1],
      garch_beta = sub$garch_beta[1],
      model_type = sub$model_type[1],
      tau_grid_id = sub$tau_grid_id[1],
      instrument_id = sub$instrument_id[1],
      m = sub$m[1],
      T = sub$T[1],
      q = sub$q[1],
      method_id = sub$method_id[1],
      threshold_type = sub$threshold_type[1],
      critical_value = sub$critical_value[1],
      empirical_size = p,
      mc_se_size = se,
      ci95_low = max(0, p - 1.96 * se),
      ci95_high = min(1, p + 1.96 * se),
      mean_stop_time_if_reject = mean(sub$stop_time[sub$overall_reject], na.rm = TRUE),
      arl0_mean = mean(sub$arl0, na.rm = TRUE),
      B_eval = B,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

summarise_detection <- function(alt_long) {
  df <- alt_long[alt_long$design_family == "alternative", , drop = FALSE]
  keys <- c("design_id", "design_family", "null_match_id", "dgp_variant_id", "model_type",
            "tau_grid_id", "instrument_id", "m", "T", "c_val", "break_frac",
            "q", "method_id", "threshold_type")
  key <- make_group_key(df, keys)
  groups <- split(seq_len(nrow(df)), key)

  out <- lapply(groups, function(idx) {
    sub <- df[idx, , drop = FALSE]
    B <- nrow(sub)
    det <- mean(sub$postbreak_detect, na.rm = TRUE)
    fa <- mean(sub$prebreak_false_alarm, na.rm = TRUE)
    rej <- mean(sub$overall_reject, na.rm = TRUE)
    data.frame(
      design_id = sub$design_id[1],
      design_family = sub$design_family[1],
      null_match_id = sub$null_match_id[1],
      dgp_variant_id = sub$dgp_variant_id[1],
      phi_y = sub$phi_y[1],
      phi_z = sub$phi_z[1],
      n3_df = sub$n3_df[1],
      n4_phi_z = sub$n4_phi_z[1],
      garch_omega = sub$garch_omega[1],
      garch_alpha = sub$garch_alpha[1],
      garch_beta = sub$garch_beta[1],
      model_type = sub$model_type[1],
      tau_grid_id = sub$tau_grid_id[1],
      instrument_id = sub$instrument_id[1],
      m = sub$m[1],
      T = sub$T[1],
      c_val = sub$c_val[1],
      break_frac = sub$break_frac[1],
      q = sub$q[1],
      method_id = sub$method_id[1],
      threshold_type = sub$threshold_type[1],
      overall_rejection = rej,
      prebreak_false_alarm = fa,
      detection_probability = det,
      ADD = mean(sub$delay, na.rm = TRUE),
      median_delay = stats::median(sub$delay, na.rm = TRUE),
      rADD = mean(sub$relative_delay, na.rm = TRUE),
      mc_se_detection = sqrt(det * (1 - det) / max(1, B)),
      B_alt = B,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

summarise_rsms_offdiag <- function(df) {
  sub <- df[df$standardizer == "RSMS" & is.finite(df$delta_off), , drop = FALSE]
  if (nrow(sub) == 0L) return(data.frame())

  keys <- c("phase", "design_id", "dgp_variant_id", "model_type", "instrument_id", "m", "T", "q", "method_id")
  key <- make_group_key(sub, keys)
  groups <- split(seq_len(nrow(sub)), key)
  out <- lapply(groups, function(idx) {
    x <- sub[idx, , drop = FALSE]
    data.frame(
      phase = x$phase[1],
      design_id = x$design_id[1],
      dgp_variant_id = x$dgp_variant_id[1],
      phi_y = x$phi_y[1],
      phi_z = x$phi_z[1],
      n3_df = x$n3_df[1],
      n4_phi_z = x$n4_phi_z[1],
      garch_omega = x$garch_omega[1],
      garch_alpha = x$garch_alpha[1],
      garch_beta = x$garch_beta[1],
      model_type = x$model_type[1],
      instrument_id = x$instrument_id[1],
      m = x$m[1],
      T = x$T[1],
      q = x$q[1],
      method_id = x$method_id[1],
      median_delta_off = stats::median(x$delta_off, na.rm = TRUE),
      p90_delta_off = as.numeric(stats::quantile(x$delta_off, 0.90, na.rm = TRUE, names = FALSE)),
      p95_delta_off = as.numeric(stats::quantile(x$delta_off, 0.95, na.rm = TRUE, names = FALSE)),
      median_lrv_offdiag_maxcorr = stats::median(x$lrv_offdiag_maxcorr, na.rm = TRUE),
      B = nrow(x),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

run_oracle_eprocess_from_hits <- function(hits, s_mon,
                                          alpha = 0.05,
                                          tau_grid = c(0.05, 0.10),
                                          theta_dict = c(-1, -0.5, 0.5, 1),
                                          theta_weights = NULL,
                                          restart_grid = 1L,
                                          restart_weights = NULL) {
  hits <- as.matrix(hits)
  s_mon <- as.numeric(s_mon)
  n_monitor <- nrow(hits)
  tau_grid <- as.numeric(tau_grid)
  theta_dict <- as.numeric(theta_dict)

  if (is.null(theta_weights)) theta_weights <- rep(1, length(theta_dict))
  theta_weights <- normalize_weights(theta_weights)
  restart_grid <- as.integer(sort(unique(restart_grid[restart_grid >= 1 & restart_grid <= n_monitor])))
  if (length(restart_grid) == 0L) restart_grid <- 1L
  if (is.null(restart_weights)) restart_weights <- 1 / restart_grid
  restart_weights <- normalize_weights(restart_weights)
  tau_weights <- normalize_weights(rep(1, length(tau_grid)))

  J <- length(tau_grid)
  L <- length(theta_dict)
  Rn <- length(restart_grid)
  log_theta_w <- log(theta_weights)
  log_restart_w <- log(restart_weights)
  log_tau_w <- log(tau_weights)
  log_components <- array(NA_real_, dim = c(J, L, Rn))
  threshold_log <- -log(alpha)
  log_path <- rep(NA_real_, n_monitor)
  stop_k <- n_monitor + 1L

  for (k in seq_len(n_monitor)) {
    new_idx <- which(restart_grid == k)
    if (length(new_idx) > 0L) {
      for (rr in new_idx) log_components[, , rr] <- 0
    }

    log_terms_all <- c()
    for (j in seq_len(J)) {
      tau <- tau_grid[j]
      I_k <- hits[k, j]
      logL_vec <- numeric(L)
      for (ell in seq_len(L)) {
        lp <- stats::qlogis(tau) + theta_dict[ell] * s_mon[k]
        p_t <- min(max(stats::plogis(lp), 1e-8), 1 - 1e-8)
        logL_vec[ell] <- I_k * log(p_t / tau) + (1 - I_k) * log((1 - p_t) / (1 - tau))
      }
      for (rr in seq_len(Rn)) {
        if (all(is.na(log_components[j, , rr]))) next
        log_components[j, , rr] <- log_components[j, , rr] + logL_vec
        log_terms_all <- c(
          log_terms_all,
          log_tau_w[j] + log_restart_w[rr] + log_theta_w + log_components[j, , rr]
        )
      }
    }

    log_path[k] <- log_sum_exp(log_terms_all)
    if (stop_k > n_monitor && is.finite(log_path[k]) && log_path[k] > threshold_log) {
      stop_k <- k
    }
  }

  list(stop_k = stop_k, log_path = log_path, log_threshold = threshold_log)
}

run_one_oracle_eprocess_replication <- function(rep_id, m, T,
                                                alpha = 0.05,
                                                tau_grid = c(0.05, 0.10),
                                                seed_base = 20260528) {
  seed <- as.integer(seed_base + 700000000L + 100000L * m + 1000L * round(10 * T) + rep_id)
  set.seed(seed)

  n_monitor <- as.integer(round(m * T))
  n <- m + n_monitor
  z <- numeric(n)
  z_prev <- 0
  for (tt in seq_len(n)) {
    z[tt] <- z_prev
    z_prev <- 0.5 * z_prev + stats::rnorm(1)
  }
  s <- make_eprocess_feature(zlag = z, feature_type = "z", m = m)
  s_mon <- s[(m + 1):(m + n_monitor)]

  hits <- matrix(NA_real_, nrow = n_monitor, ncol = length(tau_grid))
  for (j in seq_along(tau_grid)) hits[, j] <- stats::rbinom(n_monitor, 1, tau_grid[j])

  single <- run_oracle_eprocess_from_hits(
    hits = hits,
    s_mon = s_mon,
    alpha = alpha,
    tau_grid = tau_grid,
    restart_grid = 1L
  )
  bank <- run_oracle_eprocess_from_hits(
    hits = hits,
    s_mon = s_mon,
    alpha = alpha,
    tau_grid = tau_grid,
    restart_grid = make_restart_grid(n_monitor, type = "dyadic")
  )

  rows <- list()
  ii <- 1L
  for (ep in list(
    list(id = "EPROC_MIX_default", obj = single),
    list(id = "EPROC_BANK_dyadic", obj = bank)
  )) {
    metrics <- apply_threshold_to_path(
      ep$obj$log_path,
      critical_value = -log(alpha),
      n_monitor = n_monitor,
      k_star = NA_integer_,
      design_id = "E0"
    )
    rows[[ii]] <- data.frame(
      phase = "eprocess_oracle",
      rep_id = rep_id,
      seed = seed,
      design_id = "E0",
      design_family = "null",
      model_type = "quantile",
      tau_grid_id = "Q_005_010",
      instrument_id = "oracle_feature_z",
      m = m,
      T = T,
      n_monitor = n_monitor,
      q = length(tau_grid),
      method_id = ep$id,
      family = "eprocess",
      threshold_type = "eprocess_fixed",
      critical_value = -log(alpha),
      max_stat = max(ep$obj$log_path, na.rm = TRUE),
      stop_time = metrics$stop_time,
      overall_reject = metrics$overall_reject,
      arl0 = metrics$arl0,
      stringsAsFactors = FALSE
    )
    ii <- ii + 1L
  }
  do.call(rbind, rows)
}

run_oracle_eprocess_grid <- function(m_vec, T_vec, B, alpha = 0.05,
                                     seed_base = 20260528,
                                     ncores = 1L,
                                     output_path = NULL,
                                     resume = TRUE,
                                     checkpoint_batch_size = NULL) {
  jobs <- expand.grid(m = m_vec, T = T_vec, rep_id = seq_len(B))
  jobs$run_id <- paste(jobs$m, jobs$T, sep = "_")
  jobs$job_key <- paste(jobs$run_id, jobs$rep_id, sep = "|")

  if (!is.null(output_path) && !isTRUE(resume) && file.exists(output_path)) {
    unlink(output_path)
  }

  existing <- if (!is.null(output_path) && isTRUE(resume)) read_csv_if_exists(output_path) else NULL
  if (!is.null(existing) && nrow(existing) > 0L && all(c("m", "T", "rep_id") %in% names(existing))) {
    done <- unique(paste(paste(existing$m, existing$T, sep = "_"), existing$rep_id, sep = "|"))
    jobs <- jobs[!jobs$job_key %in% done, , drop = FALSE]
    message(sprintf("Resume enabled for %s: found %s completed oracle e-process jobs, %s remaining.",
                    basename(output_path), length(done), nrow(jobs)))
  }

  if (nrow(jobs) == 0L) {
    if (!is.null(existing)) return(existing)
    return(data.frame())
  }

  if (is.null(checkpoint_batch_size)) {
    checkpoint_batch_size <- max(1L, as.integer(max(1L, ncores) * 4L))
  }
  checkpoint_batch_size <- max(1L, as.integer(checkpoint_batch_size))

  project_root <- get_project_root()
  cl <- make_cluster(ncores)
  on.exit(stop_cluster(cl), add = TRUE)

  if (!is.null(cl)) {
    parallel::clusterExport(cl, varlist = c("project_root"), envir = environment())
    parallel::clusterEvalQ(cl, {
      setwd(project_root)
      options(distgc.project_root = project_root)
      source(file.path(project_root, "R", "utils.R"))
      source(file.path(project_root, "R", "weights.R"))
      source(file.path(project_root, "R", "fit_distributional_models.R"))
      source(file.path(project_root, "R", "eprocess.R"))
      source(file.path(project_root, "R", "joe_redesign_runner.R"))
      NULL
    })
  }

  batches <- split(seq_len(nrow(jobs)), ceiling(seq_len(nrow(jobs)) / checkpoint_batch_size))
  collected <- list()
  cc <- 1L

  for (bb in seq_along(batches)) {
    batch_jobs <- jobs[batches[[bb]], , drop = FALSE]
    message(sprintf("  eprocess oracle batch %s/%s: %s jobs", bb, length(batches), nrow(batch_jobs)))

    if (!is.null(cl)) {
      parallel::clusterExport(cl, varlist = c("batch_jobs", "alpha", "seed_base"), envir = environment())
      out_batch <- parallel::parLapply(cl, seq_len(nrow(batch_jobs)), function(ii) {
        run_one_oracle_eprocess_replication(
          rep_id = batch_jobs$rep_id[ii],
          m = batch_jobs$m[ii],
          T = batch_jobs$T[ii],
          alpha = alpha,
          seed_base = seed_base
        )
      })
    } else {
      out_batch <- lapply(seq_len(nrow(batch_jobs)), function(ii) {
        run_one_oracle_eprocess_replication(
          rep_id = batch_jobs$rep_id[ii],
          m = batch_jobs$m[ii],
          T = batch_jobs$T[ii],
          alpha = alpha,
          seed_base = seed_base
        )
      })
    }

    batch_df <- do.call(rbind, out_batch)
    if (!is.null(output_path)) {
      append_csv(batch_df, output_path)
    } else {
      collected[[cc]] <- batch_df
      cc <- cc + 1L
    }
  }

  if (!is.null(output_path)) {
    return(read.csv(output_path, stringsAsFactors = FALSE))
  }
  do.call(rbind, collected)
}

summarise_eprocess_validity <- function(ep_long) {
  df <- ep_long[ep_long$family == "eprocess", , drop = FALSE]
  keys <- c("phase", "design_id", "model_type", "instrument_id", "m", "T",
            "q", "method_id", "threshold_type")
  key <- make_group_key(df, keys)
  groups <- split(seq_len(nrow(df)), key)
  out <- lapply(groups, function(idx) {
    sub <- df[idx, , drop = FALSE]
    p <- mean(sub$overall_reject, na.rm = TRUE)
    B <- sum(!is.na(sub$overall_reject))
    se <- sqrt(p * (1 - p) / max(1, B))
    data.frame(
      phase = sub$phase[1],
      design_id = sub$design_id[1],
      model_type = sub$model_type[1],
      instrument_id = sub$instrument_id[1],
      m = sub$m[1],
      T = sub$T[1],
      q = sub$q[1],
      method_id = sub$method_id[1],
      threshold_type = sub$threshold_type[1],
      false_alarm = p,
      mc_se = se,
      arl0_mean = mean(sub$arl0, na.rm = TRUE),
      B = B,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

write_csv <- function(x, path) {
  dir_create(dirname(path))
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

run_joe_redesign_package <- function(package = c("pilot", "minimum", "full"),
                                     B_cal = NULL,
                                     B_eval = NULL,
                                     B_alt = NULL,
                                     ncores = NULL,
                                     reserve_cores = 3L,
                                     output_dir = NULL,
                                     include_expectile = FALSE,
                                     include_n3 = NULL,
                                     include_n4 = NULL,
                                     tau_grid_override = NULL,
                                     dgp_variant_id = "baseline",
                                     m_override = NULL,
                                     T_override = NULL,
                                     c_override = NULL,
                                     break_override = NULL,
                                     alpha = 0.05,
                                     seed_base = 20260528,
                                     prefer_quantreg = TRUE,
                                     resume = TRUE,
                                     checkpoint_batch_size = NULL) {
  package <- match.arg(package)
  if (is.null(B_cal)) B_cal <- if (package == "pilot") 20L else if (package == "minimum") 2000L else 10000L
  if (is.null(B_eval)) B_eval <- if (package == "pilot") 20L else if (package == "minimum") 1000L else 5000L
  if (is.null(B_alt)) B_alt <- if (package == "pilot") 20L else if (package == "minimum") 1000L else 3000L
  ncores <- reserve_cores_n(reserve = reserve_cores, override = ncores)

  if (is.null(output_dir)) {
    output_dir <- file.path("output", paste0("joe_redesign_", package, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  }
  if (!grepl("^(/|[A-Za-z]:[/\\])", output_dir)) output_dir <- project_path(output_dir)
  dir_create(output_dir)

  message(sprintf("JoE redesign package=%s, B_cal=%s, B_eval=%s, B_alt=%s, ncores=%s (reserve=%s), resume=%s, dgp_variant=%s, tau_grid=%s",
                  package, B_cal, B_eval, B_alt, ncores, reserve_cores, resume,
                  dgp_variant_id, ifelse(is.null(tau_grid_override), "default", tau_grid_override)))
  if (!is.null(m_override) || !is.null(T_override) || !is.null(c_override) || !is.null(break_override)) {
    message(sprintf("Design overrides: m=%s, T=%s, c=%s, break=%s",
                    ifelse(is.null(m_override), "default", paste(m_override, collapse = ",")),
                    ifelse(is.null(T_override), "default", paste(T_override, collapse = ",")),
                    ifelse(is.null(c_override), "default", paste(c_override, collapse = ",")),
                    ifelse(is.null(break_override), "default", paste(break_override, collapse = ","))))
  }
  message("Output directory: ", output_dir)

  cv <- load_critical_values(
    path_base = project_path("critical_values", "critical_values_all.csv"),
    path_weights = project_path("critical_values", "critical_values_all_weights.csv")
  )
  method_grid <- make_main_method_grid(include_eprocess = TRUE)
  design_grid <- build_design_grid(
    package = package,
    include_expectile = include_expectile,
    tau_grid_override = tau_grid_override,
    dgp_variant_id = dgp_variant_id,
    include_n3 = include_n3,
    include_n4 = include_n4,
    m_override = m_override,
    T_override = T_override,
    c_override = c_override,
    break_override = break_override
  )

  write_csv(design_grid, file.path(output_dir, "simulation_design_grid.csv"))
  write_csv(method_grid, file.path(output_dir, "method_grid.csv"))

  null_grid <- design_grid[design_grid$design_family == "null", , drop = FALSE]
  alt_grid <- design_grid[design_grid$design_family == "alternative", , drop = FALSE]

  message("Running null calibration...")
  null_cal_path <- file.path(output_dir, "null_calibration_stats_long.csv")
  null_cal <- run_replication_grid(
    grid = null_grid,
    B = B_cal,
    method_grid = method_grid,
    cv = cv,
    cv_map = NULL,
    alpha = alpha,
    seed_base = seed_base,
    phase = "calibration",
    ncores = ncores,
    prefer_quantreg = prefer_quantreg,
    output_path = null_cal_path,
    resume = resume,
    checkpoint_batch_size = checkpoint_batch_size
  )

  cv_map <- calibrate_null_cv(null_cal, alpha = alpha)
  write_csv(cv_map, file.path(output_dir, "null_empirical_critical_values.csv"))

  message("Running independent null evaluation...")
  null_eval_path <- file.path(output_dir, "null_size_evaluation_long.csv")
  null_eval <- run_replication_grid(
    grid = null_grid,
    B = B_eval,
    method_grid = method_grid,
    cv = cv,
    cv_map = cv_map,
    alpha = alpha,
    seed_base = seed_base + 300000000L,
    phase = "evaluation",
    ncores = ncores,
    prefer_quantreg = prefer_quantreg,
    output_path = null_eval_path,
    resume = resume,
    checkpoint_batch_size = checkpoint_batch_size
  )
  write_csv(summarise_size(null_eval), file.path(output_dir, "null_size_summary.csv"))

  message("Running oracle e-process calibration-null validity check...")
  oracle_ep_path <- file.path(output_dir, "eprocess_oracle_replication_long.csv")
  oracle_ep <- run_oracle_eprocess_grid(
    m_vec = sort(unique(null_grid$m)),
    T_vec = sort(unique(null_grid$T)),
    B = B_eval,
    alpha = alpha,
    seed_base = seed_base,
    ncores = ncores,
    output_path = oracle_ep_path,
    resume = resume,
    checkpoint_batch_size = checkpoint_batch_size
  )
  feasible_ep <- null_eval[null_eval$family == "eprocess", , drop = FALSE]
  eprocess_long <- rbind(
    oracle_ep[, intersect(names(oracle_ep), names(feasible_ep)), drop = FALSE],
    feasible_ep[, intersect(names(oracle_ep), names(feasible_ep)), drop = FALSE]
  )
  write_csv(eprocess_long, file.path(output_dir, "eprocess_replication_long.csv"))
  write_csv(summarise_eprocess_validity(eprocess_long), file.path(output_dir, "eprocess_validity_summary.csv"))

  message("Running alternatives with Brownian and empirical thresholds...")
  alt_long_path <- file.path(output_dir, "alternative_replication_long.csv")
  alt_long <- run_replication_grid(
    grid = alt_grid,
    B = B_alt,
    method_grid = method_grid,
    cv = cv,
    cv_map = cv_map,
    alpha = alpha,
    seed_base = seed_base + 600000000L,
    phase = "alternative",
    ncores = ncores,
    prefer_quantreg = prefer_quantreg,
    output_path = alt_long_path,
    resume = resume,
    checkpoint_batch_size = checkpoint_batch_size
  )
  write_csv(summarise_detection(alt_long), file.path(output_dir, "alternative_power_summary.csv"))

  rsms_diag <- rbind(
    summarise_rsms_offdiag(null_cal),
    summarise_rsms_offdiag(null_eval),
    summarise_rsms_offdiag(alt_long)
  )
  write_csv(
    rbind(
      null_cal[null_cal$standardizer == "RSMS", , drop = FALSE],
      null_eval[null_eval$standardizer == "RSMS", , drop = FALSE],
      alt_long[alt_long$standardizer == "RSMS", , drop = FALSE]
    ),
    file.path(output_dir, "rsms_offdiag_diagnostic.csv")
  )
  write_csv(rsms_diag, file.path(output_dir, "rsms_offdiag_summary.csv"))

  manifest <- data.frame(
    package = package,
    B_cal = B_cal,
    B_eval = B_eval,
    B_alt = B_alt,
    ncores = ncores,
    reserved_cores = reserve_cores,
    resume = resume,
    checkpoint_batch_size = ifelse(is.null(checkpoint_batch_size), NA_integer_, checkpoint_batch_size),
    dgp_variant_id = dgp_variant_id,
    tau_grid_override = ifelse(is.null(tau_grid_override), NA_character_, tau_grid_override),
    m_override = ifelse(is.null(m_override), NA_character_, paste(m_override, collapse = ",")),
    T_override = ifelse(is.null(T_override), NA_character_, paste(T_override, collapse = ",")),
    c_override = ifelse(is.null(c_override), NA_character_, paste(c_override, collapse = ",")),
    break_override = ifelse(is.null(break_override), NA_character_, paste(break_override, collapse = ",")),
    alpha = alpha,
    seed_base = seed_base,
    output_dir = output_dir,
    completed_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
  write_csv(manifest, file.path(output_dir, "run_manifest.csv"))

  invisible(list(
    output_dir = output_dir,
    design_grid = design_grid,
    method_grid = method_grid,
    cv_map = cv_map,
    null_size_summary = summarise_size(null_eval),
    alternative_power_summary = summarise_detection(alt_long),
    rsms_offdiag_summary = rsms_diag
  ))
}

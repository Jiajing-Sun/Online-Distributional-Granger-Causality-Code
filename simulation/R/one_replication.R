# ==============================================================
# one_replication.R -- one Monte Carlo replication
# ==============================================================

source_project("R", "utils.R")
source_project("R", "critical_values.R")
source_project("R", "weights.R")
source_project("R", "monitors.R")
source_project("R", "fit_distributional_models.R")
source_project("R", "dgp_online_distgc.R")
source_project("R", "eprocess.R")

make_stop_metrics <- function(stop_k, n_monitor, k_star, design) {
  reject <- stop_k <= n_monitor
  arl0 <- min(stop_k, n_monitor)

  if (is.na(k_star)) {
    return(list(reject = reject, fa = reject, det = NA, delay = NA, arl0 = arl0))
  }

  if (toupper(design) == "C1") {
    return(list(reject = reject, fa = NA, det = reject, delay = if (reject) stop_k else NA, arl0 = NA))
  }

  fa <- stop_k <= k_star
  det <- (stop_k > k_star) && (stop_k <= n_monitor)
  delay <- if (det) stop_k - k_star else NA_real_

  list(reject = reject, fa = fa, det = det, delay = delay, arl0 = NA)
}

run_one_replication <- function(rep_id,
                                cv,
                                design,
                                m,
                                T,
                                c_val = 0,
                                break_frac = 0.50,
                                burn = 200,
                                model_type = c("quantile", "expectile"),
                                tau_grid,
                                instrument_type = NULL,
                                tau_weight_scheme = c("equal", "lower_beta", "upper_beta", "center_beta"),
                                alpha = 0.05,
                                gamma_vec = c(0, 0.15),
                                cvm_weights = c("U", "Late", "Early", "Mid"),
                                include_eprocess = FALSE,
                                eprocess_tau_grid = NULL,
                                eprocess_feature = c("z", "zminus", "absz"),
                                theta_dict = c(-1, -0.5, 0.5, 1),
                                prefer_quantreg = TRUE,
                                seed_base = 20260310) {
  model_type <- match.arg(model_type)
  tau_weight_scheme <- match.arg(tau_weight_scheme)
  eprocess_feature <- match.arg(eprocess_feature)

  set.seed(as.integer(seed_base + rep_id))

  sim <- simulate_online_distgc(
    design = design,
    m = m,
    T = T,
    c_val = c_val,
    break_frac = break_frac,
    burn = burn
  )

  if (is.null(instrument_type)) instrument_type <- sim$instrument_default

  score_obj <- score_matrix_from_frozen(
    y = sim$y,
    X = sim$X,
    zlag = sim$zlag,
    m = m,
    tau_grid = tau_grid,
    model_type = model_type,
    instrument_type = instrument_type,
    tau_weight_scheme = tau_weight_scheme,
    prefer_quantreg = prefer_quantreg
  )

  psi <- score_obj$psi
  q <- ncol(psi)
  n_monitor <- sim$n_monitor

  out <- data.frame(
    rep_id = rep_id,
    design = design,
    model_type = model_type,
    instrument_type = instrument_type,
    m = m,
    T = T,
    q = q,
    c_val = c_val,
    break_frac = break_frac,
    k_star = if (is.na(sim$k_star)) NA_integer_ else sim$k_star,
    all_fits_converged = all(score_obj$converged),
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # KS-type monitors
  # ----------------------------------------------------------
  for (gamma in gamma_vec) {
    g_tag <- gsub("\\.", "", format(gamma, trim = TRUE))

    cv_ss <- get_critical_value(cv, stat = "SSMS", type = "KS", T = T, q = q, alpha = alpha, gamma = gamma)
    cv_rs <- get_critical_value(cv, stat = "RSMS", type = "KS", T = T, q = q, alpha = alpha, gamma = gamma)
    cv_ha <- get_critical_value(cv, stat = "HAC",  type = "KS", T = T, q = q, alpha = alpha, gamma = gamma)

    ss <- ssms_ks_stop(psi = psi, m = m, T = T, gamma = gamma, crit_val = cv_ss)
    rs <- rsms_ks_stop(psi = psi, m = m, T = T, gamma = gamma, crit_val = cv_rs)
    ha <- hac_ks_stop(psi = psi, m = m, T = T, gamma = gamma, crit_val = cv_ha)

    for (nm in c("SSMS", "RSMS", "HAC")) {
      stop_k <- switch(nm,
                       SSMS = ss$stop_k,
                       RSMS = rs$stop_k,
                       HAC = ha$stop_k)
      met <- make_stop_metrics(stop_k, n_monitor = n_monitor, k_star = sim$k_star, design = design)
      prefix <- paste0(nm, "_KS_g", g_tag)
      out[[paste0(prefix, "_stop")]] <- stop_k
      out[[paste0(prefix, "_rej")]] <- met$reject
      out[[paste0(prefix, "_fa")]] <- met$fa
      out[[paste0(prefix, "_det")]] <- met$det
      out[[paste0(prefix, "_delay")]] <- met$delay
      out[[paste0(prefix, "_arl0")]] <- met$arl0
    }
  }

  # ----------------------------------------------------------
  # CvM-type monitors
  # ----------------------------------------------------------
  for (w in cvm_weights) {
    w_tag <- toupper(w)

    cv_ss <- get_critical_value(cv, stat = "SSMS", type = "CvM", T = T, q = q, alpha = alpha, gamma = 0, weight = w)
    cv_rs <- get_critical_value(cv, stat = "RSMS", type = "CvM", T = T, q = q, alpha = alpha, gamma = 0, weight = w)
    cv_ha <- get_critical_value(cv, stat = "HAC",  type = "CvM", T = T, q = q, alpha = alpha, gamma = 0, weight = w)

    ss <- ssms_cvm_stop(psi = psi, m = m, T = T, weight = w, crit_val = cv_ss)
    rs <- rsms_cvm_stop(psi = psi, m = m, T = T, weight = w, crit_val = cv_rs)
    ha <- hac_cvm_stop(psi = psi, m = m, T = T, weight = w, crit_val = cv_ha)

    for (nm in c("SSMS", "RSMS", "HAC")) {
      stop_k <- switch(nm,
                       SSMS = ss$stop_k,
                       RSMS = rs$stop_k,
                       HAC = ha$stop_k)
      met <- make_stop_metrics(stop_k, n_monitor = n_monitor, k_star = sim$k_star, design = design)
      prefix <- paste0(nm, "_CvM_", w_tag)
      out[[paste0(prefix, "_stop")]] <- stop_k
      out[[paste0(prefix, "_rej")]] <- met$reject
      out[[paste0(prefix, "_fa")]] <- met$fa
      out[[paste0(prefix, "_det")]] <- met$det
      out[[paste0(prefix, "_delay")]] <- met$delay
      out[[paste0(prefix, "_arl0")]] <- met$arl0
    }
  }

  # ----------------------------------------------------------
  # Quantile e-process branch (single-start and restart-bank)
  # ----------------------------------------------------------
  if (isTRUE(include_eprocess) && identical(model_type, "quantile")) {
    if (is.null(eprocess_tau_grid)) {
      eprocess_tau_grid <- tau_grid[tau_grid <= 0.10]
      if (length(eprocess_tau_grid) == 0L) eprocess_tau_grid <- tau_grid[1:min(2L, length(tau_grid))]
    }

    single_ep <- run_quantile_eprocess(
      y = sim$y,
      X = sim$X,
      zlag = sim$zlag,
      m = m,
      n_monitor = n_monitor,
      alpha = alpha,
      tau_grid = eprocess_tau_grid,
      coeff_mat = NULL,
      feature_type = eprocess_feature,
      theta_dict = theta_dict,
      restart_grid = 1L,
      prefer_quantreg = prefer_quantreg,
      return_path = FALSE
    )

    bank_ep <- run_quantile_eprocess(
      y = sim$y,
      X = sim$X,
      zlag = sim$zlag,
      m = m,
      n_monitor = n_monitor,
      alpha = alpha,
      tau_grid = eprocess_tau_grid,
      coeff_mat = NULL,
      feature_type = eprocess_feature,
      theta_dict = theta_dict,
      restart_grid = make_restart_grid(n_monitor, type = "dyadic"),
      prefer_quantreg = prefer_quantreg,
      return_path = FALSE
    )

    for (nm in c("EPROC_MIX", "EPROC_BANK")) {
      stop_k <- switch(nm,
                       EPROC_MIX = single_ep$stop_k,
                       EPROC_BANK = bank_ep$stop_k)
      met <- make_stop_metrics(stop_k, n_monitor = n_monitor, k_star = sim$k_star, design = design)
      out[[paste0(nm, "_stop")]] <- stop_k
      out[[paste0(nm, "_rej")]] <- met$reject
      out[[paste0(nm, "_fa")]] <- met$fa
      out[[paste0(nm, "_det")]] <- met$det
      out[[paste0(nm, "_delay")]] <- met$delay
      out[[paste0(nm, "_arl0")]] <- met$arl0
    }
  }

  out
}

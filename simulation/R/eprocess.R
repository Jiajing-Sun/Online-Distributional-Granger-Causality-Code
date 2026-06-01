# ==============================================================
# eprocess.R -- quantile-hit e-process monitors
# ==============================================================

source_project("R", "utils.R")
source_project("R", "fit_distributional_models.R")

make_restart_grid <- function(n_monitor, type = c("dyadic", "single")) {
  type <- match.arg(type)
  if (type == "single") return(1L)
  out <- 2^(0:floor(log(n_monitor, base = 2)))
  as.integer(unique(out[out <= n_monitor]))
}

make_eprocess_feature <- function(zlag, feature_type = c("z", "zminus", "absz"), m) {
  feature_type <- match.arg(feature_type)
  zlag <- as.numeric(zlag)

  s <- if (feature_type == "z") {
    zlag
  } else if (feature_type == "zminus") {
    pmin(zlag, 0)
  } else {
    abs(zlag)
  }

  mu <- mean(s[1:m])
  sd0 <- stats::sd(s[1:m])
  if (!is.finite(sd0) || sd0 < 1e-8) sd0 <- 1
  (s - mu) / sd0
}

run_quantile_eprocess <- function(y, X, zlag, m, n_monitor,
                                  alpha = 0.05,
                                  tau_grid = c(0.05, 0.10),
                                  coeff_mat = NULL,
                                  feature_type = c("z", "zminus", "absz"),
                                  theta_dict = c(-1, -0.5, 0.5, 1),
                                  theta_weights = NULL,
                                  restart_grid = NULL,
                                  restart_weights = NULL,
                                  tau_weights = NULL,
                                  prefer_quantreg = TRUE,
                                  return_path = TRUE) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  zlag <- as.numeric(zlag)
  tau_grid <- as.numeric(tau_grid)
  feature_type <- match.arg(feature_type)
  theta_dict <- as.numeric(theta_dict)

  if (is.null(coeff_mat)) {
    fit <- fit_frozen_models(
      X_train = X[1:m, , drop = FALSE],
      y_train = y[1:m],
      tau_grid = tau_grid,
      model_type = "quantile",
      prefer_quantreg = prefer_quantreg
    )
    coeff_mat <- fit$coefficients
  }

  if (is.null(theta_weights)) theta_weights <- rep(1, length(theta_dict))
  theta_weights <- normalize_weights(theta_weights)

  if (is.null(restart_grid)) restart_grid <- make_restart_grid(n_monitor = n_monitor, type = "dyadic")
  restart_grid <- as.integer(sort(unique(restart_grid[restart_grid >= 1 & restart_grid <= n_monitor])))
  if (length(restart_grid) == 0L) restart_grid <- 1L

  if (is.null(restart_weights)) restart_weights <- 1 / restart_grid
  restart_weights <- normalize_weights(restart_weights)

  if (is.null(tau_weights)) tau_weights <- rep(1, length(tau_grid))
  tau_weights <- normalize_weights(tau_weights)

  s <- make_eprocess_feature(zlag = zlag, feature_type = feature_type, m = m)
  s_mon <- s[(m + 1):(m + n_monitor)]

  J <- length(tau_grid)
  L <- length(theta_dict)
  Rn <- length(restart_grid)

  # Monitoring hits for each tau
  hits <- matrix(0, nrow = n_monitor, ncol = J)
  for (j in seq_len(J)) {
    tau <- tau_grid[j]
    beta_j <- coeff_mat[, j]
    qhat <- as.vector(X %*% beta_j)
    hits[, j] <- as.numeric(y[(m + 1):(m + n_monitor)] <= qhat[(m + 1):(m + n_monitor)])
  }

  log_theta_w <- log(theta_weights)
  log_restart_w <- log(restart_weights)
  log_tau_w <- log(tau_weights)

  log_components <- array(NA_real_, dim = c(J, L, Rn))
  threshold_log <- -log(alpha)
  stop_k <- n_monitor + 1L
  rejected <- FALSE
  log_path <- if (return_path) rep(NA_real_, n_monitor) else NULL

  for (k in seq_len(n_monitor)) {
    # start new restart components at time k; they include the k-th factor
    new_idx <- which(restart_grid == k)
    if (length(new_idx) > 0L) {
      for (rr in new_idx) {
        log_components[, , rr] <- 0
      }
    }

    s_k <- s_mon[k]
    log_terms_all <- c()

    for (j in seq_len(J)) {
      tau <- tau_grid[j]
      I_k <- hits[k, j]

      logL_vec <- numeric(L)
      for (ell in seq_len(L)) {
        lp <- stats::qlogis(tau) + theta_dict[ell] * s_k
        p_t <- stats::plogis(lp)
        p_t <- min(max(p_t, 1e-8), 1 - 1e-8)
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

    logE_k <- log_sum_exp(log_terms_all)
    if (return_path) log_path[k] <- logE_k

    if (!rejected && is.finite(logE_k) && (logE_k > threshold_log)) {
      stop_k <- k
      rejected <- TRUE
      if (!return_path) break
    }
  }

  list(
    stop_k = stop_k,
    reject = rejected,
    log_threshold = threshold_log,
    log_path = log_path,
    tau_grid = tau_grid,
    theta_dict = theta_dict,
    restart_grid = restart_grid,
    feature_type = feature_type,
    coefficients = coeff_mat
  )
}

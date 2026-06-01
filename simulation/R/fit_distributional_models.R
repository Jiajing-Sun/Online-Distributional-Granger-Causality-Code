# ==============================================================
# fit_distributional_models.R -- frozen quantile / expectile fits and score construction
# ==============================================================

source_project("R", "utils.R")

check_loss <- function(u, tau) {
  u * (tau - as.numeric(u < 0))
}

fit_quantile_optim <- function(X, y, tau,
                               beta_start = NULL,
                               maxit = 4000,
                               reltol = 1e-10) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  p <- ncol(X)

  if (is.null(beta_start)) {
    beta_start <- tryCatch(qr.solve(X, y), error = function(e) rep(stats::median(y), p))
    beta_start <- as.numeric(beta_start)
    if (length(beta_start) != p) beta_start <- c(stats::median(y), rep(0, p - 1L))
  }

  obj <- function(beta) {
    u <- y - as.vector(X %*% beta)
    sum(check_loss(u, tau))
  }

  fit <- stats::optim(
    par = beta_start,
    fn = obj,
    method = "Nelder-Mead",
    control = list(maxit = maxit, reltol = reltol)
  )

  list(coefficients = as.numeric(fit$par),
       converged = (fit$convergence == 0),
       objective = fit$value,
       method = "optim")
}

fit_quantile_model <- function(X, y, tau, prefer_quantreg = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (prefer_quantreg && requireNamespace("quantreg", quietly = TRUE)) {
    fit_try <- tryCatch({
      fit <- quantreg::rq.fit.br(x = X, y = y, tau = tau)
      list(coefficients = as.numeric(fit$coefficients),
           converged = TRUE,
           objective = NA_real_,
           method = "quantreg::rq.fit.br")
    }, error = function(e) NULL)

    if (!is.null(fit_try)) return(fit_try)
  }

  fit_quantile_optim(X = X, y = y, tau = tau)
}

fit_expectile_irls <- function(X, y, tau,
                               beta_start = NULL,
                               max_iter = 500,
                               tol = 1e-8,
                               ridge = 1e-8) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  p <- ncol(X)

  if (is.null(beta_start)) {
    beta_start <- tryCatch(qr.solve(X, y), error = function(e) rep(stats::mean(y), p))
  }
  beta_old <- as.numeric(beta_start)
  if (length(beta_old) != p) beta_old <- c(stats::mean(y), rep(0, p - 1L))

  converged <- FALSE
  n_iter <- 0L

  for (it in seq_len(max_iter)) {
    n_iter <- it
    u <- y - as.vector(X %*% beta_old)
    w <- ifelse(u >= 0, tau, 1 - tau)

    XtWX <- crossprod(X, X * w)
    XtWy <- crossprod(X, y * w)
    beta_new <- as.numeric(safe_qr_solve(XtWX, XtWy, ridge = ridge))

    if (max(abs(beta_new - beta_old)) <= tol * (1 + max(abs(beta_old)))) {
      beta_old <- beta_new
      converged <- TRUE
      break
    }
    beta_old <- beta_new
  }

  list(coefficients = beta_old,
       converged = converged,
       n_iter = n_iter,
       method = "IRLS")
}

fit_frozen_models <- function(X_train, y_train, tau_grid,
                              model_type = c("quantile", "expectile"),
                              prefer_quantreg = TRUE) {
  model_type <- match.arg(model_type)
  X_train <- as.matrix(X_train)
  y_train <- as.numeric(y_train)
  tau_grid <- as.numeric(tau_grid)

  p <- ncol(X_train)
  J <- length(tau_grid)
  coef_mat <- matrix(NA_real_, nrow = p, ncol = J)
  conv <- logical(J)
  methods <- character(J)

  beta_start <- tryCatch(qr.solve(X_train, y_train), error = function(e) rep(stats::mean(y_train), p))
  beta_start <- as.numeric(beta_start)

  for (j in seq_len(J)) {
    tau <- tau_grid[j]
    fit <- if (model_type == "quantile") {
      fit_quantile_model(X_train, y_train, tau = tau, prefer_quantreg = prefer_quantreg)
    } else {
      fit_expectile_irls(X_train, y_train, tau = tau, beta_start = beta_start)
    }
    coef_mat[, j] <- fit$coefficients
    conv[j] <- isTRUE(fit$converged)
    methods[j] <- fit$method
    beta_start <- fit$coefficients
  }

  colnames(coef_mat) <- paste0("tau_", format(tau_grid, trim = TRUE))
  list(coefficients = coef_mat, converged = conv, methods = methods)
}

make_instrument_matrix <- function(zlag,
                                   type = c("z", "asym", "scale", "rich"),
                                   train_idx = NULL) {
  type <- match.arg(type)
  zlag <- as.numeric(zlag)
  if (is.null(train_idx)) train_idx <- seq_along(zlag)
  train_idx <- as.integer(train_idx)

  center_train <- function(x) {
    x - mean(x[train_idx], na.rm = TRUE)
  }

  if (type == "z") {
    H <- matrix(zlag, ncol = 1L)
    colnames(H) <- "z"
    return(H)
  }

  if (type == "asym") {
    H <- cbind(z = zlag, zminus = pmin(zlag, 0))
    return(H)
  }

  if (type == "scale") {
    H <- cbind(z = zlag, z2c = center_train(zlag^2))
    return(H)
  }

  zminus <- pmin(zlag, 0)
  H <- cbind(
    z = zlag,
    zminus = zminus,
    zminus2c = center_train(zminus^2),
    z2c = center_train(zlag^2)
  )
  H
}

make_tau_weights <- function(tau_grid, scheme = c("equal", "lower_beta", "upper_beta", "center_beta")) {
  scheme <- match.arg(scheme)
  tau_grid <- as.numeric(tau_grid)

  if (scheme == "equal") {
    w <- rep(1, length(tau_grid))
  } else if (scheme == "lower_beta") {
    w <- stats::dbeta(tau_grid, shape1 = 2, shape2 = 8)
  } else if (scheme == "upper_beta") {
    w <- stats::dbeta(tau_grid, shape1 = 8, shape2 = 2)
  } else if (scheme == "center_beta") {
    w <- stats::dbeta(tau_grid, shape1 = 5, shape2 = 5)
  }

  normalize_weights(w)
}

score_matrix_from_frozen <- function(y, X, zlag, m, tau_grid,
                                     model_type = c("quantile", "expectile"),
                                     instrument_type = c("z", "asym", "scale", "rich"),
                                     tau_weight_scheme = c("equal", "lower_beta", "upper_beta", "center_beta"),
                                     prefer_quantreg = TRUE) {
  model_type <- match.arg(model_type)
  instrument_type <- match.arg(instrument_type)
  tau_weight_scheme <- match.arg(tau_weight_scheme)

  y <- as.numeric(y)
  X <- as.matrix(X)
  zlag <- as.numeric(zlag)
  tau_grid <- as.numeric(tau_grid)

  n <- length(y)
  if (nrow(X) != n) stop("X and y must have the same number of rows.")
  if (length(zlag) != n) stop("zlag and y must have the same length.")
  if (m >= n) stop("Training size m must be strictly smaller than sample size n.")

  frozen <- fit_frozen_models(
    X_train = X[1:m, , drop = FALSE],
    y_train = y[1:m],
    tau_grid = tau_grid,
    model_type = model_type,
    prefer_quantreg = prefer_quantreg
  )
  coef_mat <- frozen$coefficients

  H <- make_instrument_matrix(zlag = zlag, type = instrument_type, train_idx = seq_len(m))
  dH <- ncol(H)
  J <- length(tau_grid)
  tau_weights <- make_tau_weights(tau_grid = tau_grid, scheme = tau_weight_scheme)

  psi <- matrix(NA_real_, nrow = n, ncol = dH * J)
  col_names <- character(dH * J)
  idx <- 1L

  for (j in seq_len(J)) {
    tau <- tau_grid[j]
    beta_j <- coef_mat[, j]
    u <- y - as.vector(X %*% beta_j)

    score_scalar <- if (model_type == "quantile") {
      tau - as.numeric(u <= 0)
    } else {
      2 * u * abs(tau - as.numeric(u <= 0))
    }

    block <- H * as.numeric(score_scalar)
    block <- tau_weights[j] * block
    psi[, idx:(idx + dH - 1L)] <- block

    block_names <- paste0(colnames(H), "_tau", format(tau, trim = TRUE))
    col_names[idx:(idx + dH - 1L)] <- block_names
    idx <- idx + dH
  }

  colnames(psi) <- col_names

  list(
    psi = psi,
    coefficients = coef_mat,
    converged = frozen$converged,
    fit_methods = frozen$methods,
    H = H,
    tau_weights = tau_weights,
    tau_grid = tau_grid
  )
}

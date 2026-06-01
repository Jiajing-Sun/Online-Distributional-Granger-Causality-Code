# ==============================================================
# dgp_online_distgc.R -- DGPs for the online distributional Granger-causality simulations
# ==============================================================

simulate_online_distgc <- function(design = c("N1", "N2", "N3", "N4",
                                              "A1", "A2", "A3", "A3A", "A3B",
                                              "A4", "A5", "C1", "E0"),
                                   m,
                                   T,
                                   c_val = 0,
                                   break_frac = 0.50,
                                   burn = 200,
                                   phi_y = 0.5,
                                   phi_z = 0.5,
                                   n3_df = 5,
                                   n4_phi_z = 0.8,
                                   garch_omega = 0.10,
                                   garch_alpha = 0.05,
                                   garch_beta = 0.85,
                                   ramp_frac = 0.20,
                                   train_contam_frac = 0.80) {
  design <- match.arg(design)
  if (design == "A3") design <- "A3A"
  n_monitor <- as.integer(round(m * T))
  n_main <- m + n_monitor

  if (!is.finite(n3_df) || n3_df <= 2) stop("n3_df must be finite and larger than 2.")
  if (!is.finite(garch_omega) || !is.finite(garch_alpha) || !is.finite(garch_beta)) {
    stop("GARCH parameters must be finite.")
  }
  if (garch_omega <= 0 || garch_alpha < 0 || garch_beta < 0 || garch_alpha + garch_beta >= 1) {
    stop("GARCH parameters must satisfy omega > 0, alpha >= 0, beta >= 0, alpha + beta < 1.")
  }

  draw_std_t <- function(n = 1L) {
    stats::rt(n, df = n3_df) / sqrt(n3_df / (n3_df - 2))
  }

  is_null <- design %in% c("N1", "N2", "N3", "N4", "E0")
  is_garch <- design %in% c("N2", "A2")
  phi_z_eff <- if (design == "N4") n4_phi_z else phi_z

  if (is_null) {
    k_star <- NA_integer_
    break_index <- Inf
  } else if (design == "C1") {
    k_star <- 0L
    break_index <- as.integer(floor(train_contam_frac * m))
  } else {
    k_star <- as.integer(floor(break_frac * n_monitor))
    break_index <- m + k_star
  }

  L_ramp <- max(1L, as.integer(floor(ramp_frac * n_monitor)))

  # ------------------------------------------------------------
  # Burn-in to obtain approximately stationary initial states
  # ------------------------------------------------------------
  y_prev <- 0
  z_prev <- 0
  sig2_y <- 1
  sig2_z <- 1
  eps_y_prev <- 0
  eps_z_prev <- 0

  for (bb in seq_len(burn + 1L)) {
    if (is_garch) {
      sig2_z <- garch_omega + garch_beta * sig2_z + garch_alpha * eps_z_prev^2
      eps_z <- sqrt(sig2_z) * stats::rnorm(1)
    } else {
      eps_z <- stats::rnorm(1)
    }
      z_curr <- phi_z_eff * z_prev + eps_z

    if (design %in% c("N1", "N4", "E0", "A1", "A5", "C1")) {
      eps_y <- stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design == "N3") {
      eps_y <- draw_std_t(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design == "N2" || design == "A2") {
      sig2_y <- garch_omega + garch_beta * sig2_y + garch_alpha * eps_y_prev^2
      eps_y <- sqrt(sig2_y) * stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design %in% c("A3A", "A3B")) {
      eps_y <- stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design == "A4") {
      eps_y <- stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else {
      stop("Unknown design in burn-in.")
    }

    y_prev <- y_curr
    z_prev <- z_curr
    eps_y_prev <- eps_y
    eps_z_prev <- eps_z
  }

  # ------------------------------------------------------------
  # Main transformed sample: rows i = 1,...,m+n_monitor
  # Each row stores y_i and lagged regressors (y_{i-1}, z_{i-1}).
  # ------------------------------------------------------------
  y <- numeric(n_main)
  ylag <- numeric(n_main)
  zlag <- numeric(n_main)
  d_path <- integer(n_main)
  beta_path <- numeric(n_main)
  sigma_y_path <- numeric(n_main)
  sigma_z_path <- numeric(n_main)

  for (i in seq_len(n_main)) {
    d_i <- if (is.finite(break_index)) as.integer(i > break_index) else 0L

    if (is_garch) {
      sig2_z <- garch_omega + garch_beta * sig2_z + garch_alpha * eps_z_prev^2
      eps_z <- sqrt(sig2_z) * stats::rnorm(1)
      sigma_z_path[i] <- sqrt(sig2_z)
    } else {
      eps_z <- stats::rnorm(1)
      sigma_z_path[i] <- 1
    }
    z_curr <- phi_z_eff * z_prev + eps_z

    beta_i <- 0
    sigma_y_i <- 1

    if (design %in% c("N1", "N4", "E0")) {
      eps_y <- stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
      sigma_y_i <- 1
    } else if (design == "N3") {
      eps_y <- draw_std_t(1)
      y_curr <- phi_y * y_prev + eps_y
      sigma_y_i <- 1
    } else if (design == "N2") {
      sig2_y <- garch_omega + garch_beta * sig2_y + garch_alpha * eps_y_prev^2
      eps_y <- sqrt(sig2_y) * stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
      sigma_y_i <- sqrt(sig2_y)
    } else if (design == "A1") {
      eps_y <- stats::rnorm(1)
      beta_i <- c_val * d_i
      y_curr <- phi_y * y_prev + beta_i * z_prev + eps_y
      sigma_y_i <- 1
    } else if (design == "A2") {
      sig2_y <- garch_omega + garch_beta * sig2_y + garch_alpha * eps_y_prev^2
      eps_y <- sqrt(sig2_y) * stats::rnorm(1)
      beta_i <- c_val * d_i
      y_curr <- phi_y * y_prev + beta_i * z_prev + eps_y
      sigma_y_i <- sqrt(sig2_y)
    } else if (design %in% c("A3A", "A3B")) {
      sigma_y_i <- sqrt(1 + c_val * d_i * (z_prev^2))
      eps_y <- sigma_y_i * stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design == "A4") {
      zminus <- min(z_prev, 0)
      sigma_y_i <- sqrt(1 + c_val * d_i * (zminus^2))
      eps_y <- sigma_y_i * stats::rnorm(1)
      y_curr <- phi_y * y_prev + eps_y
    } else if (design == "A5") {
      ramp_num <- max(i - (m + k_star), 0)
      beta_i <- c_val * min(ramp_num / L_ramp, 1)
      eps_y <- stats::rnorm(1)
      y_curr <- phi_y * y_prev + beta_i * z_prev + eps_y
      sigma_y_i <- 1
    } else if (design == "C1") {
      eps_y <- stats::rnorm(1)
      beta_i <- c_val * d_i
      y_curr <- phi_y * y_prev + beta_i * z_prev + eps_y
      sigma_y_i <- 1
    } else {
      stop("Unknown design.")
    }

    y[i] <- y_curr
    ylag[i] <- y_prev
    zlag[i] <- z_prev
    d_path[i] <- d_i
    beta_path[i] <- beta_i
    sigma_y_path[i] <- sigma_y_i

    y_prev <- y_curr
    z_prev <- z_curr
    eps_y_prev <- eps_y
    eps_z_prev <- eps_z
  }

  list(
    design = design,
    y = y,
    ylag = ylag,
    zlag = zlag,
    X = cbind(Intercept = 1, ylag = ylag),
    n_total = n_main,
    n_monitor = n_monitor,
    m = m,
    T = T,
    c_val = c_val,
    k_star = k_star,
    break_index = break_index,
    d_path = d_path,
    beta_path = beta_path,
    sigma_y_path = sigma_y_path,
    sigma_z_path = sigma_z_path,
    phi_z = phi_z_eff,
    instrument_default = if (design == "A4") {
      "asym"
    } else if (design %in% c("A3B", "N4")) {
      "scale"
    } else {
      "z"
    }
  )
}

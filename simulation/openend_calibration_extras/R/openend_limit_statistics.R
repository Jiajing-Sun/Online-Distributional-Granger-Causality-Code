
# ==============================================================
# openend_limit_statistics.R -- one-replication open-end limits
# ==============================================================

source(file.path("R", "utils_openend.R"))
source(file.path("R", "openend_weights.R"))

build_openend_meta <- function(q_max = 10L,
                               gamma_vec = c(0, 0.15),
                               weight_names = c("U", "Early", "Late", "Mid")) {
  rows <- list()
  idx <- 1L

  for (q in seq_len(q_max)) {
    for (g in gamma_vec) {
      rows[[idx]] <- data.frame(q = q, stat = "SSMS", type = "KS", gamma = g, weight = "", stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(q = q, stat = "RSMS", type = "KS", gamma = g, weight = "", stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(q = q, stat = "HAC",  type = "KS", gamma = g, weight = "", stringsAsFactors = FALSE); idx <- idx + 1L
    }
    for (w in weight_names) {
      rows[[idx]] <- data.frame(q = q, stat = "SSMS", type = "CvM", gamma = 0, weight = w, stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(q = q, stat = "RSMS", type = "CvM", gamma = 0, weight = w, stringsAsFactors = FALSE); idx <- idx + 1L
      rows[[idx]] <- data.frame(q = q, stat = "HAC",  type = "CvM", gamma = 0, weight = w, stringsAsFactors = FALSE); idx <- idx + 1L
    }
  }

  do.call(rbind, rows)
}

simulate_one_openend_rep <- function(q_max = 10L,
                                     n_train_grid = 1500L,
                                     n_open_grid = 2000L,
                                     gamma_vec = c(0, 0.15),
                                     weight_names = c("U", "Early", "Late", "Mid"),
                                     ridge = 1e-10,
                                     range_floor = 1e-8) {
  # --------------------------------------------------------------
  # Training Brownian bridge on [0,1]
  # --------------------------------------------------------------
  dB <- matrix(rnorm(n_train_grid * q_max, sd = 1 / sqrt(n_train_grid)),
               nrow = n_train_grid, ncol = q_max)
  B <- apply(dB, 2, cumsum)
  r <- (1:n_train_grid) / n_train_grid
  Z <- B[n_train_grid, ]
  B0 <- B - outer(r, Z)

  V <- crossprod(B0) / n_train_grid
  ranges <- apply(B0, 2, max) - apply(B0, 2, min)
  ranges[ranges < range_floor] <- range_floor
  inv2 <- 1 / (ranges^2)

  # --------------------------------------------------------------
  # Independent Brownian motion on [0,1] for the open-end monitor
  # --------------------------------------------------------------
  dG <- matrix(rnorm(n_open_grid * q_max, sd = 1 / sqrt(n_open_grid)),
               nrow = n_open_grid, ncol = q_max)
  G <- apply(dG, 2, cumsum)

  # Add x = 0 explicitly to improve trapezoid integration
  x <- c(0, (1:n_open_grid) / n_open_grid)
  G <- rbind(rep(0, q_max), G)

  w_mat <- sapply(weight_names, function(w) openend_weight_x(x, weight = w))
  w_mat <- as.matrix(w_mat)

  values <- numeric(0L)

  for (q in seq_len(q_max)) {
    Gq <- G[, 1:q, drop = FALSE]

    quad_h <- rowSums(Gq^2)
    quad_rs <- rowSums(sweep(Gq^2, 2, inv2[1:q], FUN = "*"))

    Vq_inv <- safe_solve(V[1:q, 1:q, drop = FALSE], ridge = ridge)
    quad_ss <- rowSums((Gq %*% Vq_inv) * Gq)

    for (g in gamma_vec) {
      denom <- x^(2 * g)
      denom[1] <- Inf  # value at x = 0 is immaterial because numerator is 0 there
      values <- c(values,
                  max(quad_ss / denom),
                  max(quad_rs / denom),
                  max(quad_h  / denom))
    }

    for (j in seq_along(weight_names)) {
      wj <- w_mat[, j]
      values <- c(values,
                  trapz_unit_interval(wj * quad_ss),
                  trapz_unit_interval(wj * quad_rs),
                  trapz_unit_interval(wj * quad_h))
    }
  }

  values
}

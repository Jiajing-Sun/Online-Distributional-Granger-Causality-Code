# ==============================================================
# monitors.R -- SSMS / RSMS monitoring (KS and CvM)
# ==============================================================
# Input: score matrix psi_t (n_total x q), with first m rows = training.
# Output: stopping time in monitoring index k \in {1,...,mT} (or mT+1 if no alarm).
#
# IMPORTANT: In the paper, the boundary exponent gamma enters in the DENOMINATOR:
#   M(k) = quad(k) / { m (1+k/m)^2 * (k/(k+m))^{2 gamma} }.
# So gamma>0 inflates the statistic for early k (since k/(k+m) < 1).

source_project("R", "utils.R")
source_project("R", "weights.R")

compute_ssms_sequence <- function(psi, m, T, gamma = 0, ridge = 1e-10) {
  psi <- as.matrix(psi)
  q <- ncol(psi)
  mT <- as.integer(round(m * T))

  # Training mean
  psi_bar <- colMeans(psi[1:m, , drop = FALSE])
  phi <- sweep(psi, 2, psi_bar, FUN = "-")

  phi_train <- phi[1:m, , drop = FALSE]
  phi_mon <- phi[(m + 1):(m + mT), , drop = FALSE]

  # Self-normalizer D_m = (1/m^2) sum_{t=1}^m S_t S_t'
  S_train <- apply(phi_train, 2, cumsum)  # m x q
  Dm <- crossprod(S_train) / (m^2)
  Dm_inv <- safe_solve(Dm, ridge = ridge)

  # Monitoring partial sums
  S_mon <- apply(phi_mon, 2, cumsum)  # mT x q
  quad <- rowSums((S_mon %*% Dm_inv) * S_mon)

  k_vec <- 1:mT
  ratio <- k_vec / (k_vec + m)
  denom <- m * (1 + k_vec / m)^2 * (ratio)^(2 * gamma)

  M <- quad / denom
  return(list(M = M, psi_bar = psi_bar, Dm = Dm))
}

ldl_prewhitener <- function(Sigma0, ridge = 1e-8) {
  Sigma0 <- as.matrix(Sigma0)
  q <- nrow(Sigma0)
  if (q != ncol(Sigma0)) stop("Sigma0 must be square.")

  A <- (Sigma0 + t(Sigma0)) / 2
  scale0 <- mean(diag(A), na.rm = TRUE)
  if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
  ridge_abs <- ridge * scale0
  A <- A + ridge_abs * diag(q)

  L <- diag(q)
  D <- numeric(q)

  for (j in seq_len(q)) {
    prev <- if (j > 1L) seq_len(j - 1L) else integer(0L)
    d_j <- A[j, j]
    if (length(prev) > 0L) {
      d_j <- d_j - sum((L[j, prev]^2) * D[prev])
    }
    if (!is.finite(d_j) || d_j <= ridge_abs) d_j <- ridge_abs
    D[j] <- d_j

    if (j < q) {
      for (i in (j + 1L):q) {
        num <- A[i, j]
        if (length(prev) > 0L) {
          num <- num - sum(L[i, prev] * L[j, prev] * D[prev])
        }
        L[i, j] <- num / D[j]
      }
    }
  }

  L_inv <- forwardsolve(L, diag(q), upper.tri = FALSE)
  H <- diag(1 / sqrt(D), q, q) %*% L_inv
  rownames(H) <- colnames(H) <- colnames(Sigma0)

  list(
    H = H,
    L = L,
    D = D,
    min_ldl_diag = min(D),
    sigma0_condition = tryCatch(kappa(A), error = function(e) NA_real_)
  )
}

matrix_offdiag <- function(A) {
  A <- as.matrix(A)
  diag(A) <- 0
  A
}

compute_delta_off <- function(phi_train,
                              H_whiten = NULL,
                              hac_kernel = "bartlett",
                              hac_bw_rule = "m13",
                              hac_bw_const = 1,
                              hac_bandwidth = NULL,
                              ridge = 1e-8) {
  phi_train <- as.matrix(phi_train)
  Sigma0 <- crossprod(phi_train) / nrow(phi_train)

  ldl <- NULL
  if (is.null(H_whiten)) {
    ldl <- ldl_prewhitener(Sigma0, ridge = ridge)
    H_whiten <- ldl$H
  }

  hac <- estimate_hac_lrv(
    phi_train = phi_train,
    kernel = hac_kernel,
    bandwidth = hac_bandwidth,
    bw_rule = hac_bw_rule,
    bw_const = hac_bw_const,
    ridge = ridge
  )

  transformed <- H_whiten %*% hac$Sigma_hat %*% t(H_whiten)
  denom <- norm(transformed, type = "F")
  if (!is.finite(denom) || denom <= .Machine$double.eps) denom <- .Machine$double.eps
  delta_off <- norm(matrix_offdiag(transformed), type = "F") / denom

  Dv <- sqrt(pmax(diag(transformed), .Machine$double.eps))
  corr <- transformed / outer(Dv, Dv)
  maxcorr <- max(abs(matrix_offdiag(corr)), na.rm = TRUE)

  list(
    delta_off = as.numeric(delta_off),
    lrv_offdiag_maxcorr = as.numeric(maxcorr),
    transformed_lrv = transformed,
    hac_bandwidth = hac$h,
    sigma0_condition = if (!is.null(ldl)) ldl$sigma0_condition else NA_real_,
    min_ldl_diag = if (!is.null(ldl)) ldl$min_ldl_diag else NA_real_
  )
}

ssms_ks_stop <- function(psi, m, T, gamma, crit_val, ridge = 1e-10) {
  seq <- compute_ssms_sequence(psi, m, T, gamma = gamma, ridge = ridge)
  M <- seq$M
  mT <- length(M)
  k_stop <- which(M > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(M), M = M,
              delta_off = seq$delta_off,
              lrv_offdiag_maxcorr = seq$lrv_offdiag_maxcorr,
              sigma0_condition = seq$sigma0_condition,
              min_ldl_diag = seq$min_ldl_diag,
              rsms_hac_bandwidth = seq$rsms_hac_bandwidth))
}

ssms_cvm_stop <- function(psi, m, T, weight, crit_val, ridge = 1e-10) {
  # CvM uses gamma = 0 in the paper
  seq <- compute_ssms_sequence(psi, m, T, gamma = 0, ridge = ridge)
  M <- seq$M
  mT <- length(M)
  k_vec <- 1:mT
  w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = weight)
  I <- cumsum(w_vec * M) / m
  k_stop <- which(I > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(I), I = I, M = M,
              delta_off = seq$delta_off,
              lrv_offdiag_maxcorr = seq$lrv_offdiag_maxcorr,
              sigma0_condition = seq$sigma0_condition,
              min_ldl_diag = seq$min_ldl_diag,
              rsms_hac_bandwidth = seq$rsms_hac_bandwidth))
}

compute_rsms_sequence <- function(psi, m, T, gamma = 0,
                                  ridge = 1e-8,
                                  range_floor = 1e-8) {
  psi <- as.matrix(psi)
  q <- ncol(psi)
  mT <- as.integer(round(m * T))

  # Training mean
  psi_bar <- colMeans(psi[1:m, , drop = FALSE])
  phi <- sweep(psi, 2, psi_bar, FUN = "-")

  phi_train <- phi[1:m, , drop = FALSE]
  phi_mon <- phi[(m + 1):(m + mT), , drop = FALSE]

  # Lag-0 LDL prewhitening (training). This is a coordinate transformation,
  # not a long-run covariance estimator.
  Sigma0 <- crossprod(phi_train) / m  # q x q
  ldl <- ldl_prewhitener(Sigma0, ridge = ridge)
  phi_w <- phi %*% t(ldl$H)

  phi_w_train <- phi_w[1:m, , drop = FALSE]
  phi_w_mon <- phi_w[(m + 1):(m + mT), , drop = FALSE]

  # Range normalizer on training Brownian bridge
  S_train <- apply(phi_w_train, 2, cumsum)  # m x q
  # Brownian bridge adjustment: B(t) = S(t) - (t/m) S(m)
  t_vec <- (1:m) / m
  S_end <- matrix(rep(S_train[m, ], each = m), nrow = m, ncol = q)
  B_train <- S_train - S_end * t_vec

  ranges <- apply(B_train, 2, function(x) max(x) - min(x))
  ranges[ranges < range_floor] <- range_floor

  inv2 <- m / (ranges^2)  # diag entries of R_m^{-2}

  # Monitoring partial sums (whitened)
  S_mon <- apply(phi_w_mon, 2, cumsum)  # mT x q
  quad <- rowSums(sweep(S_mon^2, 2, inv2, FUN = "*"))

  k_vec <- 1:mT
  ratio <- k_vec / (k_vec + m)
  denom <- m * (1 + k_vec / m)^2 * (ratio)^(2 * gamma)

  M <- quad / denom

  diag <- compute_delta_off(
    phi_train = phi_train,
    H_whiten = ldl$H,
    hac_kernel = "bartlett",
    hac_bw_rule = "m13",
    hac_bw_const = 1,
    ridge = ridge
  )

  return(list(M = M, psi_bar = psi_bar, Sigma0 = Sigma0,
              inv2 = inv2, H_whiten = ldl$H,
              min_ldl_diag = ldl$min_ldl_diag,
              sigma0_condition = ldl$sigma0_condition,
              delta_off = diag$delta_off,
              lrv_offdiag_maxcorr = diag$lrv_offdiag_maxcorr,
              rsms_hac_bandwidth = diag$hac_bandwidth))
}

rsms_ks_stop <- function(psi, m, T, gamma, crit_val,
                         ridge = 1e-8, range_floor = 1e-8) {
  seq <- compute_rsms_sequence(psi, m, T, gamma = gamma,
                               ridge = ridge, range_floor = range_floor)
  M <- seq$M
  mT <- length(M)
  k_stop <- which(M > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(M), M = M))
}

rsms_cvm_stop <- function(psi, m, T, weight, crit_val, gamma = 0,
                          ridge = 1e-8, range_floor = 1e-8) {
  seq <- compute_rsms_sequence(psi, m, T, gamma = gamma,
                               ridge = ridge, range_floor = range_floor)
  M <- seq$M
  mT <- length(M)
  k_vec <- 1:mT
  w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = weight)
  I <- cumsum(w_vec * M) / m
  k_stop <- which(I > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(I), I = I, M = M))
}

# ==============================================================
# HAC baseline monitor
#   - Standardize by a training-window HAC/LRV estimator
#   - See Appendix \ref{sec:HACbaseline} in the paper
#
# Notation (matches the paper):
#   phi_t = psi_t - \bar\psi_m
#   S_m(k) = \sum_{t=m+1}^{m+k} phi_t
#   \widehat\Sigma_{m,h} = \widehat\Gamma_0 + \sum_{\ell=1}^h K(\ell/h) (\widehat\Gamma_\ell + \widehat\Gamma_\ell')
# ==============================================================

hac_kernel_weight <- function(x, kernel = c("bartlett")) {
  # x in [0,1]
  kernel <- match.arg(tolower(kernel), c("bartlett", "parzen", "quadratic_spectral", "qs"))
  u <- abs(x)
  if (kernel == "bartlett") {
    # Bartlett / Newey-West: K(u) = 1 - |u| on |u|<=1
    return(pmax(0, 1 - u))
  }
  if (kernel == "parzen") {
    out <- numeric(length(u))
    idx1 <- u <= 0.5
    idx2 <- u > 0.5 & u <= 1
    out[idx1] <- 1 - 6 * u[idx1]^2 + 6 * u[idx1]^3
    out[idx2] <- 2 * (1 - u[idx2])^3
    return(out)
  }
  if (kernel %in% c("quadratic_spectral", "qs")) {
    out <- numeric(length(u))
    near0 <- u < 1e-12
    out[near0] <- 1
    v <- u[!near0]
    x0 <- 6 * pi * v / 5
    out[!near0] <- 25 / (12 * pi^2 * v^2) * (sin(x0) / x0 - cos(x0))
    return(out)
  }
  stop("Unknown HAC kernel: ", kernel)
}

pick_hac_bandwidth <- function(m, rule = c("m13", "nw"), c = 1) {
  # Simple, portable bandwidth rules.
  # - m13: floor(c * m^{1/3})
  # - nw : Newey-West style: floor(c * 4 * (m/100)^{2/9})
  rule <- match.arg(rule)
  if (rule == "m13") {
    h <- floor(c * (m^(1/3)))
  } else if (rule == "nw") {
    h <- floor(c * 4 * (m / 100)^(2/9))
  }
  h <- as.integer(max(0, h))
  return(h)
}

estimate_hac_lrv <- function(phi_train,
                             kernel = "bartlett",
                             bandwidth = NULL,
                             bw_rule = "m13",
                             bw_const = 1,
                             ridge = 1e-10) {
  # phi_train: m x q matrix
  phi_train <- as.matrix(phi_train)
  m <- nrow(phi_train)
  q <- ncol(phi_train)

  if (is.null(bandwidth)) {
    h <- pick_hac_bandwidth(m = m, rule = bw_rule, c = bw_const)
  } else {
    h <- as.integer(max(0, round(bandwidth)))
  }
  if (h >= m) h <- m - 1L

  # Gamma_0
  Sigma_hat <- crossprod(phi_train) / m

  if (h >= 1L) {
    for (ell in 1:h) {
      w_ell <- hac_kernel_weight(ell / h, kernel = kernel)
      if (w_ell == 0) next

      A <- phi_train[(ell + 1):m, , drop = FALSE]
      B <- phi_train[1:(m - ell), , drop = FALSE]
      Gamma_ell <- crossprod(A, B) / m
      Sigma_hat <- Sigma_hat + w_ell * (Gamma_ell + t(Gamma_ell))
    }
  }

  Sigma_inv <- safe_solve(Sigma_hat, ridge = ridge)
  return(list(Sigma_hat = Sigma_hat, Sigma_inv = Sigma_inv, h = h))
}

compute_hac_sequence <- function(psi, m, T, gamma = 0,
                                 kernel = "bartlett",
                                 bandwidth = NULL,
                                 bw_rule = "m13",
                                 bw_const = 1,
                                 ridge = 1e-10) {
  psi <- as.matrix(psi)
  q <- ncol(psi)
  mT <- as.integer(round(m * T))

  psi_bar <- colMeans(psi[1:m, , drop = FALSE])
  phi <- sweep(psi, 2, psi_bar, FUN = "-")

  phi_train <- phi[1:m, , drop = FALSE]
  phi_mon <- phi[(m + 1):(m + mT), , drop = FALSE]

  hac <- estimate_hac_lrv(phi_train,
                          kernel = kernel,
                          bandwidth = bandwidth,
                          bw_rule = bw_rule,
                          bw_const = bw_const,
                          ridge = ridge)
  Sigma_inv <- hac$Sigma_inv

  S_mon <- apply(phi_mon, 2, cumsum)  # mT x q
  quad <- rowSums((S_mon %*% Sigma_inv) * S_mon)

  k_vec <- 1:mT
  ratio <- k_vec / (k_vec + m)
  denom <- m * (1 + k_vec / m)^2 * (ratio)^(2 * gamma)

  M <- quad / denom
  return(list(M = M, psi_bar = psi_bar,
              Sigma_hat = hac$Sigma_hat, h = hac$h))
}

hac_ks_stop <- function(psi, m, T, gamma, crit_val,
                        kernel = "bartlett",
                        bandwidth = NULL,
                        bw_rule = "m13",
                        bw_const = 1,
                        ridge = 1e-10) {
  seq <- compute_hac_sequence(psi, m, T, gamma = gamma,
                              kernel = kernel,
                              bandwidth = bandwidth,
                              bw_rule = bw_rule,
                              bw_const = bw_const,
                              ridge = ridge)
  M <- seq$M
  mT <- length(M)
  k_stop <- which(M > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(M), M = M,
              Sigma_hat = seq$Sigma_hat, h = seq$h))
}

hac_cvm_stop <- function(psi, m, T, weight, crit_val,
                         kernel = "bartlett",
                         bandwidth = NULL,
                         bw_rule = "m13",
                         bw_const = 1,
                         ridge = 1e-10) {
  # CvM uses gamma = 0
  seq <- compute_hac_sequence(psi, m, T, gamma = 0,
                              kernel = kernel,
                              bandwidth = bandwidth,
                              bw_rule = bw_rule,
                              bw_const = bw_const,
                              ridge = ridge)
  M <- seq$M
  mT <- length(M)
  k_vec <- 1:mT
  w_vec <- make_cvm_weight(k_vec, m = m, T = T, weight = weight)
  I <- cumsum(w_vec * M) / m
  k_stop <- which(I > crit_val)[1]
  if (is.na(k_stop)) k_stop <- mT + 1L
  return(list(stop_k = k_stop, reject = (k_stop <= mT),
              max_stat = max(I), I = I, M = M,
              Sigma_hat = seq$Sigma_hat, h = seq$h))
}

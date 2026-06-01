
# ==============================================================
# utils_openend.R -- small helpers for open-end critical values
# ==============================================================

safe_solve <- function(A, ridge = 1e-10) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A))
  q <- nrow(A)
  A2 <- A + ridge * diag(q)
  out <- tryCatch({
    chol2inv(chol(A2))
  }, error = function(e) {
    solve(A2)
  })
  out
}

make_cluster <- function(ncores) {
  if (ncores <= 1L) return(NULL)
  parallel::makeCluster(ncores, type = "PSOCK")
}

stop_cluster <- function(cl) {
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }
}

trapz_unit_interval <- function(y) {
  # Simple trapezoid rule on an equally spaced grid over [0,1].
  n <- length(y)
  if (n < 2L) return(0)
  dx <- 1 / (n - 1)
  dx * (0.5 * y[1] + sum(y[2:(n - 1)]) + 0.5 * y[n])
}

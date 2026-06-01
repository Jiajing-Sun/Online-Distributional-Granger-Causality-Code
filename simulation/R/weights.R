# ==============================================================
# weights.R -- CvM weights w(s) on [0,T]
# ==============================================================
# We use normalized time tau = s/T in [0,1].
# All weights are scaled to satisfy int_0^T w(s) ds = T (unit average weight).

make_cvm_weight <- function(k_vec, m, T, weight = "U") {
  # k_vec: monitoring indices 1,2,...,mT (integers)
  # m: training size
  # T: horizon (e.g., 1,2,5,10)
  # weight: "U", "Late", "Early", "Mid" (case-insensitive)

  w <- toupper(as.character(weight))
  tau <- (k_vec / m) / T  # tau = s/T = (k/m)/T = k/(mT)

  if (w %in% c("U", "UNIFORM", "CONST", "CONSTANT", "ONE", "WU", "W_U", "1")) {
    return(rep(1, length(k_vec)))
  }
  if (w %in% c("LATE", "W_LATE")) {
    return(2 * tau)
  }
  if (w %in% c("EARLY", "W_EARLY")) {
    return(2 * (1 - tau))
  }
  if (w %in% c("MID", "W_MID")) {
    return(6 * tau * (1 - tau))
  }

  stop("Unknown CvM weight. Use one of: 'U', 'Late', 'Early', 'Mid'.")
}

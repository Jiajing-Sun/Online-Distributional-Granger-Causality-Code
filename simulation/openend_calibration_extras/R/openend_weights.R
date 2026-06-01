
# ==============================================================
# openend_weights.R -- admissible open-end CvM weights
# ==============================================================

openend_weight_x <- function(x, weight = "U") {
  w <- toupper(as.character(weight)[1])

  if (w %in% c("U", "UNIFORM", "WU", "W_U", "1")) {
    return(rep(1, length(x)))
  }
  if (w %in% c("EARLY", "W_EARLY")) {
    return(2 * (1 - x))
  }
  if (w %in% c("LATE", "W_LATE")) {
    return(2 * x)
  }
  if (w %in% c("MID", "W_MID")) {
    return(6 * x * (1 - x))
  }

  stop("Unknown open-end weight. Use one of: U, Early, Late, Mid.")
}

openend_weight_s <- function(s, weight = "U") {
  w <- toupper(as.character(weight)[1])

  if (w %in% c("U", "UNIFORM", "WU", "W_U", "1")) {
    return((1 + s)^(-2))
  }
  if (w %in% c("EARLY", "W_EARLY")) {
    return(2 * (1 + s)^(-3))
  }
  if (w %in% c("LATE", "W_LATE")) {
    return(2 * s * (1 + s)^(-3))
  }
  if (w %in% c("MID", "W_MID")) {
    return(6 * s * (1 + s)^(-4))
  }

  stop("Unknown open-end weight. Use one of: U, Early, Late, Mid.")
}

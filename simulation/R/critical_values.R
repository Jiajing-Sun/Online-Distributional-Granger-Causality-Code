# ==============================================================
# critical_values.R -- load, validate, and look up critical values
# ==============================================================

source_project("R", "utils.R")

normalize_cv_table <- function(df, table_name = "base") {
  if (is.null(df)) return(NULL)

  required <- c("stat", "type", "T", "q", "alpha", "critical_value")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0L) {
    stop(sprintf("Critical-values %s table is missing required columns: %s",
                 table_name, paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  out <- df
  out$stat <- toupper(trimws(as.character(out$stat)))

  type_upper <- toupper(trimws(as.character(out$type)))
  out$type <- ifelse(type_upper == "CVM", "CvM",
                     ifelse(type_upper == "KS", "KS", trimws(as.character(out$type))))

  if (!"gamma" %in% names(out)) {
    out$gamma <- NA_real_
  } else {
    out$gamma <- as.numeric(out$gamma)
  }

  if ("weight_name" %in% names(out)) {
    out$weight_name <- trimws(as.character(out$weight_name))
  }

  out$T <- as.numeric(out$T)
  out$q <- as.integer(out$q)
  out$alpha <- as.numeric(out$alpha)
  out$critical_value <- as.numeric(out$critical_value)

  out
}

make_cv_key <- function(df, cols) {
  pieces <- lapply(cols, function(cc) {
    x <- df[[cc]]
    if (is.numeric(x)) {
      x <- ifelse(is.na(x), "NA",
                  format(x, digits = 16, trim = TRUE, scientific = FALSE))
    } else {
      x <- ifelse(is.na(x), "NA", trimws(as.character(x)))
    }
    as.character(x)
  })
  do.call(paste, c(pieces, sep = "|"))
}

validate_unique_keys <- function(df, cols, table_name) {
  keys <- make_cv_key(df, cols)
  dup <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  if (any(dup)) {
    bad <- unique(keys[dup])
    stop(sprintf("Duplicate keys found in %s critical-values table. First duplicate key: %s",
                 table_name, bad[1]),
         call. = FALSE)
  }
  invisible(TRUE)
}

validate_critical_value_tables <- function(base, weights = NULL) {
  if (is.null(base)) stop("Base critical-values table is NULL.", call. = FALSE)

  if (any(!is.finite(base$critical_value))) {
    stop("Base critical-values table contains non-finite critical values.", call. = FALSE)
  }
  validate_unique_keys(base, c("stat", "type", "T", "gamma", "q", "alpha"), "base")

  if (!is.null(weights)) {
    if (any(toupper(weights$type) != "CVM")) {
      stop("Weights critical-values table must contain CvM rows only.", call. = FALSE)
    }
    if (!"weight_name" %in% names(weights)) {
      stop("Weights critical-values table is missing the 'weight_name' column.", call. = FALSE)
    }
    if (any(!is.finite(weights$critical_value))) {
      stop("Weights critical-values table contains non-finite critical values.", call. = FALSE)
    }
    validate_unique_keys(weights, c("stat", "type", "weight_name", "T", "gamma", "q", "alpha"), "weights")
  }

  invisible(TRUE)
}

load_critical_values <- function(path_base = file.path("critical_values", "critical_values_all.csv"),
                                 path_weights = file.path("critical_values", "critical_values_all_weights.csv")) {
  if (!file.exists(path_base)) {
    stop("Cannot find base critical values CSV at: ", path_base, call. = FALSE)
  }

  base <- read.csv(path_base, stringsAsFactors = FALSE)
  base <- normalize_cv_table(base, table_name = "base")

  weights <- NULL
  if (!is.null(path_weights) && file.exists(path_weights)) {
    weights <- read.csv(path_weights, stringsAsFactors = FALSE)
    weights <- normalize_cv_table(weights, table_name = "weights")
    weights <- weights[toupper(weights$type) == "CVM", , drop = FALSE]
  }

  validate_critical_value_tables(base = base, weights = weights)
  list(base = base, weights = weights)
}

match_numeric_scalar <- function(x, target, tol = 1e-10) {
  (!is.na(x)) & (!is.na(target)) & (abs(x - target) < tol)
}

refine_by_gamma <- function(sub, gamma, context, tol = 1e-10) {
  if (!"gamma" %in% names(sub) || nrow(sub) == 0L) return(sub)

  has_non_na_gamma <- any(!is.na(sub$gamma))
  if (!has_non_na_gamma) {
    return(sub)
  }

  if (is.na(gamma)) {
    stop(sprintf("A gamma value is required for %s because the critical-values table indexes this case by gamma.",
                 context),
         call. = FALSE)
  }

  sub <- sub[match_numeric_scalar(sub$gamma, gamma, tol = tol), , drop = FALSE]
  sub
}

finalize_cv_match <- function(sub, context) {
  if (nrow(sub) != 1L) {
    stop(sprintf("Critical value not uniquely found for %s. Matches=%s. Please check the CSV tables.",
                 context, nrow(sub)),
         call. = FALSE)
  }

  cv_out <- as.numeric(sub$critical_value[1])
  if (length(cv_out) != 1L || !is.finite(cv_out)) {
    stop(sprintf("Non-finite critical value encountered for %s. Please check the CSV tables.",
                 context),
         call. = FALSE)
  }
  cv_out
}

get_critical_value <- function(cv,
                               stat = c("SSMS", "RSMS", "HAC"),
                               type = c("KS", "CvM"),
                               T,
                               q,
                               alpha = 0.05,
                               gamma = NA,
                               weight = "U") {
  stat <- toupper(as.character(stat)[1])
  type0 <- as.character(type)[1]
  type_upper <- toupper(type0)

  if (type_upper == "KS") {
    df <- cv$base
    if (is.null(df)) stop("Base critical values table is NULL.", call. = FALSE)

    sub <- df[
      (df$stat == stat) &
      (toupper(df$type) == "KS") &
      (df$T == T) &
      (df$q == q) &
      (abs(df$alpha - alpha) < 1e-12),
      , drop = FALSE
    ]

    context <- sprintf("(%s,%s,T=%s,q=%s,alpha=%s,gamma=%s)",
                       stat, type0, T, q, alpha, gamma)
    sub <- refine_by_gamma(sub, gamma = gamma, context = context)
    return(finalize_cv_match(sub, context = context))
  }

  if (type_upper == "CVM") {
    w_upper <- toupper(as.character(weight)[1])
    use_base <- w_upper %in% c("U", "UNIFORM", "CONST", "CONSTANT", "ONE", "WU", "W_U", "1")

    if (use_base) {
      df <- cv$base
      if (is.null(df)) stop("Base critical values table is NULL.", call. = FALSE)

      sub <- df[
        (df$stat == stat) &
        (toupper(df$type) == "CVM") &
        (df$T == T) &
        (df$q == q) &
        (abs(df$alpha - alpha) < 1e-12),
        , drop = FALSE
      ]

      context <- sprintf("(%s,CvM[U],T=%s,q=%s,alpha=%s,gamma=%s)",
                         stat, T, q, alpha, gamma)
      sub <- refine_by_gamma(sub, gamma = gamma, context = context)
      return(finalize_cv_match(sub, context = context))
    }

    dfw <- cv$weights
    if (is.null(dfw)) {
      stop("Requested weighted CvM critical value but the weights table was not loaded.",
           call. = FALSE)
    }

    sub <- dfw[
      (dfw$stat == stat) &
      (toupper(dfw$type) == "CVM") &
      (dfw$T == T) &
      (dfw$q == q) &
      (abs(dfw$alpha - alpha) < 1e-12) &
      (toupper(dfw$weight_name) == w_upper),
      , drop = FALSE
    ]

    context <- sprintf("(%s,CvM[%s],T=%s,q=%s,alpha=%s,gamma=%s)",
                       stat, w_upper, T, q, alpha, gamma)
    sub <- refine_by_gamma(sub, gamma = gamma, context = context)
    return(finalize_cv_match(sub, context = context))
  }

  stop("Unknown type: ", type0, ". Use 'KS' or 'CvM'.", call. = FALSE)
}

# ==============================================================
# method_catalog.R -- method registries and pretty names
# ==============================================================

get_method_catalog <- function(which = c("paper", "all"), model_type = c("quantile", "expectile")) {
  which <- match.arg(which)
  model_type <- match.arg(model_type)

  ks_methods <- c(
    "SSMS_KS_g0", "SSMS_KS_g015",
    "RSMS_KS_g0", "RSMS_KS_g015",
    "HAC_KS_g0",  "HAC_KS_g015"
  )
  cvm_methods <- c(
    "SSMS_CvM_U", "SSMS_CvM_EARLY", "SSMS_CvM_MID", "SSMS_CvM_LATE",
    "RSMS_CvM_U", "RSMS_CvM_EARLY", "RSMS_CvM_MID", "RSMS_CvM_LATE",
    "HAC_CvM_U",  "HAC_CvM_EARLY",  "HAC_CvM_MID",  "HAC_CvM_LATE"
  )
  e_methods <- c("EPROC_MIX", "EPROC_BANK")

  if (which == "paper") {
    methods <- c("SSMS_KS_g0", "RSMS_KS_g0", "HAC_KS_g0", "RSMS_CvM_LATE", "HAC_CvM_LATE")
    if (identical(model_type, "quantile")) methods <- c(methods, e_methods)
    return(methods)
  }

  methods <- c(ks_methods, cvm_methods)
  if (identical(model_type, "quantile")) methods <- c(methods, e_methods)
  methods
}

pretty_method_name <- function(x) {
  out <- x
  out <- gsub("_CvM_", "-CvM-", out, fixed = TRUE)
  out <- gsub("_KS_", "-KS-", out, fixed = TRUE)
  out <- gsub("SSMS", "SSMS", out, fixed = TRUE)
  out <- gsub("RSMS", "RSMS", out, fixed = TRUE)
  out <- gsub("HAC", "HAC", out, fixed = TRUE)
  out <- gsub("EPROC_MIX", "E-process mix", out, fixed = TRUE)
  out <- gsub("EPROC_BANK", "E-process multi-start", out, fixed = TRUE)
  out <- gsub("-g0$", " (gamma=0)", out)
  out <- gsub("-g015$", " (gamma=0.15)", out)
  out <- gsub("-U$", " (U)", out)
  out <- gsub("-EARLY$", " (Early)", out)
  out <- gsub("-MID$", " (Mid)", out)
  out <- gsub("-LATE$", " (Late)", out)
  out
}

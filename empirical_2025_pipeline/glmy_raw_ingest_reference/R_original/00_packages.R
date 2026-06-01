# ============================================================
# 00_packages.R
# Utility to install/load required packages.
# ============================================================

ensure_packages <- function(pkgs, repos = getOption("repos")) {
  stopifnot(is.character(pkgs))
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing package: ", p)
      install.packages(p, repos = repos)
    }
  }
  invisible(TRUE)
}

load_packages <- function(pkgs) {
  ensure_packages(pkgs)
  for (p in pkgs) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(TRUE)
}

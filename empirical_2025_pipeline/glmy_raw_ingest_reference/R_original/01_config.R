# ============================================================
# 01_config.R
# Project configuration + path helpers.
#
# This project supports TWO common folder layouts.
#
# Layout A (raw folders at code root):
#   ./bitfinex/
#   ./bitstamp/
#   ./btc-options/
#   ./Coingecko/
#   ./BTC TX Network/   (or Btctxnetwork/)
#   ./data/             (created for outputs)
#
# Layout B (raw folders inside ./data):
#   ./data/bitfinex/
#   ./data/bitstamp/
#   ./data/btc-options/
#   ./data/Coingecko/
#   ./data/BTC TX Network/
#   ./data/{processed,graphs,glmy,cpd,figures}/ (created for outputs)
#
# In both cases, outputs are always written under ONE "data" folder
# (never "data/data").
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Detect BTC TX folder variants (some systems dislike spaces)
.detect_btc_tx_dir <- function(raw_root) {
  cand <- c(
    file.path(raw_root, "BTC TX Network"),
    file.path(raw_root, "Btctxnetwork"),
    file.path(raw_root, "BTC_TX_Network"),
    file.path(raw_root, "btc_tx_network"),
    file.path(raw_root, "BTC-TX-Network")
  )
  cand[dir.exists(cand)][1] %||% file.path(raw_root, "BTC TX Network")
}

# Heuristic: choose where the raw folders live (either code_root or code_root/data)
.detect_raw_root <- function(code_root) {
  cand <- c(code_root, file.path(code_root, "data"))
  expected <- c("bitfinex", "bitstamp", "btc-options", "Coingecko", "BTC TX Network", "Btctxnetwork")
  scores <- vapply(cand, function(root) sum(dir.exists(file.path(root, expected))), numeric(1))
  cand[which.max(scores)]
}

# Main config
get_config <- function(project_root = getwd()) {
  code_root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)

  raw_root <- .detect_raw_root(code_root)

  # Where to write outputs:
  # - If the code_root itself is named "data", treat that as the output root.
  # - Else if raw lives in code_root/data, write outputs into the same folder (avoid data/data).
  # - Else write outputs into code_root/data.
  if (basename(code_root) == "data") {
    out_root <- code_root
  } else if (identical(raw_root, file.path(code_root, "data"))) {
    out_root <- raw_root
  } else {
    out_root <- file.path(code_root, "data")
  }

  cfg <- list(
    project_root = code_root,
    raw_root = raw_root,

    # output root (everything we create)
    out_root = out_root,

    # raw data roots
    raw = list(
      bitfinex  = file.path(raw_root, "bitfinex"),
      bitstamp  = file.path(raw_root, "bitstamp"),
      options   = file.path(raw_root, "btc-options"),
      coingecko = file.path(raw_root, "Coingecko"),
      btc_tx    = .detect_btc_tx_dir(raw_root)
    ),

    # subfolders under out_root
    out = list(
      processed = file.path(out_root, "processed"),
      graphs    = file.path(out_root, "graphs"),
      glmy      = file.path(out_root, "glmy"),
      cpd       = file.path(out_root, "cpd"),
      figures   = file.path(out_root, "figures")
    )
  )

  cfg
}

ensure_output_dirs <- function(cfg) {
  stopifnot(is.list(cfg), !is.null(cfg$out_root))
  dirs <- c(cfg$out_root, unlist(cfg$out, use.names = FALSE))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

# Helper to create module-specific subdirs
module_dirs <- function(cfg, module_name) {
  stopifnot(is.character(module_name), length(module_name) == 1)
  out <- list(
    processed = file.path(cfg$out$processed, module_name),
    graphs    = file.path(cfg$out$graphs, module_name),
    glmy      = file.path(cfg$out$glmy, module_name),
    cpd       = file.path(cfg$out$cpd, module_name),
    figures   = file.path(cfg$out$figures, module_name)
  )
  for (d in out) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  out
}

log_msg <- function(..., .prefix = "[GLMY-CRYPTO]") {
  message(.prefix, " ", sprintf(...))
}

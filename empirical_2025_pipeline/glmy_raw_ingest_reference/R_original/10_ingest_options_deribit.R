# ============================================================
# 10_ingest_options_deribit.R
# Parse Deribit BTC options trades and build bucketed IV return
# matrix for spillover/GLMY analysis.
#
# This is the "node construction" described in your draft:
#   node = maturity_bin × abs_log_moneyness_bin × {C,P}
#   aggregate to bars (default 15-min), compute weighted mean IV
#   then log-IV returns for spillover estimation.
# ============================================================

source(file.path("R", "00_packages.R"))
load_packages(c("data.table", "lubridate", "stringi"))
source(file.path("R", "parse_deribit_trades.R"))

# Floor timestamps to bar bins (fast, numeric)
floor_time_bin <- function(ts_posix, bar_minutes) {
  bsec <- as.numeric(bar_minutes) * 60
  as.POSIXct(floor(as.numeric(ts_posix) / bsec) * bsec,
             origin = "1970-01-01", tz = "UTC")
}

# Build bucketed node series from Deribit trades
build_options_bucket_panel <- function(
    trade_paths,
    bar_minutes = 15L,
    keep_top_nodes = 24L,
    expiry_breaks_days = c(0, 7, 30, 90, 180, 365, Inf),
    expiry_labels = c("0-7d","7-30d","30-90d","90-180d","180-365d","365d+"),
    mny_breaks_abslog = c(0, 0.02, 0.05, 0.10, Inf),
    mny_labels = c("ATM","Near","OTM","Far"),
    iv_floor = 1e-6
) {
  dt <- load_deribit_trade_csv(trade_paths)

  # dedup (helpful if files overlap)
  dt <- unique(dt, by = c("timestamp_ms","instrument_name","iv_dec","amount","index_price"))

  # Deribit options expire around 08:00 UTC (bucketing convention)
  dt[, expiry_dt := as.POSIXct(paste(expiry_date, "08:00:00"), tz = "UTC")]
  dt[, tte_days := as.numeric(difftime(expiry_dt, datetime_utc, units = "days"))]
  dt <- dt[is.finite(tte_days) & tte_days >= 0]

  dt[, spot := index_price]
  dt[, log_mny := log(strike / spot)]
  dt <- dt[is.finite(log_mny)]

  dt[, expiry_bucket := cut(tte_days, breaks = expiry_breaks_days, labels = expiry_labels,
                            include.lowest = TRUE, right = TRUE)]
  dt[, mny_bucket := cut(abs(log_mny), breaks = mny_breaks_abslog, labels = mny_labels,
                         include.lowest = TRUE, right = FALSE)]

  dt <- dt[!is.na(expiry_bucket) & !is.na(mny_bucket) & cp_flag %in% c("C","P")]
  dt[, node := paste(expiry_bucket, mny_bucket, cp_flag, sep = "|")]

  # Aggregate to bars: weighted mean IV + volume
  dt[, time_bin := floor_time_bin(datetime_utc, bar_minutes)]

  agg <- dt[, .(
    iv_w = weighted.mean(iv_dec, w = amount, na.rm = TRUE),
    vol = sum(amount, na.rm = TRUE),
    n_trades = .N
  ), by = .(time_bin, node)]

  # Select top nodes by total volume
  top_nodes <- agg[, .(totvol = sum(vol, na.rm = TRUE)), by = node][order(-totvol)]
  top_nodes <- top_nodes[1:min(keep_top_nodes, .N), node]
  agg <- agg[node %in% top_nodes]

  all_times <- seq(min(agg$time_bin), max(agg$time_bin), by = paste(bar_minutes, "mins"))
  panel <- merge(
    data.table::CJ(time_bin = all_times, node = sort(unique(agg$node))),
    agg, by = c("time_bin", "node"), all.x = TRUE
  )
  panel[is.na(vol), vol := 0]

  # wide IV (time x node)
  iv_wide <- data.table::dcast(panel, time_bin ~ node, value.var = "iv_w")
  data.table::setorder(iv_wide, time_bin)

  node_cols <- setdiff(names(iv_wide), "time_bin")

  # fill missing IV: LOCF -> NOCB -> median
  for (cc in node_cols) {
    v <- iv_wide[[cc]]
    v <- data.table::nafill(v, type = "locf")
    v <- data.table::nafill(v, type = "nocb")
    if (anyNA(v)) v[is.na(v)] <- median(v, na.rm = TRUE)
    iv_wide[[cc]] <- v
  }

  X <- as.matrix(iv_wide[, ..node_cols])
  X <- log(pmax(X, iv_floor))
  R <- apply(X, 2, diff)  # log-IV returns

  list(
    times = iv_wide$time_bin[-1],
    nodes = node_cols,
    returns = R,
    raw_panel = iv_wide
  )
}

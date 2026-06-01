# Build the hourly BTC--Deribit panel for the April 2025--March 2026 application.

source(file.path(dirname(sys.frame(1)$ofile), "00_config.R"))

emp2025_discover_deribit_files <- function(cfg) {
  files <- list.files(cfg$deribit_dir, pattern = "^[0-9]{8}-[0-9]{8}\\.csv$",
                      full.names = TRUE)
  if (length(files)) {
    bn <- basename(files)
    starts <- as.POSIXct(substr(bn, 1, 8), format = "%Y%m%d", tz = "UTC")
    ends <- as.POSIXct(substr(bn, 10, 17), format = "%Y%m%d", tz = "UTC") + 23 * 3600
    keep <- is.finite(starts) & is.finite(ends) &
      ends >= cfg$sample_start - 7 * 86400 & starts <= cfg$sample_end
    files <- files[keep]
  }
  extra <- cfg$extra_deribit_files
  if (length(extra)) {
    missing <- extra[!file.exists(extra)]
    if (length(missing)) stop("extra_deribit_files not found: ", paste(missing, collapse = ", "))
    files <- c(files, extra)
  }
  files <- unique(normalizePath(sort(files), winslash = "/", mustWork = TRUE))
  if (!length(files)) stop("No Deribit files found for requested sample window.")
  files
}

emp2025_floor_hour <- function(x) {
  as.POSIXct(floor(as.numeric(x) / 3600) * 3600, origin = "1970-01-01", tz = "UTC")
}

emp2025_panel_feature_cols <- function() {
  list(
    x_cols = c("ret_1", "ret_mom", "rv_ann_24", "range_1", "log_quote_volume"),
    z_cols = c("iv_rv_spread", "skew_proxy", "term_slope", "put_call_imbalance", "activity_log"),
    amount_cols = c("trade_count", "total_amount", "near_amount", "near_put_amount", "near_call_amount"),
    iv_cols = c("near_atm_iv", "mid_atm_iv", "near_otm_put_iv", "near_otm_call_iv",
                "near_index_price", "all_index_price")
  )
}

emp2025_finalize_panel <- function(panel) {
  cols <- emp2025_panel_feature_cols()
  keep <- c("hour_end", "open", "high", "low", "close", cols$x_cols, "y", "asset_simple_ret",
            cols$amount_cols, cols$iv_cols, cols$z_cols)
  panel <- panel[, ..keep]
  panel <- panel[stats::complete.cases(panel[, c("y", "asset_simple_ret", cols$x_cols, cols$z_cols), with = FALSE])]
  data.table::setorder(panel, hour_end)
  attr(panel, "x_cols") <- cols$x_cols
  attr(panel, "z_cols") <- cols$z_cols
  panel
}

emp2025_read_panel_csv <- function(path) {
  emp2025_require_packages(c("data.table"))
  panel <- data.table::fread(path)
  if (!"hour_end" %in% names(panel)) stop("panel_csv must contain an hour_end column: ", path)
  panel[, hour_end := as.POSIXct(hour_end, tz = "UTC")]
  emp2025_finalize_panel(panel)
}

emp2025_build_synthetic_panel <- function(cfg) {
  emp2025_require_packages(c("data.table"))
  set.seed(cfg$seed + 7L)
  n <- max(cfg$training_size + cfg$monitor_size + 4L * cfg$refit_every + 72L, 900L)
  hour_end <- seq(cfg$sample_start, by = "hour", length.out = n)
  innov <- stats::rnorm(n + 1L, sd = 0.012)
  ret_all <- as.numeric(stats::filter(innov, filter = 0.25, method = "recursive"))
  close_all <- 50000 * exp(cumsum(ret_all))
  ret_1 <- ret_all[seq_len(n)]
  close <- close_all[seq_len(n)]
  open <- c(close[1], head(close, -1L))
  high <- pmax(open, close) * exp(abs(stats::rnorm(n, sd = 0.003)))
  low <- pmin(open, close) * exp(-abs(stats::rnorm(n, sd = 0.003)))
  quote_volume <- exp(18 + stats::rnorm(n, sd = 0.35))
  latent_vol <- pmax(0.35, sqrt(365 * 24) * data.table::frollapply(ret_1, 24L, stats::sd, align = "right"))
  latent_vol[!is.finite(latent_vol)] <- stats::median(latent_vol, na.rm = TRUE)
  near_atm_iv <- latent_vol + stats::rnorm(n, sd = 0.05)
  mid_atm_iv <- near_atm_iv + 0.02 + stats::rnorm(n, sd = 0.02)
  near_otm_put_iv <- near_atm_iv + 0.06 + stats::rnorm(n, sd = 0.025)
  near_otm_call_iv <- near_atm_iv + 0.03 + stats::rnorm(n, sd = 0.025)
  total_amount <- exp(6 + 0.8 * abs(ret_1) / stats::sd(ret_1) + stats::rnorm(n, sd = 0.4))
  near_put_amount <- total_amount * stats::runif(n, 0.25, 0.55)
  near_call_amount <- total_amount * stats::runif(n, 0.25, 0.55)
  panel <- data.table::data.table(
    hour_end = hour_end, open = open, high = high, low = low, close = close,
    ret_1 = ret_1,
    asset_simple_ret = close_all[2:(n + 1L)] / close_all[seq_len(n)] - 1,
    y = ret_all[2:(n + 1L)],
    ret_mom = data.table::frollsum(ret_1, n = 24L, align = "right"),
    rv_ann_24 = sqrt(365 * 24) * data.table::frollapply(ret_1, n = 24L, FUN = stats::sd, align = "right"),
    range_1 = log(high / low),
    log_quote_volume = log1p(quote_volume),
    trade_count = stats::rpois(n, lambda = 80),
    total_amount = total_amount,
    near_amount = total_amount * stats::runif(n, 0.65, 0.95),
    near_put_amount = near_put_amount,
    near_call_amount = near_call_amount,
    near_atm_iv = near_atm_iv,
    mid_atm_iv = mid_atm_iv,
    near_otm_put_iv = near_otm_put_iv,
    near_otm_call_iv = near_otm_call_iv,
    near_index_price = close,
    all_index_price = close
  )
  panel[, put_call_imbalance := (near_call_amount - near_put_amount) / (near_call_amount + near_put_amount)]
  panel[, activity_log := log1p(total_amount)]
  panel[, skew_proxy := near_otm_put_iv - near_otm_call_iv]
  panel[, term_slope := mid_atm_iv - near_atm_iv]
  panel[, iv_rv_spread := near_atm_iv - rv_ann_24]
  emp2025_finalize_panel(panel)
}

emp2025_fetch_binance_klines <- function(cfg) {
  emp2025_require_packages(c("data.table", "jsonlite"))
  cache <- file.path(cfg$output_dir, "cache", "binance_btcusdt_1h_2025.csv")
  if (file.exists(cache) && !isTRUE(cfg$force)) return(data.table::fread(cache))

  start_fetch <- cfg$sample_start - 48 * 3600
  end_fetch <- cfg$sample_end + 2 * 3600
  start_ms <- as.numeric(start_fetch) * 1000
  end_ms <- as.numeric(end_fetch) * 1000
  endpoint <- "https://data-api.binance.vision/api/v3/klines"

  out <- list()
  cur <- start_ms
  i <- 1L
  while (cur <= end_ms) {
    url <- paste0(endpoint, "?symbol=", cfg$symbol, "&interval=1h&limit=1000",
                  "&startTime=", format(cur, scientific = FALSE, trim = TRUE),
                  "&endTime=", format(end_ms, scientific = FALSE, trim = TRUE))
    raw <- jsonlite::fromJSON(url)
    if (length(raw) == 0L) break
    dt <- data.table::as.data.table(raw)
    out[[i]] <- dt
    last_open <- as.numeric(dt[[1]][nrow(dt)])
    next_cur <- last_open + 3600 * 1000
    if (!is.finite(next_cur) || next_cur <= cur) break
    cur <- next_cur
    i <- i + 1L
    Sys.sleep(0.05)
  }
  if (!length(out)) stop("No Binance hourly data returned.")

  x <- data.table::rbindlist(out, fill = TRUE)
  data.table::setnames(x, names(x)[1:12], c(
    "open_time_ms", "open", "high", "low", "close", "volume",
    "close_time_ms", "quote_volume", "n_trades",
    "taker_buy_base", "taker_buy_quote", "ignore"
  ))
  num_cols <- setdiff(names(x), "ignore")
  for (cc in num_cols) x[, (cc) := as.numeric(get(cc))]
  x[, hour_end := as.POSIXct(close_time_ms / 1000, origin = "1970-01-01", tz = "UTC")]
  x[, hour_end := emp2025_floor_hour(hour_end)]
  x <- unique(x, by = "hour_end")
  data.table::setorder(x, hour_end)
  emp2025_write_csv(x, cache)
  x
}

emp2025_build_spot_features <- function(spot, cfg) {
  spot <- data.table::copy(spot)
  data.table::setorder(spot, hour_end)
  spot[, ret_1 := log(close / data.table::shift(close, 1))]
  spot[, asset_simple_ret := data.table::shift(close, type = "lead") / close - 1]
  spot[, y := data.table::shift(ret_1, type = "lead")]
  spot[, ret_mom := data.table::frollsum(ret_1, n = 24, align = "right", na.rm = FALSE)]
  spot[, rv_ann_24 := sqrt(365 * 24) * data.table::frollapply(ret_1, n = 24, FUN = stats::sd, align = "right")]
  spot[, range_1 := log(high / low)]
  spot[, log_quote_volume := log1p(quote_volume)]

  end_predictor <- cfg$sample_end
  if (isTRUE(cfg$keep_outcomes_inside_sample)) end_predictor <- cfg$sample_end - 3600
  spot[hour_end >= cfg$sample_start & hour_end <= end_predictor]
}

emp2025_aggregate_deribit_one <- function(file, cfg) {
  emp2025_require_packages(c("data.table", "lubridate", "stringi"))
  parse_file <- emp2025_path("glmy_raw_ingest_reference", "R_original", "parse_deribit_trades.R")
  source(parse_file, local = TRUE)
  dt <- load_deribit_trade_csv(file, auto_discover = FALSE)
  dt <- dt[datetime_utc >= cfg$sample_start - 7 * 86400 & datetime_utc <= cfg$sample_end]
  if (!nrow(dt)) return(data.table::data.table())

  dt[, hour_end := emp2025_floor_hour(datetime_utc)]
  dt[, expiry_dt := as.POSIXct(paste(expiry_date, "08:00:00"), tz = "UTC")]
  dt[, tte_days := as.numeric(difftime(expiry_dt, datetime_utc, units = "days"))]
  dt <- dt[is.finite(tte_days) & tte_days >= 0 & is.finite(index_price) & index_price > 0]
  dt[, log_mny := log(strike / index_price)]
  dt <- dt[is.finite(log_mny)]

  dt[, near := tte_days <= 30]
  dt[, mid := tte_days > 30 & tte_days <= 90]
  dt[, atm := abs(log_mny) <= 0.05]
  dt[, otm_put := cp_flag == "P" & log_mny < -0.03 & log_mny >= -0.25]
  dt[, otm_call := cp_flag == "C" & log_mny > 0.03 & log_mny <= 0.25]

  wmean_num <- function(x, w) sum(x * w, na.rm = TRUE)
  wsum <- function(w) sum(w, na.rm = TRUE)

  out <- dt[, .(
    trade_count = .N,
    total_amount = sum(amount, na.rm = TRUE),
    near_amount = sum(amount[near], na.rm = TRUE),
    near_put_amount = sum(amount[near & cp_flag == "P"], na.rm = TRUE),
    near_call_amount = sum(amount[near & cp_flag == "C"], na.rm = TRUE),
    near_atm_iv_num = wmean_num(iv_dec[near & atm], amount[near & atm]),
    near_atm_iv_den = wsum(amount[near & atm]),
    mid_atm_iv_num = wmean_num(iv_dec[mid & atm], amount[mid & atm]),
    mid_atm_iv_den = wsum(amount[mid & atm]),
    near_otm_put_iv_num = wmean_num(iv_dec[near & otm_put], amount[near & otm_put]),
    near_otm_put_iv_den = wsum(amount[near & otm_put]),
    near_otm_call_iv_num = wmean_num(iv_dec[near & otm_call], amount[near & otm_call]),
    near_otm_call_iv_den = wsum(amount[near & otm_call]),
    near_index_price_num = wmean_num(index_price[near], amount[near]),
    near_index_price_den = wsum(amount[near]),
    all_index_price_num = wmean_num(index_price, amount),
    all_index_price_den = wsum(amount)
  ), by = hour_end]
  out[, source_file := basename(file)]
  out
}

emp2025_build_deribit_hourly <- function(cfg) {
  cache <- file.path(cfg$output_dir, "cache", "deribit_hourly_2025.csv")
  if (file.exists(cache) && !isTRUE(cfg$force)) return(data.table::fread(cache))
  files <- emp2025_discover_deribit_files(cfg)
  if (isTRUE(cfg$smoke)) files <- files[1]
  parts <- lapply(files, emp2025_aggregate_deribit_one, cfg = cfg)
  agg <- data.table::rbindlist(parts, fill = TRUE)
  if (!nrow(agg)) stop("No Deribit rows after aggregation.")
  if ("source_file" %in% names(agg)) {
    data.table::setorder(agg, hour_end, -trade_count, source_file)
    agg <- agg[, .SD[1L], by = hour_end]
    agg[, source_file := NULL]
  }
  sum_cols <- setdiff(names(agg), "hour_end")
  agg <- agg[, lapply(.SD, sum, na.rm = TRUE), by = hour_end, .SDcols = sum_cols]
  ratio <- function(num, den) ifelse(is.finite(den) & den > 0, num / den, NA_real_)
  agg[, near_atm_iv := ratio(near_atm_iv_num, near_atm_iv_den)]
  agg[, mid_atm_iv := ratio(mid_atm_iv_num, mid_atm_iv_den)]
  agg[, near_otm_put_iv := ratio(near_otm_put_iv_num, near_otm_put_iv_den)]
  agg[, near_otm_call_iv := ratio(near_otm_call_iv_num, near_otm_call_iv_den)]
  agg[, near_index_price := ratio(near_index_price_num, near_index_price_den)]
  agg[, all_index_price := ratio(all_index_price_num, all_index_price_den)]
  keep <- c("hour_end", "trade_count", "total_amount", "near_amount", "near_put_amount",
            "near_call_amount", "near_atm_iv", "mid_atm_iv", "near_otm_put_iv",
            "near_otm_call_iv", "near_index_price", "all_index_price")
  agg <- agg[, ..keep]
  data.table::setorder(agg, hour_end)
  emp2025_write_csv(agg, cache)
  agg
}

emp2025_build_panel <- function(cfg) {
  emp2025_prepare_dirs(cfg)
  panel_path <- file.path(cfg$output_dir, "panel", "panel_hourly_2025.csv")
  if (nzchar(cfg$panel_csv)) return(emp2025_read_panel_csv(cfg$panel_csv))
  if (isTRUE(cfg$smoke) && isTRUE(cfg$smoke_synthetic)) {
    panel <- emp2025_build_synthetic_panel(cfg)
    emp2025_write_csv(panel, panel_path)
    cols <- emp2025_panel_feature_cols()
    emp2025_write_csv(data.table::data.table(x_cols = cols$x_cols), file.path(cfg$output_dir, "panel", "x_cols.csv"))
    emp2025_write_csv(data.table::data.table(z_cols = cols$z_cols), file.path(cfg$output_dir, "panel", "z_cols.csv"))
    return(panel)
  }
  if (file.exists(panel_path) && !isTRUE(cfg$force)) return(data.table::fread(panel_path))

  spot <- emp2025_build_spot_features(emp2025_fetch_binance_klines(cfg), cfg)
  opt <- emp2025_build_deribit_hourly(cfg)

  all_hours <- data.table::data.table(
    hour_end = seq(cfg$sample_start, if (cfg$keep_outcomes_inside_sample) cfg$sample_end - 3600 else cfg$sample_end,
                   by = "hour")
  )
  panel <- merge(all_hours, spot, by = "hour_end", all.x = TRUE)
  panel <- merge(panel, opt, by = "hour_end", all.x = TRUE)
  data.table::setorder(panel, hour_end)

  amount_cols <- c("trade_count", "total_amount", "near_amount", "near_put_amount", "near_call_amount")
  for (cc in amount_cols) panel[is.na(get(cc)), (cc) := 0]
  iv_cols <- c("near_atm_iv", "mid_atm_iv", "near_otm_put_iv", "near_otm_call_iv",
               "near_index_price", "all_index_price")
  for (cc in iv_cols) panel[, (cc) := data.table::nafill(get(cc), type = "locf")]

  panel[, put_call_imbalance := data.table::fifelse(
    near_put_amount + near_call_amount > 0,
    (near_call_amount - near_put_amount) / (near_call_amount + near_put_amount),
    NA_real_
  )]
  panel[, activity_log := log1p(total_amount)]
  panel[, trade_count_log := log1p(trade_count)]
  panel[, skew_proxy := near_otm_put_iv - near_otm_call_iv]
  panel[, term_slope := mid_atm_iv - near_atm_iv]
  panel[, iv_rv_spread := near_atm_iv - rv_ann_24]

  panel <- emp2025_finalize_panel(panel)
  cols <- emp2025_panel_feature_cols()
  emp2025_write_csv(panel, panel_path)
  emp2025_write_csv(data.table::data.table(x_cols = cols$x_cols), file.path(cfg$output_dir, "panel", "x_cols.csv"))
  emp2025_write_csv(data.table::data.table(z_cols = cols$z_cols), file.path(cfg$output_dir, "panel", "z_cols.csv"))
  panel
}

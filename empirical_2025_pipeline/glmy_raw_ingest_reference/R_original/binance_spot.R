# Binance spot price helper (public API)
#
# Uses:
#   GET https://api.binance.com/api/v3/klines?symbol=BTCUSDT&interval=1d&startTime=...&endTime=...&limit=1000
# Docs: https://developers.binance.com/docs/binance-spot-api-docs/rest-api/market-data-endpoints

get_binance_daily_close <- function(start_date,
                                   end_date,
                                   cache_path = "data/btc_spot_binance.csv",
                                   symbol = "BTCUSDT",
                                   interval = "1d",
                                   # Binance provides a dedicated endpoint for public market data.
                                   # If this is blocked in your environment, you can override via the endpoint arg.
                                   endpoint = "https://data-api.binance.vision/api/v3/klines",
                                   verbose = TRUE) {

  stopifnot(inherits(start_date, "Date"), inherits(end_date, "Date"))
  if (end_date < start_date) stop("end_date must be >= start_date")

  # If cached file exists and covers range, reuse
  if (file.exists(cache_path)) {
    cached <- tryCatch(data.table::fread(cache_path), error = function(e) NULL)
    if (!is.null(cached) && all(c("date", "close") %in% names(cached))) {
      cached[, date := as.Date(date)]
      if (min(cached$date, na.rm = TRUE) <= start_date && max(cached$date, na.rm = TRUE) >= end_date) {
        out <- cached[date >= start_date & date <= end_date]
        data.table::setorder(out, date)
        return(out)
      }
    }
  }

  # Fetch from Binance in (up to) 1000-day chunks
  start_ms <- as.numeric(as.POSIXct(start_date, tz = "UTC")) * 1000
  end_ms   <- as.numeric(as.POSIXct(end_date + 1, tz = "UTC")) * 1000 - 1

  out_list <- list()
  i <- 1
  cur_ms <- start_ms
  repeat {
    if (cur_ms > end_ms) break
    if (verbose) message(sprintf("  Binance request %d (from %s)...", i, as.character(as.POSIXct(cur_ms/1000, origin = "1970-01-01", tz = "UTC"))))

    resp <- httr::GET(
      url = endpoint,
      query = list(
        symbol = symbol,
        interval = interval,
        startTime = format(cur_ms, scientific = FALSE, trim = TRUE),
        endTime = format(end_ms, scientific = FALSE, trim = TRUE),
        limit = 1000
      ),
      httr::user_agent("relevant-hypothesis-btc-empirical/1.0")
    )

    if (httr::status_code(resp) != 200) {
      msg <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
      stop(sprintf("Binance API request failed (HTTP %s). Response: %s", httr::status_code(resp), msg))
    }

    txt_json <- httr::content(resp, as = "text", encoding = "UTF-8")
    dat <- jsonlite::fromJSON(txt_json, simplifyVector = TRUE)
    if (length(dat) == 0) break

    # Binance returns a list of arrays. Columns:
    # 0 open time, 1 open, 2 high, 3 low, 4 close, 5 volume, 6 close time, ...
    # jsonlite may simplify Binance's array-of-arrays into a matrix/data.frame; handle both cases.
    m <- if (is.list(dat) && !is.data.frame(dat)) do.call(rbind, dat) else dat
    m <- data.table::as.data.table(m)
    if (ncol(m) < 7) stop("Unexpected Binance klines payload shape")

    tmp <- data.table::data.table(
      open_time_ms  = as.numeric(m[[1]]),
      open          = as.numeric(m[[2]]),
      high          = as.numeric(m[[3]]),
      low           = as.numeric(m[[4]]),
      close         = as.numeric(m[[5]]),
      volume        = as.numeric(m[[6]]),
      close_time_ms = as.numeric(m[[7]])
    )
    tmp[, date := as.Date(as.POSIXct(open_time_ms/1000, origin = "1970-01-01", tz = "UTC"))]
    out_list[[i]] <- tmp

    last_open <- max(tmp$open_time_ms, na.rm = TRUE)
    if (!is.finite(last_open)) break
    cur_ms <- last_open + 1
    i <- i + 1

    # Safety: if Binance returns less than limit and last date already beyond end
    if (nrow(tmp) < 1000 && max(tmp$date, na.rm = TRUE) >= end_date) break
  }

  out <- data.table::rbindlist(out_list, use.names = TRUE, fill = TRUE)
  out <- unique(out, by = "open_time_ms")
  out <- out[date >= start_date & date <= end_date]
  data.table::setorder(out, date)

  # Save cache (full fetched range)
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(out[, .(date, close, open, high, low, volume, open_time_ms, close_time_ms)], cache_path)

  return(out[, .(date, close, open, high, low, volume, open_time_ms, close_time_ms)])
}

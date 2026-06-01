# Parse Deribit-style BTC options trade CSV (s,t,trade)
#
# The CSV you provided has:
#   - t: millisecond timestamp (also duplicated inside the JSON)
#   - trade: JSON string with fields like instrument_name, iv, amount, index_price
#
# This loader:
#   - reads the CSV with data.table::fread
#   - extracts only the fields we need using fast regex (stringi)
#   - parses instrument_name like "BTC-30JAN26-80000-P" into expiry/strike/type
#
# UPDATE (extended sample):
#   load_deribit_trade_csv() now supports combining multiple files seamlessly.
#   - path can be a CSV file, a directory, a ZIP file, or a vector of these.
#   - If a SINGLE CSV/ZIP is provided, the loader will also auto-discover
#     other trade CSV/ZIP files in the same folder and bind them together.

parse_instrument_name <- function(instr) {
  # Expected: BTC-30JAN26-80000-P
  #           BTC-16JAN26-85000-C
  # Returns data.table with columns: expiry_date, strike, cp_flag

  # Capture groups: day, mon, year, strike, flag
  m <- stringi::stri_match_first_regex(
    instr,
    "^BTC-([0-9]{1,2})([A-Z]{3})([0-9]{2})-([0-9]+(?:\\.[0-9]+)?)-([CP])$"
  )
  day   <- as.integer(m[, 2])
  mon_s <- m[, 3]
  yy    <- as.integer(m[, 4])
  strike <- as.numeric(m[, 5])
  cp    <- m[, 6]

  mon_map <- c(JAN=1, FEB=2, MAR=3, APR=4, MAY=5, JUN=6,
               JUL=7, AUG=8, SEP=9, OCT=10, NOV=11, DEC=12)
  mon <- unname(mon_map[mon_s])
  # assume 20xx
  year <- 2000L + yy
  expiry_date <- as.Date(sprintf("%04d-%02d-%02d", year, mon, day), tz = "UTC")

  data.table::data.table(
    expiry_date = expiry_date,
    strike = strike,
    cp_flag = cp
  )
}


# -----------------------------
# Helpers: multi-file / zip IO
# -----------------------------

.is_zip <- function(p) {
  is.character(p) && length(p) == 1 && grepl("\\.zip$", p, ignore.case = TRUE)
}

.is_csv <- function(p) {
  is.character(p) && length(p) == 1 && grepl("\\.csv$", p, ignore.case = TRUE)
}

.is_trade_csv <- function(csv_path, keep_cols = c("t", "trade")) {
  dt0 <- tryCatch(
    data.table::fread(csv_path, nrows = 0, showProgress = FALSE),
    error = function(e) NULL
  )
  !is.null(dt0) && all(keep_cols %in% names(dt0))
}

.zip_list_csv_members <- function(zip_path) {
  zlist <- tryCatch(utils::unzip(zip_path, list = TRUE), error = function(e) NULL)
  if (is.null(zlist) || nrow(zlist) == 0) return(character())

  mem <- zlist$Name
  mem <- mem[grepl("\\.csv$", mem, ignore.case = TRUE)]
  mem <- mem[!grepl("^__MACOSX/", mem)]
  mem <- mem[!grepl("/\\._", mem)]
  mem
}

.fread_zip_member <- function(zip_path, member, select, nrows) {
  # Try streaming via system unzip (fast, no temp file)
  cmd <- sprintf("unzip -p %s %s", shQuote(zip_path), shQuote(member))
  dt <- tryCatch(
    data.table::fread(cmd = cmd, select = select, nrows = nrows, showProgress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(dt)) return(dt)

  # Fallback: extract to tempdir (portable)
  tmpdir <- tempfile(pattern = "unz_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  utils::unzip(zip_path, files = member, exdir = tmpdir)
  extracted <- file.path(tmpdir, member)

  if (!file.exists(extracted)) {
    found <- list.files(tmpdir, recursive = TRUE, full.names = TRUE,
                        pattern = paste0(basename(member), "$"))
    if (length(found) == 0) {
      stop("Failed to extract '", member, "' from zip: ", zip_path)
    }
    extracted <- found[1]
  }

  data.table::fread(extracted, select = select, nrows = nrows, showProgress = FALSE)
}

.zip_trade_members <- function(zip_path, keep_cols = c("t", "trade")) {
  mem <- .zip_list_csv_members(zip_path)
  if (length(mem) == 0) return(character())

  good <- logical(length(mem))
  for (i in seq_along(mem)) {
    m <- mem[i]
    dt0 <- tryCatch(
      .fread_zip_member(zip_path, m, select = NULL, nrows = 0),
      error = function(e) NULL
    )
    good[i] <- !is.null(dt0) && all(keep_cols %in% names(dt0))
  }
  mem[good]
}

.resolve_trade_sources <- function(path,
                                  keep_cols = c("t", "trade"),
                                  auto_discover = TRUE) {
  paths <- as.character(path)
  paths <- paths[nzchar(paths)]

  # Allow multiple inputs via a single string (useful for Sys.setenv(OPTIONS_CSV=...))
  # Example: "data/options_trades.csv;data/part-00000-...csv.zip"
  if (length(paths) == 1 && !file.exists(paths) && !dir.exists(paths) && grepl("[,;]", paths)) {
    pieces <- trimws(unlist(strsplit(paths, "[,;]")))
    pieces <- pieces[nzchar(pieces)]
    if (length(pieces) > 0) paths <- pieces
  }
  if (length(paths) == 0) stop("Empty 'path' provided to load_deribit_trade_csv().")

  # Collect candidate CSV/ZIP files.
  candidates <- character()
  for (p in paths) {
    if (dir.exists(p)) {
      candidates <- c(candidates,
                      list.files(p, pattern = "\\.(csv|zip)$", full.names = TRUE, ignore.case = TRUE))
    } else if (file.exists(p)) {
      candidates <- c(candidates, p)

      # If the user passed a single file, also load other trade files in that same folder.
      if (auto_discover && length(paths) == 1) {
        candidates <- c(candidates,
                        list.files(dirname(p), pattern = "\\.(csv|zip)$", full.names = TRUE, ignore.case = TRUE))
      }
    } else {
      stop("Options CSV/ZIP (or directory) not found: ", p)
    }
  }

  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = TRUE))
  candidates <- sort(candidates)

  sources <- list()
  for (f in candidates) {
    if (.is_zip(f)) {
      mem <- .zip_trade_members(f, keep_cols = keep_cols)
      if (length(mem) > 0) {
        for (m in mem) {
          sources[[length(sources) + 1]] <- list(type = "zip", zip = f, member = m)
        }
      }
    } else if (.is_csv(f)) {
      if (.is_trade_csv(f, keep_cols = keep_cols)) {
        sources[[length(sources) + 1]] <- list(type = "csv", file = f)
      }
    }
  }

  if (length(sources) == 0) {
    stop(
      "No Deribit-style trade CSVs found in: ", paste(paths, collapse = ", "),
      "\nExpected a CSV with columns: ", paste(keep_cols, collapse = ", ")
    )
  }
  sources
}

.parse_trade_dt <- function(dt) {
  # Some exports (e.g. Spark) escape inner quotes as \" inside a quoted CSV field.
  # data.table::fread can often handle this automatically, but we defensively
  # unescape \" -> " so the JSON keys match our regex extraction.
  if ("trade" %in% names(dt)) {
    dt[, trade := gsub('\\\"', '"', trade, fixed = TRUE)]
  }

  # timestamp in ms (prefer column t)
  dt[, timestamp_ms := as.numeric(t)]
  dt[, datetime_utc := as.POSIXct(timestamp_ms / 1000, origin = "1970-01-01", tz = "UTC")]
  dt[, date := as.Date(datetime_utc, tz = "UTC")]

  # Fast regex extraction from the JSON string in column 'trade'
  # We deliberately avoid jsonlite::fromJSON row-by-row (too slow for large files).
  tr <- dt[["trade"]]

  # instrument_name
  instr <- stringi::stri_match_first_regex(tr, '"instrument_name":"([^"]+)"')[, 2]
  iv    <- stringi::stri_match_first_regex(tr, '"iv":([0-9.]+)')[, 2]
  amt   <- stringi::stri_match_first_regex(tr, '"amount":([0-9.]+)')[, 2]
  idxp  <- stringi::stri_match_first_regex(tr, '"index_price":([0-9.]+)')[, 2]

  dt[, instrument_name := instr]
  dt[, iv := as.numeric(iv)]
  dt[, amount := as.numeric(amt)]
  dt[, index_price := as.numeric(idxp)]

  # parse instrument
  parsed <- parse_instrument_name(dt$instrument_name)
  dt <- cbind(dt, parsed)

  # basic cleaning
  dt <- dt[is.finite(iv) & is.finite(strike) & !is.na(expiry_date)]
  dt[, iv_dec := iv / 100]  # Deribit-style IV is typically in percent
  dt <- dt[iv_dec > 0 & iv_dec < 5]
  dt[, amount := data.table::fifelse(is.finite(amount) & amount > 0, amount, 1.0)]

  # keep only necessary columns
  dt <- dt[, .(
    date,
    datetime_utc,
    timestamp_ms,
    instrument_name,
    expiry_date,
    strike,
    cp_flag,
    iv_dec,
    amount,
    index_price
  )]

  data.table::setorder(dt, date)
  dt
}


load_deribit_trade_csv <- function(path,
                                   max_rows = Inf,
                                   keep_cols = c("t", "trade"),
                                   auto_discover = TRUE) {
  sources <- .resolve_trade_sources(path, keep_cols = keep_cols, auto_discover = auto_discover)

  message(sprintf("Found %d trade source(s).", length(sources)))

  dts <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    src <- sources[[i]]

    nrows <- if (is.finite(max_rows)) max_rows else -1

    if (identical(src$type, "csv")) {
      message("Reading: ", src$file)
      raw <- data.table::fread(
        src$file,
        select = keep_cols,
        nrows = nrows,
        showProgress = FALSE
      )
    } else if (identical(src$type, "zip")) {
      message("Reading from zip: ", src$zip, " :: ", src$member)
      raw <- .fread_zip_member(
        zip_path = src$zip,
        member   = src$member,
        select   = keep_cols,
        nrows    = nrows
      )
    } else {
      next
    }

    dts[[i]] <- .parse_trade_dt(raw)
  }

  dt <- data.table::rbindlist(dts, use.names = TRUE, fill = TRUE)
  data.table::setorder(dt, date)
  dt
}

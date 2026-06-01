# Trading-overlay summaries and robustness diagnostics for the 2025 empirical run.

source(file.path(dirname(sys.frame(1)$ofile), "00_config.R"))
emp2025_require_packages(c("data.table"))

emp2025_strategy_returns <- function(panel, alarm_table, method, side, hold, cost_bps = 0,
                                     short_funding_bps_per_hour = 0) {
  n <- nrow(panel)
  pos <- rep(1, n)
  mth <- method
  alarms <- alarm_table[method == mth & !is.na(trade_start_idx), trade_start_idx]
  alarms <- alarms[alarms >= 1 & alarms <= n]
  if (length(alarms)) {
    for (a in alarms) {
      end <- min(n, a + hold - 1L)
      pos[a:end] <- if (side == "flat") 0 else -1
    }
  }
  turnover_vec <- c(0, abs(diff(pos)))
  cost <- (cost_bps / 10000) * turnover_vec
  funding <- if (side == "short") (short_funding_bps_per_hour / 10000) * as.numeric(pos < 0) else 0
  ret <- pos * panel$asset_simple_ret - cost - funding
  data.table::data.table(hour_end = panel$hour_end, method = method, strategy_side = side,
                         hold_period = hold, position = pos, strategy_return = ret,
                         benchmark_return = panel$asset_simple_ret,
                         turnover = turnover_vec, cost = cost, short_funding = funding)
}

emp2025_perf <- function(r) {
  r <- as.numeric(r)
  gross <- cumprod(1 + data.table::fifelse(is.na(r), 0, r))
  total <- tail(gross, 1) - 1
  ann_factor <- 365 * 24
  mu <- mean(r, na.rm = TRUE) * ann_factor
  vol <- stats::sd(r, na.rm = TRUE) * sqrt(ann_factor)
  downside <- stats::sd(pmin(r, 0), na.rm = TRUE) * sqrt(ann_factor)
  dd <- gross / cummax(gross) - 1
  data.table::data.table(
    ann_return = mu, ann_vol = vol,
    sharpe = ifelse(vol > 0, mu / vol, NA_real_),
    sortino = ifelse(downside > 0, mu / downside, NA_real_),
    max_drawdown = min(dd, na.rm = TRUE),
    calmar = ifelse(min(dd, na.rm = TRUE) < 0, mu / abs(min(dd, na.rm = TRUE)), NA_real_),
    total_return = total, ending_gross = tail(gross, 1)
  )
}

emp2025_make_strategy_outputs <- function(panel, alarm_table, cfg,
                                          cost_bps = 0, short_funding_bps_per_hour = 0) {
  methods <- unique(alarm_table$method)
  rets <- list()
  rr <- 1L
  for (mm in methods) {
    for (side in cfg$strategy_sides) {
      for (hold in cfg$hold_periods) {
        rets[[rr]] <- emp2025_strategy_returns(panel, alarm_table, mm, side, hold,
                                               cost_bps, short_funding_bps_per_hour)
        rr <- rr + 1L
      }
    }
  }
  ret_long <- data.table::rbindlist(rets)
  summ <- ret_long[, c(emp2025_perf(strategy_return),
                       .(turnover = sum(turnover, na.rm = TRUE) / .N,
                         n_alarm_entries = alarm_table[method == .BY$method & !is.na(trade_start_idx), .N])),
                   by = .(method, strategy_side, hold_period)]
  list(strategy_returns = ret_long, strategy_summary = summ)
}

emp2025_alarm_quality <- function(panel, alarm_table, cfg) {
  event_cut <- stats::quantile(panel$asset_simple_ret, probs = 0.10, na.rm = TRUE)
  event <- panel$asset_simple_ret <= event_cut
  rows <- list()
  rr <- 1L
  for (mm in unique(alarm_table$method)) {
    for (side in cfg$strategy_sides) {
      for (hold in cfg$hold_periods) {
        alarms <- alarm_table[method == mm & !is.na(trade_start_idx), trade_start_idx]
        hit_hours <- rep(FALSE, nrow(panel))
        if (length(alarms)) for (a in alarms) hit_hours[a:min(nrow(panel), a + hold - 1L)] <- TRUE
        precision <- if (length(alarms)) mean(vapply(alarms, function(a) {
          any(event[a:min(nrow(panel), a + hold - 1L)], na.rm = TRUE)
        }, logical(1))) else NA_real_
        rows[[rr]] <- data.table::data.table(
          method = mm, strategy_side = side, hold_period = hold,
          precision = precision,
          event_capture = if (sum(event, na.rm = TRUE) > 0) sum(hit_hours & event, na.rm = TRUE) / sum(event, na.rm = TRUE) else NA_real_
        )
        rr <- rr + 1L
      }
    }
  }
  data.table::rbindlist(rows)
}

emp2025_no_lookahead_audit <- function(alarm_table) {
  alarm_table[, .(
    window_id, method,
    training_before_monitoring = train_end_idx < monitor_start_idx,
    alarm_inside_monitoring = is.na(alarm_abs_idx) | (alarm_abs_idx >= monitor_start_idx & alarm_abs_idx <= monitor_end_idx),
    trade_after_alarm = is.na(trade_start_idx) | trade_start_idx > alarm_abs_idx,
    train_end_time, monitor_start_time, monitor_end_time, alarm_time, trade_start_time
  )]
}

emp2025_alarm_sparsity <- function(panel, alarm_table, strategy_summary, quality) {
  out <- alarm_table[, .(
    n_windows_tested = data.table::uniqueN(window_id),
    n_alarm_windows = data.table::uniqueN(window_id[!is.na(alarm_abs_idx)]),
    n_alarm_entries = sum(!is.na(alarm_abs_idx)),
    median_first_alarm_hour = as.numeric(stats::median(alarm_k, na.rm = TRUE)),
    mean_stat_cv_ratio = mean(max_path_over_threshold, na.rm = TRUE),
    max_stat_cv_ratio = max(max_path_over_threshold, na.rm = TRUE)
  ), by = .(method, family, stat, gamma, weight)]
  flat <- strategy_summary[strategy_side == "flat" & hold_period == 24,
                           .(method, total_hours_flat_proxy = round(turnover * nrow(panel)))]
  short <- strategy_summary[strategy_side == "short" & hold_period == 24,
                            .(method, total_hours_short_proxy = round(turnover * nrow(panel)))]
  q <- quality[strategy_side == "flat" & hold_period == 24, .(method, event_capture)]
  Reduce(function(a, b) merge(a, b, by = "method", all.x = TRUE), list(out, flat, short, q))
}

emp2025_transaction_cost_sensitivity <- function(panel, alarm_table, cfg) {
  rows <- list()
  rr <- 1L
  for (cost in cfg$transaction_cost_bps) {
    for (fund in cfg$short_funding_bps_per_hour) {
      x <- emp2025_make_strategy_outputs(panel, alarm_table, cfg, cost, fund)$strategy_summary
      x[, `:=`(cost_bps = cost, short_funding_bps_per_hour = fund)]
      rows[[rr]] <- x
      rr <- rr + 1L
    }
  }
  data.table::rbindlist(rows)
}

emp2025_mc_p_ge <- function(vals, observed) {
  vals <- vals[is.finite(vals)]
  if (!length(vals) || !is.finite(observed)) return(NA_real_)
  (1 + sum(vals >= observed)) / (length(vals) + 1)
}

emp2025_random_alarm_placebo <- function(panel, alarm_table, cfg) {
  set.seed(cfg$seed + 19L)
  B <- if (isTRUE(cfg$smoke)) min(50L, cfg$n_random_alarm) else cfg$n_random_alarm
  valid <- seq_len(nrow(panel))
  observed <- emp2025_make_strategy_outputs(panel, alarm_table, cfg)$strategy_summary
  rows <- list()
  rr <- 1L
  for (ii in seq_len(nrow(observed))) {
    row <- observed[ii]
    n_alarm <- row$n_alarm_entries
    vals_ret <- vals_sortino <- vals_dd <- numeric(B)
    for (b in seq_len(B)) {
      fake <- if (n_alarm > 0) {
        data.table::data.table(method = row$method, trade_start_idx = sample(valid, n_alarm))
      } else {
        data.table::data.table(method = character(), trade_start_idx = integer())
      }
      fake_ret <- emp2025_strategy_returns(panel, fake, row$method, row$strategy_side, row$hold_period)
      perf <- emp2025_perf(fake_ret$strategy_return)
      vals_ret[b] <- perf$total_return
      vals_sortino[b] <- perf$sortino
      vals_dd[b] <- perf$max_drawdown
    }
    rows[[rr]] <- data.table::data.table(
      method = row$method, strategy_side = row$strategy_side, hold_period = row$hold_period,
      n_alarm_entries = n_alarm, B = B,
      observed_total_return = row$total_return,
      observed_sortino = row$sortino,
      observed_max_drawdown = row$max_drawdown,
      p_placebo_return = emp2025_mc_p_ge(vals_ret, row$total_return),
      p_placebo_sortino = emp2025_mc_p_ge(vals_sortino, row$sortino),
      p_placebo_drawdown = emp2025_mc_p_ge(vals_dd, row$max_drawdown),
      pvalue_correction = "(1 + exceedances) / (1 + B)"
    )
    rr <- rr + 1L
  }
  data.table::rbindlist(rows)
}

emp2025_reality_check_random_alarm <- function(panel, alarm_table, cfg, strategy_summary = NULL) {
  set.seed(cfg$seed + 23L)
  B <- if (isTRUE(cfg$smoke)) min(50L, cfg$n_random_alarm) else cfg$n_random_alarm
  observed <- if (is.null(strategy_summary)) {
    emp2025_make_strategy_outputs(panel, alarm_table, cfg)$strategy_summary
  } else {
    data.table::copy(strategy_summary)
  }
  keys <- observed[, .(method, strategy_side, hold_period, total_return, n_alarm_entries)]
  fam <- unique(alarm_table[, .(method, family)])
  keys <- merge(keys, fam, by = "method", all.x = TRUE)
  valid <- seq_len(nrow(panel))
  method_counts <- unique(keys[, .(method, n_alarm_entries)])
  rand <- matrix(NA_real_, nrow = B, ncol = nrow(keys))

  for (b in seq_len(B)) {
    fake_by_method <- vector("list", nrow(method_counts))
    names(fake_by_method) <- method_counts$method
    for (jj in seq_len(nrow(method_counts))) {
      mm <- method_counts$method[jj]
      n_alarm <- method_counts$n_alarm_entries[jj]
      fake_by_method[[mm]] <- if (n_alarm > 0) {
        data.table::data.table(method = mm, trade_start_idx = sample(valid, n_alarm))
      } else {
        data.table::data.table(method = character(), trade_start_idx = integer())
      }
    }
    for (ii in seq_len(nrow(keys))) {
      k <- keys[ii]
      fake_ret <- emp2025_strategy_returns(panel, fake_by_method[[k$method]], k$method,
                                           k$strategy_side, k$hold_period)
      rand[b, ii] <- emp2025_perf(fake_ret$strategy_return)$total_return
    }
  }

  rand_mean <- colMeans(rand, na.rm = TRUE)
  observed_excess <- keys$total_return - rand_mean
  make_scope <- function(label, idx) {
    max_obs <- max(observed_excess[idx], na.rm = TRUE)
    max_boot <- apply(sweep(rand[, idx, drop = FALSE], 2, rand_mean[idx], "-"), 1, max, na.rm = TRUE)
    data.table::data.table(
      family = label,
      n_comparisons = length(idx),
      observed_max_excess_return = max_obs,
      p_reality_check_return = emp2025_mc_p_ge(max_boot, max_obs),
      B = B,
      adjustment = "Shared matched-random-alarm max-statistic reality check"
    )
  }
  rows <- lapply(sort(unique(keys$family)), function(ff) make_scope(ff, which(keys$family == ff)))
  rows[[length(rows) + 1L]] <- make_scope("ALL", seq_len(nrow(keys)))
  data.table::rbindlist(rows, fill = TRUE)
}

emp2025_circular_boot <- function(x, block, n) {
  starts <- sample.int(length(x), ceiling(n / block), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) ((s - 1L + seq_len(block) - 1L) %% length(x)) + 1L))
  x[idx[seq_len(n)]]
}

emp2025_bootstrap_uncertainty <- function(strategy_returns, cfg) {
  set.seed(cfg$seed + 31L)
  B <- if (isTRUE(cfg$smoke)) min(50L, cfg$n_bootstrap) else cfg$n_bootstrap
  rows <- list()
  rr <- 1L
  keys <- unique(strategy_returns[, .(method, strategy_side, hold_period)])
  for (ii in seq_len(nrow(keys))) {
    k <- keys[ii]
    r <- strategy_returns[k, on = .(method, strategy_side, hold_period)]$strategy_return
    bench <- strategy_returns[k, on = .(method, strategy_side, hold_period)]$benchmark_return
    for (block in cfg$bootstrap_blocks) {
      vals <- replicate(B, {
        rb <- emp2025_circular_boot(r, block, length(r))
        bb <- emp2025_circular_boot(bench, block, length(bench))
        c(total_return = emp2025_perf(rb)$total_return,
          mean_hourly_excess_return = mean(rb - bb, na.rm = TRUE),
          sortino = emp2025_perf(rb)$sortino,
          max_drawdown = emp2025_perf(rb)$max_drawdown,
          return_difference_vs_benchmark = emp2025_perf(rb)$total_return - emp2025_perf(bb)$total_return)
      })
      qs <- t(apply(vals, 1, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
      rows[[rr]] <- data.table::data.table(k, block_length = block, metric = rownames(qs),
                                           ci_low = qs[, 1], median = qs[, 2], ci_high = qs[, 3], B = B)
      rr <- rr + 1L
    }
  }
  data.table::rbindlist(rows)
}

emp2025_multiple_testing <- function(random_alarm_placebo, strategy_summary, alarm_table,
                                     reality_check = NULL) {
  fam <- unique(alarm_table[, .(method, family)])
  rand <- merge(random_alarm_placebo, fam, by = "method", all.x = TRUE)
  p_rows <- rand[, .(
    n_comparisons = .N,
    min_matched_alarm_p = min(p_placebo_return, na.rm = TRUE),
    bonferroni_p_return = min(1, .N * min(p_placebo_return, na.rm = TRUE)),
    adjustment = "Bonferroni over matched-random-alarm p-values"
  ), by = family]
  all_p <- rand[, .(
    family = "ALL",
    n_comparisons = .N,
    min_matched_alarm_p = min(p_placebo_return, na.rm = TRUE),
    bonferroni_p_return = min(1, .N * min(p_placebo_return, na.rm = TRUE)),
    adjustment = "Bonferroni over all detector-rule p-values"
  )]
  out <- data.table::rbindlist(list(p_rows, all_p), fill = TRUE)
  if (!is.null(reality_check) && nrow(reality_check)) {
    rc <- reality_check[, .(family, reality_check_p_return = p_reality_check_return,
                            observed_max_excess_return, reality_check_B = B,
                            reality_check_adjustment = adjustment)]
    out <- merge(out, rc, by = "family", all.x = TRUE)
  }
  out
}

emp2025_make_placebo_z <- function(panel, type, cfg = NULL) {
  z_cols <- c("iv_rv_spread", "skew_proxy", "term_slope", "put_call_imbalance", "activity_log")
  if (type == "lagged_spot_only") {
    out <- panel[, .(spot_lag1 = ret_1, spot_mom = ret_mom, spot_rv = rv_ann_24,
                     spot_range = range_1, spot_volume = log_quote_volume)]
    return(out)
  }
  if (type == "delayed_option_168h") {
    out <- panel[, ..z_cols]
    for (cc in names(out)) out[, (cc) := data.table::shift(get(cc), 168L)]
    return(out)
  }
  if (type == "block_permuted_week") {
    seed <- if (!is.null(cfg) && !is.null(cfg$seed)) cfg$seed + 41L else 20250601L
    set.seed(seed)
    out <- data.table::copy(panel[, ..z_cols])
    week <- as.integer(as.Date(panel$hour_end) - min(as.Date(panel$hour_end))) %/% 7L
    split_idx <- split(seq_len(nrow(panel)), week)
    source_idx <- unlist(split_idx[sample(names(split_idx))], use.names = FALSE)
    source_idx <- source_idx[seq_len(nrow(panel))]
    return(out[source_idx])
  }
  stop("Unknown placebo type: ", type)
}

emp2025_summarize_placebo <- function(placebo_alarm) {
  placebo_alarm[, .(
    n_windows = data.table::uniqueN(window_id),
    n_alarm_windows = data.table::uniqueN(window_id[!is.na(alarm_abs_idx)]),
    median_alarm_time = as.numeric(suppressWarnings(stats::median(alarm_k, na.rm = TRUE))),
    mean_stat_cv_ratio = mean(max_path_over_threshold, na.rm = TRUE),
    max_stat_cv_ratio = max(max_path_over_threshold, na.rm = TRUE)
  ), by = .(placebo_type = sample_label, method, stat, gamma, weight)]
}

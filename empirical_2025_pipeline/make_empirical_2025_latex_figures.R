#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
parse_named_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  named <- grep("^--[A-Za-z0-9_.-]+=", args, value = TRUE)
  for (x in named) {
    key <- sub("^--([^=]+)=.*$", "\\1", x)
    val <- sub("^--[^=]+=", "", x)
    out[[key]] <- val
  }
  out
}
path_arg <- function(x, default, base = script_dir) {
  if (is.null(x) || !nzchar(x)) return(default)
  if (grepl("^/", x)) return(normalizePath(x, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(base, x), winslash = "/", mustWork = FALSE)
}
args <- parse_named_args()
out <- path_arg(args$output_dir, file.path(root, "empirical_2025_pipeline", "output_apr2025_mar2026_ncv20000"))
tables_dir <- path_arg(args$tables_dir, file.path(root, "tables"))
fig_dir <- path_arg(args$fig_dir, file.path(root, "figures"))
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

panel <- fread(file.path(out, "panel", "panel_hourly_2025.csv"))
alarm <- fread(file.path(out, "main", "alarm_table.csv"))
strategy <- fread(file.path(out, "main", "strategy_summary.csv"))
strategy_ret <- fread(file.path(out, "main", "strategy_returns_long.csv"))
quality <- fread(file.path(out, "main", "alarm_quality_summary.csv"))
sparsity <- fread(file.path(out, "main", "alarm_sparsity_summary.csv"))
baseline <- fread(file.path(out, "diagnostics", "baseline_calibration_monitor.csv"))
placebo <- fread(file.path(out, "placebo", "placebo_monitoring_summary.csv"))
tc <- fread(file.path(out, "robustness", "transaction_cost_sensitivity.csv"))
rap <- fread(file.path(out, "robustness", "random_alarm_placebo.csv"))
mt <- fread(file.path(out, "robustness", "multiple_testing_adjustment.csv"))
boot <- fread(file.path(out, "robustness", "bootstrap_strategy_uncertainty.csv"))
cv10 <- fread(file.path(out, "cache", "critical_values_q10_m2160_h336.csv"))
meta <- fread(file.path(out, "run_metadata.csv"))

for (cc in c("hour_end")) panel[, (cc) := as.POSIXct(get(cc), tz = "UTC")]
for (cc in c("alarm_time", "trade_start_time", "train_start_time", "train_end_time",
             "monitor_start_time", "monitor_end_time")) {
  if (cc %in% names(alarm)) alarm[, (cc) := as.POSIXct(get(cc), tz = "UTC")]
}
if ("hour_end" %in% names(strategy_ret)) strategy_ret[, hour_end := as.POSIXct(hour_end, tz = "UTC")]

pct <- function(x, digits = 1) ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), 100 * x))
num <- function(x, digits = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), x))
intc <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
tex_escape <- function(x) {
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x
}
write_lines <- function(path, x) {
  writeLines(x, con = path, useBytes = TRUE)
}

pretty_method <- function(x) {
  out <- x
  out <- sub("^EProc_Adaptive$", "Adaptive e-process", out)
  out <- sub("^EProc_MultiStart$", "Multiple-start e-process", out)
  out <- sub("^EProc_Mix$", "Mixture e-process", out)
  out <- sub("^SSMS_KS_g0$", "SSMS-KS ($\\\\gamma=0$)", out)
  out <- sub("^SSMS_KS_g015$", "SSMS-KS ($\\\\gamma=0.15$)", out)
  out <- sub("^RSMS_KS_g0$", "RSMS-KS ($\\\\gamma=0$)", out)
  out <- sub("^RSMS_KS_g015$", "RSMS-KS ($\\\\gamma=0.15$)", out)
  out <- sub("^HAC_KS_g0$", "HAC-KS ($\\\\gamma=0$)", out)
  out <- sub("^HAC_KS_g015$", "HAC-KS ($\\\\gamma=0.15$)", out)
  out <- sub("_CvM_U$", "-CvM (Uniform)", out)
  out <- sub("_CvM_Early$", "-CvM (Early)", out)
  out <- sub("_CvM_Mid$", "-CvM (Mid)", out)
  out <- sub("_CvM_Late$", "-CvM (Late)", out)
  out
}
plot_method <- function(x) {
  out <- pretty_method(x)
  out <- gsub("\\$\\\\gamma=([0-9.]+)\\$", "gamma=\\1", out)
  out
}
method_family <- function(x) {
  fifelse(grepl("^EProc", x), "E-process",
  fifelse(grepl("^SSMS_KS", x), "SSMS-KS",
  fifelse(grepl("^RSMS_KS", x), "RSMS-KS",
  fifelse(grepl("^HAC_KS", x), "HAC-KS",
  fifelse(grepl("^SSMS_CvM", x), "SSMS-CvM",
  fifelse(grepl("^RSMS_CvM", x), "RSMS-CvM", "HAC-CvM"))))))
}
method_family_root <- function(x) {
  fifelse(grepl("^EProc", x), "EProc",
  fifelse(grepl("^SSMS", x), "SSMS",
  fifelse(grepl("^RSMS", x), "RSMS", "HAC")))
}
rule_class <- function(side, hold) {
  paste0(ifelse(side == "flat", "Flat", "Short"), " ", hold, "h")
}

perf <- function(r) {
  r <- as.numeric(r)
  gross <- cumprod(1 + fifelse(is.na(r), 0, r))
  total <- tail(gross, 1) - 1
  ann_factor <- 365 * 24
  mu <- mean(r, na.rm = TRUE) * ann_factor
  vol <- sd(r, na.rm = TRUE) * sqrt(ann_factor)
  downside <- sd(pmin(r, 0), na.rm = TRUE) * sqrt(ann_factor)
  dd <- gross / cummax(gross) - 1
  data.table(
    ann_return = mu,
    ann_vol = vol,
    sharpe = ifelse(vol > 0, mu / vol, NA_real_),
    sortino = ifelse(downside > 0, mu / downside, NA_real_),
    max_drawdown = min(dd, na.rm = TRUE),
    total_return = total,
    ending_gross = tail(gross, 1)
  )
}
bench <- perf(panel$asset_simple_ret)

strategy[, method_label := pretty_method(method)]
strategy[, detector_family := method_family(method)]
strategy[, rule_class := rule_class(strategy_side, hold_period)]
quality[, method_label := pretty_method(method)]
sparsity[, method_label := pretty_method(method)]
sparsity[, detector_family := method_family(method)]
alarm[, method_label := pretty_method(method)]
alarm[, detector_family := method_family(method)]
rap[, method_label := pretty_method(method)]
rap[, rule_class := rule_class(strategy_side, hold_period)]
tc[, rule_class := rule_class(strategy_side, hold_period)]

leader <- merge(strategy, quality, by = c("method", "strategy_side", "hold_period"), all.x = TRUE)
leader[, `:=`(
  method_label = pretty_method(method),
  detector_family = method_family(method),
  rule_class = rule_class(strategy_side, hold_period)
)]
leader[, precision_rank_value := fifelse(is.na(precision), -Inf, precision)]
leader[, capture_rank_value := fifelse(is.na(event_capture), -Inf, event_capture)]
leader[, `:=`(
  rank_sortino = frank(-sortino, ties.method = "average"),
  rank_drawdown = frank(-max_drawdown, ties.method = "average"),
  rank_precision = frank(-precision_rank_value, ties.method = "average"),
  rank_capture = frank(-capture_rank_value, ties.method = "average")
), by = rule_class]
leader[, practical_score := rowMeans(.SD), .SDcols = c("rank_sortino", "rank_drawdown", "rank_precision", "rank_capture")]

# Metadata macros
panel_start <- min(panel$hour_end)
panel_end <- max(panel$hour_end)
duration_months <- as.numeric(difftime(panel_end, panel_start, units = "days")) / 30.4375
sample_label <- sprintf("%s--%s", format(panel_start, "%Y-%m-%d"), format(panel_end, "%Y-%m-%d"))
meta_lines <- c(
  sprintf("%% Auto-generated empirical metadata from the %s empirical run", sample_label),
  sprintf("\\newcommand{\\EmpPanelRows}{%s}", intc(nrow(panel))),
  sprintf("\\newcommand{\\EmpPanelRowsPlain}{%d}", nrow(panel)),
  sprintf("\\newcommand{\\EmpPanelStart}{%s}", format(panel_start, "%Y-%m-%d %H:%M:%S")),
  sprintf("\\newcommand{\\EmpPanelEnd}{%s}", format(panel_end, "%Y-%m-%d %H:%M:%S")),
  sprintf("\\newcommand{\\EmpPanelStartDate}{%s}", format(panel_start, "%Y-%m-%d")),
  sprintf("\\newcommand{\\EmpPanelEndDate}{%s}", format(panel_end, "%Y-%m-%d")),
  sprintf("\\newcommand{\\EmpRollingWindows}{%s}", intc(meta$n_windows[1])),
  sprintf("\\newcommand{\\EmpPanelDuration}{%.1f}", duration_months),
  sprintf("\\newcommand{\\EmpCvSims}{%s}", intc(meta$n_cv_sims[1])),
  sprintf("\\newcommand{\\EmpRandomAlarmSims}{%s}", intc(meta$n_random_alarm[1])),
  sprintf("\\newcommand{\\EmpBootstrapSims}{%s}", intc(meta$n_bootstrap[1])),
  sprintf("\\newcommand{\\EmpBenchmarkTotal}{%s}", pct(bench$total_return)),
  sprintf("\\newcommand{\\EmpBenchmarkSortino}{%s}", num(bench$sortino, 3)),
  sprintf("\\newcommand{\\EmpRealityPAll}{%s}", num(mt[family == "ALL", reality_check_p_return], 3))
)
write_lines(file.path(tables_dir, "empirical_meta.tex"), meta_lines)

# Descriptive table
desc_vars <- list(
  "Panel A: response and spot-only baseline controls" = c(
    "Next-hour BTC return $Y_t$" = "y",
    "Lagged BTC return" = "ret_1",
    "24-hour return momentum" = "ret_mom",
    "24-hour realized volatility (ann.)" = "rv_ann_24",
    "Hourly high-low range" = "range_1",
    "Log spot quote volume" = "log_quote_volume"
  ),
  "Panel B: Deribit signals entering the omitted block $Z_{t-1}$" = c(
    "Call-minus-put imbalance" = "put_call_imbalance",
    "Log option activity" = "activity_log",
    "Option skew proxy" = "skew_proxy",
    "Term-structure slope" = "term_slope",
    "IV-RV spread" = "iv_rv_spread"
  ),
  "Panel C: auxiliary implied-volatility levels used to construct the signals" = c(
    "Near-dated ATM IV" = "near_atm_iv",
    "Mid-dated ATM IV" = "mid_atm_iv"
  )
)
desc_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Descriptive statistics for the BTC spot--Deribit option estimation panel}",
  "\\label{tab:empirical_desc}",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lrrrrrr}",
  "\\toprule",
  "Variable & $N$ & Mean & SD & P5 & Median & P95 \\\\",
  "\\midrule"
)
for (panel_name in names(desc_vars)) {
  desc_rows <- c(desc_rows, sprintf("\\multicolumn{7}{l}{\\textit{%s}} \\\\", panel_name))
  vars <- desc_vars[[panel_name]]
  for (lab in names(vars)) {
    v <- panel[[vars[[lab]]]]
    qs <- quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
    desc_rows <- c(desc_rows, sprintf(
      "%s & %d & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
      lab, sum(is.finite(v)), mean(v, na.rm = TRUE), sd(v, na.rm = TRUE),
      qs[1], qs[2], qs[3]
    ))
  }
  desc_rows <- c(desc_rows, "\\midrule")
}
desc_rows[length(desc_rows)] <- "\\bottomrule"
desc_rows <- c(desc_rows,
  "\\end{tabular}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.94\\textwidth}",
  "\\footnotesize\\textit{Notes:} The final complete-case estimation panel contains \\EmpPanelRows\\ hourly observations from \\EmpPanelStart\\ to \\EmpPanelEnd\\ UTC. The spot leg uses hourly BTCUSDT candles from the Binance Spot API. The option leg uses BTC option transactions from the Blockchain Research Center accessible-data service (``BTC Deribit Transactions / BTC Deribit Option Data''; \\url{https://blockchain-research-center.com/}), constructed from transaction-level records covering the empirical sample. Panel B lists the five variables included in the omitted block $Z_{t-1}$; Panel C lists auxiliary implied-volatility levels used to form the IV-RV spread and term-slope signals.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_desc.tex"), desc_rows)

# Rules table
rules_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Alarm-to-trade rules used in the empirical illustration}",
  "\\label{tab:empirical_rules}",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\begin{tabularx}{\\textwidth}{@{}l>{\\raggedright\\arraybackslash}p{3.5cm}>{\\raggedright\\arraybackslash}p{2.7cm}>{\\raggedright\\arraybackslash}X@{}}",
  "\\toprule",
  "Rule class & Exposure after an alarm & Hold length & Interpretation \\\\",
  "\\midrule",
  "Flat 6h & Move BTC exposure from $+1$ to $0$ & 6 observed hours & Tactical de-risking rule for brief downside-warning episodes. It asks whether the alarm is useful as a short-lived risk filter rather than as a directional short signal. \\\\",
  "Flat 24h & Move BTC exposure from $+1$ to $0$ & 24 observed hours & Persistent de-risking rule. It is designed for alarms that identify downside conditions lasting beyond a few hours. \\\\",
  "Short 6h & Move BTC exposure from $+1$ to $-1$ & 6 observed hours & Aggressive short-horizon overlay. It tests whether an alarm should be interpreted as an immediate bearish signal rather than merely as a cue to step aside. \\\\",
  "Short 24h & Move BTC exposure from $+1$ to $-1$ & 24 observed hours & Persistent bearish overlay. It is the most demanding rule because it requires alarm information to carry over into a longer short position. \\\\",
  "\\bottomrule",
  "\\end{tabularx}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.94\\textwidth}",
  "\\footnotesize\\textit{Notes:} Trading begins at the next hourly observation after the first threshold crossing within a rolling monitoring window. If multiple alarms overlap, the alarm state is extended rather than stacked. All returns and drawdown measures are computed on the complete-case hourly panel from \\EmpPanelStartDate\\ to \\EmpPanelEndDate.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_rules.tex"), rules_rows)

# Main method comparison table
wide <- dcast(strategy[, .(method, method_label, strategy_side, hold_period, total_return, sortino)],
              method + method_label ~ strategy_side + hold_period,
              value.var = c("total_return", "sortino"))
order_methods <- c("BENCHMARK", sort(unique(strategy$method)))
wide <- wide[match(sort(unique(strategy$method)), method)]
wide[, method_label := pretty_method(method)]
for (side in c("flat", "short")) for (h in c(6, 24)) {
  ctot <- paste0("total_return_", side, "_", h)
  csor <- paste0("sortino_", side, "_", h)
  mx <- max(wide[[ctot]], na.rm = TRUE)
  wide[, (paste0(ctot, "_fmt")) := fifelse(abs(get(ctot) - mx) < 1e-12,
                                           paste0("\\textbf{", pct(get(ctot)), "}"), pct(get(ctot)))]
  wide[, (paste0(csor, "_fmt")) := num(get(csor), 3)]
}
bench_row <- data.table(method_label = "Buy-and-hold benchmark")
for (side in c("flat", "short")) for (h in c(6, 24)) {
  bench_row[[paste0("total_return_", side, "_", h, "_fmt")]] <- pct(bench$total_return)
  bench_row[[paste0("sortino_", side, "_", h, "_fmt")]] <- num(bench$sortino, 3)
}
method_order <- c("EProc_Adaptive", "EProc_MultiStart", "EProc_Mix",
                  "SSMS_KS_g0", "SSMS_KS_g015", "RSMS_KS_g0", "RSMS_KS_g015", "HAC_KS_g0", "HAC_KS_g015",
                  "SSMS_CvM_U", "SSMS_CvM_Early", "SSMS_CvM_Mid", "SSMS_CvM_Late",
                  "RSMS_CvM_U", "RSMS_CvM_Early", "RSMS_CvM_Mid", "RSMS_CvM_Late",
                  "HAC_CvM_U", "HAC_CvM_Early", "HAC_CvM_Mid", "HAC_CvM_Late")
wide <- wide[match(method_order, method)]
method_table_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  sprintf("\\caption{Detector-by-detector performance across the four alarm-to-trade rules, %s}", sample_label),
  "\\label{tab:empirical_method_compare}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lrrrrrrrr}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Flat 6h} & \\multicolumn{2}{c}{Flat 24h} & \\multicolumn{2}{c}{Short 6h} & \\multicolumn{2}{c}{Short 24h} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Method & Total (\\%) & Sortino & Total (\\%) & Sortino & Total (\\%) & Sortino & Total (\\%) & Sortino \\\\",
  "\\midrule",
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
          bench_row$method_label,
          bench_row$total_return_flat_6_fmt, bench_row$sortino_flat_6_fmt,
          bench_row$total_return_flat_24_fmt, bench_row$sortino_flat_24_fmt,
          bench_row$total_return_short_6_fmt, bench_row$sortino_short_6_fmt,
          bench_row$total_return_short_24_fmt, bench_row$sortino_short_24_fmt),
  "\\midrule"
)
for (i in seq_len(nrow(wide))) {
  if (i %in% c(4, 10, 14, 18)) method_table_rows <- c(method_table_rows, "\\midrule")
  method_table_rows <- c(method_table_rows, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    wide$method_label[i],
    wide$total_return_flat_6_fmt[i], wide$sortino_flat_6_fmt[i],
    wide$total_return_flat_24_fmt[i], wide$sortino_flat_24_fmt[i],
    wide$total_return_short_6_fmt[i], wide$sortino_short_6_fmt[i],
    wide$total_return_short_24_fmt[i], wide$sortino_short_24_fmt[i]
  ))
}
method_table_rows <- c(method_table_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize\\textit{Notes:} Total return is the cumulative sample return earned by the corresponding alarm-to-trade rule over the complete-case hourly panel from \\EmpPanelStartDate\\ to \\EmpPanelEndDate. Sortino is computed from hourly strategy returns. Boldface marks the highest total return within a rule class. Abbreviations: U = Uniform, E = Early, M = Mid, L = Late, and the number in parentheses for KS rules is the boundary exponent $\\gamma$.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_method_compare.tex"), method_table_rows)

# Robustness table
best_rule <- leader[which.max(total_return)]
best_tc <- tc[method == best_rule$method & strategy_side == best_rule$strategy_side & hold_period == best_rule$hold_period]
tc_5 <- best_tc[cost_bps == 5 & short_funding_bps_per_hour == 0.25][1]
tc_10 <- best_tc[cost_bps == 10 & short_funding_bps_per_hour == 0.50][1]
mt_all <- mt[family == "ALL"][1]
mt_best_family <- mt[family == method_family_root(best_rule$method)][1]
rap_best <- rap[method == best_rule$method & strategy_side == best_rule$strategy_side & hold_period == best_rule$hold_period][1]
robust_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  sprintf("\\caption{Robustness diagnostics for the BTC--Deribit empirical illustration, %s}", sample_label),
  "\\label{tab:empirical_robustness}",
  "\\small",
  "\\setlength{\\tabcolsep}{5pt}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llr}",
  "\\toprule",
  "Diagnostic & Quantity & Value \\\\",
  "\\midrule",
  sprintf("Benchmark & Buy-and-hold total return (\\%%) & %s \\\\", pct(bench$total_return)),
  sprintf("Best raw overlay & %s, %s & %s \\\\", best_rule$method_label, best_rule$rule_class, pct(best_rule$total_return)),
  sprintf("Matched random alarms & p-value for best raw overlay & %.3f \\\\", rap_best$p_placebo_return),
  sprintf("Family-level multiple testing & Max-statistic p-value for %s detector family & %.3f \\\\", method_family_root(best_rule$method), mt_best_family$reality_check_p_return),
  sprintf("All-rule multiple testing & Max-statistic p-value over 84 detector-rule comparisons & %.3f \\\\", mt_all$reality_check_p_return),
  sprintf("Transaction costs & Best overlay, 5 bps turnover cost and 0.25 bps/hour short funding (\\%%) & %s \\\\", pct(tc_5$total_return)),
  sprintf("Transaction costs & Best overlay, 10 bps turnover cost and 0.50 bps/hour short funding (\\%%) & %s \\\\", pct(tc_10$total_return)),
  sprintf("Baseline-only monitor & Rejections among baseline-only detector-windows & %s of %s \\\\", intc(sum(baseline$reject, na.rm = TRUE)), intc(nrow(baseline))),
  sprintf("Placebo monitors & Alarm windows: lagged spot / block-permuted option / delayed option & %s / %s / %s \\\\",
          intc(placebo[placebo_type == "lagged_spot_only", sum(n_alarm_windows)]),
          intc(placebo[placebo_type == "block_permuted_week", sum(n_alarm_windows)]),
          intc(placebo[placebo_type == "delayed_option_168h", sum(n_alarm_windows)])),
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize\\textit{Notes:} The matched-random-alarm p-value compares the observed overlay with random alarm placements that match the same alarm count and holding period. The max-statistic multiple-testing p-values use shared random-alarm draws and the maximum excess-return statistic across detector-rule comparisons. The diagnostics are intended to qualify the economic interpretation of alarms rather than to validate a trading strategy.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_robustness.tex"), robust_rows)

# Rank-based summary tables
best_by_rule <- leader[order(practical_score), .SD[1], by = rule_class]
best_by_rule <- best_by_rule[match(c("Flat 6h", "Flat 24h", "Short 6h", "Short 24h"), rule_class)]
perf_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Best configuration under the rank-based summary index within each alarm-to-trade class}",
  "\\label{tab:empirical_perfclass}",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrrrrr}",
  "\\toprule",
  "Rule class & Best-ranked method & Ann.\\ return (\\%) & Total (\\%) & Sortino & Max DD (\\%) & Alarm hit (\\%) & Tail coverage (\\%) & Alarms \\\\",
  "\\midrule",
  sprintf("Buy-and-hold benchmark & No-alarm / always long & %s & %s & %s & %s & -- & -- & 0 \\\\",
          pct(bench$ann_return), pct(bench$total_return), num(bench$sortino, 3), pct(bench$max_drawdown))
)
for (i in seq_len(nrow(best_by_rule))) {
  r <- best_by_rule[i]
  perf_rows <- c(perf_rows, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %d \\\\",
    r$rule_class, r$method_label, pct(r$ann_return), pct(r$total_return), num(r$sortino, 3),
    pct(r$max_drawdown), pct(r$precision), pct(r$event_capture), r$n_alarm_entries
  ))
}
perf_rows <- c(perf_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize\\textit{Notes:} The rank-based summary index is the equally weighted average of within-class ranks for Sortino, maximum drawdown, alarm hit rate, and tail-event coverage. Alarm hit rate is the fraction of alarms whose holding interval contains at least one lower-tail BTC return. Tail-event coverage is the fraction of sample lower-tail return hours covered by alarm-induced holding intervals. Lower values of the rank index indicate a better within-rule ranking. Total and annualized returns are reported as descriptive sample quantities.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_perfclass.tex"), perf_rows)

family <- leader[, .(
  ann_return = mean(ann_return, na.rm = TRUE),
  total_return = mean(total_return, na.rm = TRUE),
  sortino = mean(sortino, na.rm = TRUE),
  max_drawdown = mean(max_drawdown, na.rm = TRUE),
  precision = mean(precision, na.rm = TRUE),
  event_capture = mean(event_capture, na.rm = TRUE),
  avg_alarms = mean(n_alarm_entries, na.rm = TRUE)
), by = .(rule_class, detector_family)]
family_order <- c("E-process", "SSMS-KS", "RSMS-KS", "HAC-KS", "SSMS-CvM", "RSMS-CvM", "HAC-CvM")
family <- family[order(match(rule_class, c("Flat 6h", "Flat 24h", "Short 6h", "Short 24h")),
                       match(detector_family, family_order))]
fam_rows <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Average performance by detector family and trading rule}",
  "\\label{tab:empirical_family}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{llrrrrrrr}",
  "\\toprule",
  "Rule class & Family & Ann.\\ return (\\%) & Total (\\%) & Sortino & Max DD (\\%) & Alarm hit (\\%) & Tail coverage (\\%) & Avg. alarms \\\\",
  "\\midrule"
)
last_rule <- ""
for (i in seq_len(nrow(family))) {
  r <- family[i]
  if (nzchar(last_rule) && r$rule_class != last_rule) fam_rows <- c(fam_rows, "\\midrule")
  fam_rows <- c(fam_rows, sprintf(
    "%s & %s & %s & %s & %s & %s & %s & %s & %.1f \\\\",
    ifelse(r$rule_class == last_rule, "", r$rule_class), r$detector_family,
    pct(r$ann_return), pct(r$total_return), num(r$sortino, 3), pct(r$max_drawdown),
    pct(r$precision), pct(r$event_capture), r$avg_alarms
  ))
  last_rule <- r$rule_class
}
fam_rows <- c(fam_rows,
  "\\bottomrule",
  "\\end{tabular}%",
  "}",
  "\\vspace{0.2cm}",
  "\\begin{minipage}{0.96\\textwidth}",
  "\\footnotesize\\textit{Notes:} Each entry averages over all methods inside the named detector family and rule class. Alarm hit rate and tail-event coverage summarize whether alarm windows overlap with lower-tail realized BTC-return episodes; they are descriptive diagnostics rather than formal tests.",
  "\\end{minipage}",
  "\\end{table}"
)
write_lines(file.path(tables_dir, "tab_empirical_family.tex"), fam_rows)

full_table <- function(side, path, label, caption) {
  x <- leader[strategy_side == side][order(hold_period, practical_score)]
  rows <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{llrrrrrrrrr}",
    "\\toprule",
    "Hold & Method & Ann.\\ return (\\%) & Total (\\%) & Ending gross & Sortino & Max DD (\\%) & Alarm hit (\\%) & Tail cov. (\\%) & Alarms & Composite rank \\\\",
    "\\midrule"
  )
  for (h in c(6L, 24L)) {
    if (h == 24L) rows <- c(rows, "\\midrule")
    rows <- c(rows, sprintf("\\multicolumn{11}{l}{\\textit{Panel %s: %d-hour holding period}} \\\\", ifelse(h == 6L, "A", "B"), h))
    subx <- x[hold_period == h]
    for (i in seq_len(nrow(subx))) {
      r <- subx[i]
      rows <- c(rows, sprintf(
        "%dh & %s & %s & %s & %.3f & %s & %s & %s & %s & %d & %.3f \\\\",
        h, r$method_label, pct(r$ann_return), pct(r$total_return), r$ending_gross,
        num(r$sortino, 3), pct(r$max_drawdown), pct(r$precision), pct(r$event_capture),
        r$n_alarm_entries, r$practical_score
      ))
    }
  }
  rows <- c(rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.15cm}",
    "\\begin{minipage}{0.96\\textwidth}",
    "\\footnotesize\\textit{Notes:} Composite rank denotes the rank-based summary index, defined as the equally weighted average of within-class ranks for Sortino, maximum drawdown, alarm hit rate, and tail-event coverage. Alarm hit rate is the fraction of alarms whose holding interval contains at least one lower-tail BTC return. Tail-event coverage is the fraction of sample lower-tail return hours covered by alarm-induced holding intervals. Rows are sorted by this index within each holding-period panel.",
    "\\end{minipage}",
    "\\end{table}"
  )
  write_lines(file.path(tables_dir, path), rows)
}
full_table("flat", "tab_empirical_flat_full_panel.tex", "tab:empirical_flat_full_panel",
           "Full empirical results: flat-on-alarm rules")
full_table("short", "tab_empirical_short_full_panel.tex", "tab:empirical_short_full_panel",
           "Full empirical results: short-on-alarm rules")

# Also refresh legacy single-rule table files for consistency.
for (side in c("flat", "short")) for (h in c(6L, 24L)) {
  file <- sprintf("tab_empirical_%s%d_full.tex", side, h)
  subx <- leader[strategy_side == side & hold_period == h][order(practical_score)]
  rows <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{Full empirical results: %s-on-alarm, %d-hour hold}", side, h),
    sprintf("\\label{tab:empirical_%s%d_full}", side, h),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{p{3.0cm}rrrrrrrrr}",
    "\\toprule",
    "Method & Ann.\\ return (\\%) & Total (\\%) & Ending gross & Sortino & Max DD (\\%) & Alarm hit (\\%) & Tail cov. (\\%) & Alarms & Composite rank \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(subx))) {
    r <- subx[i]
    rows <- c(rows, sprintf(
      "%s & %s & %s & %.3f & %s & %s & %s & %s & %d & %.3f \\\\",
      r$method_label, pct(r$ann_return), pct(r$total_return), r$ending_gross,
      num(r$sortino, 3), pct(r$max_drawdown), pct(r$precision), pct(r$event_capture),
      r$n_alarm_entries, r$practical_score
    ))
  }
  rows <- c(rows, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}")
  write_lines(file.path(tables_dir, file), rows)
}

# Figures
theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = base_size + 1),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom",
          axis.title = element_text(size = base_size))
}
family_cols <- c("E-process" = "#6B7280", "SSMS" = "#2B6CB0", "RSMS" = "#D97706", "HAC" = "#2F855A")
family_shapes <- c("E-process" = 16, "SSMS" = 17, "RSMS" = 15, "HAC" = 3)
method_cols <- c(
  "EProc_Adaptive" = "#4B5563", "EProc_MultiStart" = "#6B7280", "EProc_Mix" = "#9CA3AF",
  "SSMS_KS_g0" = "#1D4E89", "SSMS_KS_g015" = "#2B6CB0",
  "SSMS_CvM_U" = "#63B3ED", "SSMS_CvM_Early" = "#4299E1",
  "SSMS_CvM_Mid" = "#3182CE", "SSMS_CvM_Late" = "#2C5282",
  "RSMS_KS_g0" = "#B45309", "RSMS_KS_g015" = "#D97706",
  "RSMS_CvM_U" = "#F6AD55", "RSMS_CvM_Early" = "#ED8936",
  "RSMS_CvM_Mid" = "#DD6B20", "RSMS_CvM_Late" = "#9C4221",
  "HAC_KS_g0" = "#276749", "HAC_KS_g015" = "#2F855A",
  "HAC_CvM_U" = "#68D391", "HAC_CvM_Early" = "#48BB78",
  "HAC_CvM_Mid" = "#38A169", "HAC_CvM_Late" = "#22543D"
)

feat_map <- data.table(
  variable = c("ret_mom", "rv_ann_24", "close", "iv_rv_spread",
               "activity_log", "skew_proxy", "put_call_imbalance", "term_slope"),
  label = c("24-hour momentum", "24-hour realized volatility", "BTC spot price",
            "IV-RV spread", "Log option activity", "Skew proxy",
            "Call-minus-put imbalance", "Term slope")
)
feat_long <- melt(panel[, c("hour_end", feat_map$variable), with = FALSE], id.vars = "hour_end",
                  variable.name = "variable", value.name = "value")
feat_long <- merge(feat_long, feat_map, by = "variable", all.x = TRUE)
p_feat <- ggplot(feat_long, aes(hour_end, value)) +
  geom_line(color = "black", linewidth = 0.25, na.rm = TRUE) +
  facet_wrap(~label, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = NULL) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  theme_paper(9)
ggsave(file.path(fig_dir, "fig_empirical_feature_panels.png"), p_feat, width = 10, height = 7.5, dpi = 300)

alarm_plot <- alarm[!is.na(alarm_time)]
alarm_plot[, detector_group := fifelse(grepl("^EProc", method), "E-process",
                                fifelse(grepl("^SSMS", method), "SSMS",
                                fifelse(grepl("^RSMS", method), "RSMS", "HAC")))]
alarm_plot[, plot_label := plot_method(method)]
method_levels <- alarm_plot[, .N, by = plot_label][order(N, plot_label), plot_label]
alarm_plot[, plot_label := factor(plot_label, levels = method_levels)]
price_plot <- ggplot(panel, aes(hour_end, close)) +
  geom_line(color = "black", linewidth = 0.28) +
  scale_y_continuous(labels = label_dollar()) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(x = NULL, y = "BTC price") +
  theme_paper(9) +
  theme(legend.position = "none")
raster <- ggplot(alarm_plot, aes(alarm_time, plot_label, color = detector_group, shape = detector_group)) +
  geom_point(size = 0.95, alpha = 0.9, stroke = 0.25) +
  scale_color_manual(values = family_cols) +
  scale_shape_manual(values = family_shapes) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(x = NULL, y = NULL, color = NULL, shape = NULL) +
  theme_paper(8)
ggsave(file.path(fig_dir, "fig_empirical_alarm_raster.png"), price_plot / raster + plot_layout(heights = c(1, 2.6)),
       width = 9.5, height = 7.6, dpi = 300)

heat <- alarm[, .(
  alarm = as.integer(any(!is.na(alarm_abs_idx))),
  alarm_frac = if (any(!is.na(alarm_k))) min(alarm_k, na.rm = TRUE) / max(monitor_end_idx - monitor_start_idx + 1L, na.rm = TRUE) else NA_real_
), by = .(method, window_id)]
all_methods <- unique(alarm[, .(method)])
all_windows <- data.table(window_id = sort(unique(alarm$window_id)))
heat <- merge(CJ(method = all_methods$method, window_id = all_windows$window_id),
              all_methods, by = "method", all.x = TRUE)
heat <- merge(heat, alarm[, .(
  alarm = as.integer(any(!is.na(alarm_abs_idx))),
  alarm_frac = if (any(!is.na(alarm_k))) min(alarm_k, na.rm = TRUE) / max(monitor_end_idx - monitor_start_idx + 1L, na.rm = TRUE) else NA_real_
), by = .(method, window_id)],
              by = c("method", "window_id"), all.x = TRUE)
heat[is.na(alarm), alarm := 0L]
heat[alarm == 0L, alarm_frac := NA_real_]
heat[, detector_group := fifelse(grepl("^EProc", method), "E-process",
                          fifelse(grepl("^SSMS", method), "SSMS",
                          fifelse(grepl("^RSMS", method), "RSMS", "HAC")))]
heat[, plot_label := plot_method(method)]
heat[, plot_label := factor(plot_label, levels = rev(plot_method(unique(method))))]
p_heat <- ggplot(heat, aes(window_id, plot_label, fill = alarm_frac)) +
  geom_tile(color = "grey85", linewidth = 0.05) +
  scale_fill_gradientn(colors = c("white", "#C6DBEF", "#6BAED6", "#08519C"),
                       na.value = "grey94", limits = c(0, 1),
                       breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                       labels = c("0.00", "0.25", "0.50", "0.75", "1.00")) +
  labs(x = "Rolling-window index", y = NULL, fill = "Alarm time\nwithin window\n(light=early)") +
  theme_paper(8) +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "fig_empirical_alarm_heatmap.png"), p_heat, width = 9.5, height = 6.5, dpi = 300)

flat24 <- leader[strategy_side == "flat" & hold_period == 24][order(practical_score)][1:12]
flat24[, label := sprintf("Ret %s%% | S %s | DD %s%% | Hit %s%% | Cov %s%%",
                          pct(total_return), num(sortino, 2), pct(max_drawdown),
                          pct(precision), pct(event_capture))]
flat24[, plot_label := plot_method(method)]
flat24[, plot_label := factor(plot_label, levels = rev(plot_label))]
flat24[, detector_group := fifelse(grepl("^EProc", method), "E-process",
                            fifelse(grepl("^SSMS", method), "SSMS",
                            fifelse(grepl("^RSMS", method), "RSMS", "HAC")))]
p_top <- ggplot(flat24, aes(practical_score, plot_label, fill = detector_group)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = label), hjust = -0.02, size = 2.25) +
  scale_fill_manual(values = family_cols) +
  scale_x_continuous(limits = c(0, max(flat24$practical_score, na.rm = TRUE) + 8),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Composite rank (lower is better)", y = NULL, fill = NULL) +
  coord_cartesian(clip = "off") +
  theme_paper(9) +
  theme(legend.position = "none",
        plot.margin = margin(5.5, 30, 5.5, 5.5))
ggsave(file.path(fig_dir, "fig_empirical_top_methods_flat24.png"), p_top, width = 10.5, height = 5.1, dpi = 300)

ret_flat24 <- strategy_ret[strategy_side == "flat" & hold_period == 24]
ret_flat24[, plot_label := plot_method(method)]
ret_flat24[, `:=`(
  strategy_equity = cumprod(1 + fifelse(is.na(strategy_return), 0, strategy_return)),
  benchmark_equity = cumprod(1 + fifelse(is.na(benchmark_return), 0, benchmark_return))
), by = method]
ret_flat24[, plot_label := factor(plot_label, levels = plot_method(method_order))]
ret_equity_long <- melt(
  ret_flat24,
  id.vars = c("hour_end", "method", "plot_label"),
  measure.vars = c("benchmark_equity", "strategy_equity"),
  variable.name = "series", value.name = "equity"
)
ret_equity_long[, series := fifelse(series == "benchmark_equity", "Buy-and-hold", "Alarm rule")]
ret_equity_long[, color_key := fifelse(series == "Buy-and-hold", "Buy-and-hold", method)]
p_equity <- ggplot(ret_equity_long, aes(hour_end, equity)) +
  geom_line(aes(color = color_key, linetype = series), linewidth = 0.32) +
  scale_color_manual(values = c("Buy-and-hold" = "black", method_cols)) +
  scale_linetype_manual(values = c("Buy-and-hold" = "solid", "Alarm rule" = "longdash")) +
  facet_wrap(~plot_label, ncol = 3) +
  scale_x_datetime(date_breaks = "4 months", date_labels = "%b") +
  labs(x = NULL, y = "Gross return", linetype = NULL) +
  theme_paper(7) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  guides(color = "none")
ggsave(file.path(fig_dir, "fig_empirical_strategy_equity.png"), p_equity, width = 9.5, height = 12.2, dpi = 300)

selected_equity <- strategy[order(-total_return)][seq_len(min(4L, .N)),
  .(method, strategy_side, hold_period)]
selected_equity[, selected_label := sprintf("%s, %s %dh",
                                            plot_method(method),
                                            fifelse(strategy_side == "short", "Short", "Flat"),
                                            hold_period)]
ret_selected <- merge(strategy_ret, selected_equity,
                      by = c("method", "strategy_side", "hold_period"),
                      all = FALSE)
ret_selected[, `:=`(
  strategy_equity = cumprod(1 + fifelse(is.na(strategy_return), 0, strategy_return)),
  benchmark_equity = cumprod(1 + fifelse(is.na(benchmark_return), 0, benchmark_return))
), by = .(method, strategy_side, hold_period)]
ret_selected_long <- melt(
  ret_selected,
  id.vars = c("hour_end", "method", "selected_label"),
  measure.vars = c("benchmark_equity", "strategy_equity"),
  variable.name = "series", value.name = "equity"
)
ret_selected_long[, series := fifelse(series == "benchmark_equity", "Buy-and-hold", "Alarm overlay")]
ret_selected_long[, selected_label := factor(selected_label, levels = selected_equity$selected_label)]
p_equity_selected <- ggplot(ret_selected_long, aes(hour_end, equity)) +
  geom_line(aes(color = series, linetype = series), linewidth = 0.42) +
  scale_color_manual(values = c("Buy-and-hold" = "black", "Alarm overlay" = "#1F78B4")) +
  scale_linetype_manual(values = c("Buy-and-hold" = "solid", "Alarm overlay" = "longdash")) +
  facet_wrap(~selected_label, ncol = 2) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(x = NULL, y = "Gross return", color = NULL, linetype = NULL) +
  theme_paper(9) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_empirical_equity_selected.png"), p_equity_selected,
       width = 10.5, height = 6.0, dpi = 300)

message("Wrote empirical tables and figures.")

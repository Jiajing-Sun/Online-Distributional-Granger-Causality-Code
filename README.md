# Online Monitoring of Distributional Granger Causality: Replication Code

This repository contains R code for the paper:

> Jiajing Sun, Abderrahim Taamouti, and Yongmiao Hong, "Online Monitoring of Distributional Granger Causality".

The repository is intended for private replication and review of the numerical work. It contains code only; raw Deribit option-transaction data and generated outputs are not included.

## Repository Contents

```text
simulation/
  R/                           Core simulation, DGP, monitoring, critical-value, and e-process functions
  run_*.R                      Scripts used to produce simulation tables and figures
  openend_calibration_extras/  Open-end monitoring critical-value utilities

empirical_2025_pipeline/
  R/                           Empirical panel, monitoring, and robustness functions
  run_empirical_2025.R         Main script for the BTC--Deribit application
  make_empirical_2025_latex_figures.R
                               Builds LaTeX tables and figures from empirical outputs
  glmy_raw_ingest_reference/   Minimal copied parser reference for raw Deribit/Binance ingestion
```

## Data Availability

The empirical application uses:

- hourly BTCUSDT spot candles from Binance, fetched by the code; and
- BTC option transaction records from Deribit, accessed through the Blockchain Research Center data service.

The Deribit transaction files are proprietary third-party data and are not redistributed here. The empirical code can read either a local raw-data directory or explicit extra CSV file paths.

Default empirical sample:

```text
2025-04-01 00:00:00 UTC to 2026-03-31 23:00:00 UTC
```

The current local setup used an almost complete twelve-month BTC--Deribit panel with small remaining missing segments in the raw transaction coverage.

## R Requirements

The scripts use base R plus the following packages:

```r
install.packages(c(
  "data.table",
  "jsonlite",
  "httr",
  "quantreg",
  "ggplot2",
  "scales",
  "patchwork"
))
```

`rstudioapi` is optional and is used only to infer script locations when running interactively in RStudio.

## Empirical Code

Run from the repository root:

```bash
Rscript empirical_2025_pipeline/run_empirical_2025.R \
  --deribit_dir="/path/to/BTC Deribit Transactions" \
  --sample_start="2025-04-01 00:00:00" \
  --sample_end="2026-03-31 23:00:00" \
  --n_cv_sims=20000 \
  --ncores=4
```

If some raw Deribit files are outside the main directory, pass them explicitly:

```bash
Rscript empirical_2025_pipeline/run_empirical_2025.R \
  --deribit_dir="/path/to/BTC Deribit Transactions" \
  --extra_deribit_files="/path/to/file1.csv,/path/to/file2.csv" \
  --sample_start="2025-04-01 00:00:00" \
  --sample_end="2026-03-31 23:00:00" \
  --n_cv_sims=20000 \
  --ncores=4
```

Alternatively, set `DERIBIT_DIR` in the shell instead of passing `--deribit_dir`.

The default output directory is:

```text
empirical_2025_pipeline/output_apr2025_mar2026_ncv20000/
```

For a lightweight check that does not require raw data, use synthetic smoke mode:

```bash
Rscript empirical_2025_pipeline/run_empirical_2025.R \
  --smoke=true \
  --smoke_synthetic=true \
  --n_cv_sims=25 \
  --n_random_alarm=50 \
  --n_bootstrap=50 \
  --ncores=2 \
  --force=true
```

## Simulation Code

The simulation scripts are in `simulation/`. The main entry points are:

```bash
Rscript simulation/run_all_paper_sims.R
Rscript simulation/run_table1_null.R
Rscript simulation/run_table2_abrupt.R
Rscript simulation/run_table3_tail_eprocess.R
Rscript simulation/run_table4_gradual.R
Rscript simulation/run_appendix_training_contam.R
```

The finite-sample critical-value assessment is run with:

```bash
Rscript simulation/run_joe_redesign_sims.R
```

## Reproducibility Notes

- All scripts are written to avoid storing raw transaction data in the repository.
- Empirical results can vary if the local Deribit raw files differ from the data deliveries used for the manuscript.
- The empirical critical values for the BTC--Deribit horizon are simulated because the empirical monitoring horizon is outside the main tabulated asymptotic grid.
- Long simulations should be run with an explicit `--ncores` or script-level core setting so that the machine remains usable.

## License / Sharing

This repository is currently for private replication and review. Please do not redistribute proprietary raw Deribit transaction data.

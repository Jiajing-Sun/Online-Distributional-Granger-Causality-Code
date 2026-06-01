# Empirical BTC--Deribit Code Notes

This directory contains the code for the BTC--Deribit empirical illustration reported in the manuscript. The default sample is April 1, 2025 through March 31, 2026, with the final predictor hour trimmed so the return outcome remains inside the sample.

The intended Deribit input files are the raw CSV deliveries covering this sample. Depending on the local download layout, these can be supplied either through `--deribit_dir=<folder>` or explicitly through `--extra_deribit_files=<comma-separated files>`. The code discovers files named `YYYYMMDD-YYYYMMDD.csv` that overlap the requested sample window.

If `--deribit_dir` is omitted, the code reads the `DERIBIT_DIR` environment variable and otherwise falls back to `data/raw_deribit`.

```text
/path/to/BTC Deribit Transactions/20250301-20250430.csv
/path/to/BTC Deribit Transactions/20250501-20250630.csv
/path/to/BTC Deribit Transactions/20250701-20250831.csv
/path/to/BTC Deribit Transactions/20250901-20251031.csv
/path/to/BTC Deribit Transactions/20251101-20251231.csv
... plus any supplementary January--March 2026 delivery files passed through `--extra_deribit_files`
```

The manuscript uses one clean twelve-month empirical sample, not the earlier full-calendar-year 2025 pilot. The default Monte Carlo critical-value count is `n_cv_sims = 20000`, matching the empirical tables and figures currently reported in the paper.

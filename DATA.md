# Data Notes

This repository does not include raw Deribit option-transaction data.

The empirical BTC--Deribit application was built from:

- Binance BTCUSDT hourly spot candles, which the code can fetch from Binance; and
- BTC option transaction records from Deribit, accessed through the Blockchain Research Center data service.

The Deribit data are third-party/proprietary and must be obtained separately. To run the empirical code, provide the local Deribit data location through:

```bash
--deribit_dir="/path/to/BTC Deribit Transactions"
```

or pass specific files through:

```bash
--extra_deribit_files="/path/to/file1.csv,/path/to/file2.csv"
```

The manuscript sample is April 1, 2025 through March 31, 2026. The final local run used an almost complete one-year panel with small remaining missing segments in the raw Deribit coverage.

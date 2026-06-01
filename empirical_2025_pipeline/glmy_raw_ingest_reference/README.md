# GLMY Raw-Ingest Reference

This folder contains copied reference scripts from the earlier GLMY project. The original GLMY directories were not modified.

The scripts were copied from a separate GLMY project directory for parser reference only. The original GLMY project directory is not part of this repository.

Copied files:

```text
R_original/00_packages.R
R_original/01_config.R
R_original/parse_deribit_trades.R
R_original/10_ingest_options_deribit.R
R_original/binance_spot.R
```

These files are kept as reference material for reading raw Deribit and Binance spot data. They should not be sourced directly into the JoE empirical pipeline without adaptation, because the GLMY scripts contain project-specific paths and node-construction choices. The adapted JoE pipeline currently uses the April 2025--March 2026 BTC--Deribit sample.

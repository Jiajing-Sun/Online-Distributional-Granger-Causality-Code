# ==============================================================
# simulation_runner.R -- Monte Carlo driver
# ==============================================================

source_project("R", "utils.R")
source_project("R", "critical_values.R")
source_project("R", "one_replication.R")

summarise_one_config <- function(df_one) {
  keep_cols <- intersect(
    c("design", "model_type", "instrument_type", "m", "T", "q", "c_val", "break_frac", "k_star"),
    names(df_one)
  )
  out <- df_one[1, keep_cols, drop = FALSE]

  rej_cols <- grep("_rej$", names(df_one), value = TRUE)
  methods <- sub("_rej$", "", rej_cols)

  for (mth in methods) {
    rej <- df_one[[paste0(mth, "_rej")]]
    fa <- if (paste0(mth, "_fa") %in% names(df_one)) df_one[[paste0(mth, "_fa")]] else NA
    det <- if (paste0(mth, "_det") %in% names(df_one)) df_one[[paste0(mth, "_det")]] else NA
    delay <- if (paste0(mth, "_delay") %in% names(df_one)) df_one[[paste0(mth, "_delay")]] else NA
    arl0 <- if (paste0(mth, "_arl0") %in% names(df_one)) df_one[[paste0(mth, "_arl0")]] else NA

    out[[paste0(mth, "_rej_rate")]] <- mean(rej, na.rm = TRUE)
    out[[paste0(mth, "_fa_rate")]] <- if (all(is.na(fa))) NA_real_ else mean(fa, na.rm = TRUE)
    out[[paste0(mth, "_det_rate")]] <- if (all(is.na(det))) NA_real_ else mean(det, na.rm = TRUE)
    out[[paste0(mth, "_avg_delay")]] <- if (all(is.na(delay))) NA_real_ else mean(delay, na.rm = TRUE)
    out[[paste0(mth, "_med_delay")]] <- if (all(is.na(delay))) NA_real_ else stats::median(delay, na.rm = TRUE)
    out[[paste0(mth, "_arl0_mean")]] <- if (all(is.na(arl0))) NA_real_ else mean(arl0, na.rm = TRUE)
  }

  out
}

run_simulation_grid <- function(config_grid,
                                model_type = c("quantile", "expectile"),
                                tau_grid,
                                instrument_type = "z",
                                tau_weight_scheme = "equal",
                                include_eprocess = FALSE,
                                eprocess_tau_grid = NULL,
                                eprocess_feature = "z",
                                nrep = 200,
                                ncores = 1L,
                                alpha = 0.05,
                                gamma_vec = c(0, 0.15),
                                cvm_weights = c("U", "Late", "Early", "Mid"),
                                prefer_quantreg = TRUE,
                                seed_base = 20260310,
                                output_dir = file.path("output", "default_run"),
                                file_stub = "sim") {
  model_type <- match.arg(model_type)
  project_root <- get_project_root()
  if (!grepl("^(/|[A-Za-z]:[/\\])", output_dir)) output_dir <- project_path(output_dir)
  dir_create(output_dir)

  cv <- load_critical_values(
    path_base = project_path("critical_values", "critical_values_all.csv"),
    path_weights = project_path("critical_values", "critical_values_all_weights.csv")
  )

  config_grid <- as.data.frame(config_grid, stringsAsFactors = FALSE)
  all_summary <- vector("list", nrow(config_grid))

  cl <- make_cluster(ncores)
  on.exit(stop_cluster(cl), add = TRUE)

  if (!is.null(cl)) {
    parallel::clusterExport(
      cl,
      varlist = c("project_root", "run_one_replication"),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      setwd(project_root)
      options(distgc.project_root = project_root)
      source(file.path(project_root, "R", "utils.R"))
      source(file.path(project_root, "R", "critical_values.R"))
      source(file.path(project_root, "R", "weights.R"))
      source(file.path(project_root, "R", "monitors.R"))
      source(file.path(project_root, "R", "fit_distributional_models.R"))
      source(file.path(project_root, "R", "dgp_online_distgc.R"))
      source(file.path(project_root, "R", "eprocess.R"))
      source(file.path(project_root, "R", "one_replication.R"))
      NULL
    })
  }

  for (cc in seq_len(nrow(config_grid))) {
    conf <- config_grid[cc, , drop = FALSE]

    inst_this <- if ("instrument_type" %in% names(conf) && !is.na(conf$instrument_type)) conf$instrument_type else instrument_type
    include_ep_this <- if ("include_eprocess" %in% names(conf) && !is.na(conf$include_eprocess)) isTRUE(conf$include_eprocess) else include_eprocess
    ep_feature_this <- if ("eprocess_feature" %in% names(conf) && !is.na(conf$eprocess_feature)) conf$eprocess_feature else eprocess_feature

    message("--------------------------------------------------")
    message(sprintf(
      "Config %s / %s: design=%s, model=%s, m=%s, T=%s, c=%s, break=%s, instrument=%s",
      cc, nrow(config_grid),
      conf$design, model_type, conf$m, conf$T, conf$c_val, conf$break_frac, inst_this
    ))

    rep_ids <- seq_len(nrep)

    if (!is.null(cl)) {
      parallel::clusterExport(
        cl,
        varlist = c("cv", "conf", "model_type", "tau_grid", "inst_this", "tau_weight_scheme",
                    "alpha", "gamma_vec", "cvm_weights", "include_ep_this", "eprocess_tau_grid",
                    "ep_feature_this", "prefer_quantreg", "seed_base", "cc"),
        envir = environment()
      )

      df_list <- parallel::parLapply(cl, rep_ids, function(rid) {
        run_one_replication(
          rep_id = rid,
          cv = cv,
          design = conf$design,
          m = conf$m,
          T = conf$T,
          c_val = conf$c_val,
          break_frac = conf$break_frac,
          model_type = model_type,
          tau_grid = tau_grid,
          instrument_type = inst_this,
          tau_weight_scheme = tau_weight_scheme,
          alpha = alpha,
          gamma_vec = gamma_vec,
          cvm_weights = cvm_weights,
          include_eprocess = include_ep_this,
          eprocess_tau_grid = eprocess_tau_grid,
          eprocess_feature = ep_feature_this,
          prefer_quantreg = prefer_quantreg,
          seed_base = seed_base + 100000L * cc
        )
      })
    } else {
      df_list <- lapply(rep_ids, function(rid) {
        run_one_replication(
          rep_id = rid,
          cv = cv,
          design = conf$design,
          m = conf$m,
          T = conf$T,
          c_val = conf$c_val,
          break_frac = conf$break_frac,
          model_type = model_type,
          tau_grid = tau_grid,
          instrument_type = inst_this,
          tau_weight_scheme = tau_weight_scheme,
          alpha = alpha,
          gamma_vec = gamma_vec,
          cvm_weights = cvm_weights,
          include_eprocess = include_ep_this,
          eprocess_tau_grid = eprocess_tau_grid,
          eprocess_feature = ep_feature_this,
          prefer_quantreg = prefer_quantreg,
          seed_base = seed_base + 100000L * cc
        )
      })
    }

    raw_df <- do.call(rbind, df_list)
    raw_path <- file.path(output_dir, sprintf("%s_raw_%02d.csv", file_stub, cc))
    write.csv(raw_df, raw_path, row.names = FALSE)

    all_summary[[cc]] <- summarise_one_config(raw_df)
  }

  summary_df <- do.call(rbind, all_summary)
  sum_path <- file.path(output_dir, sprintf("%s_summary.csv", file_stub))
  write.csv(summary_df, sum_path, row.names = FALSE)
  invisible(summary_df)
}

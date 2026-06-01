
# ==============================================================
# openend_critical_values.R -- Monte Carlo open-end critical values
# ==============================================================

source(file.path("R", "utils_openend.R"))
source(file.path("R", "openend_limit_statistics.R"))

simulate_openend_critical_values <- function(q_max = 10L,
                                             gamma_vec = c(0, 0.15),
                                             weight_names = c("U", "Early", "Late", "Mid"),
                                             alpha_levels = c(0.05, 0.10),
                                             nrep = 5000L,
                                             n_train_grid = 1500L,
                                             n_open_grid = 2000L,
                                             ridge = 1e-10,
                                             range_floor = 1e-8,
                                             ncores = 1L,
                                             seed = 13579L,
                                             verbose = TRUE) {
  meta <- build_openend_meta(q_max = q_max,
                             gamma_vec = gamma_vec,
                             weight_names = weight_names)

  seeds <- seed + seq_len(nrep)

  if (verbose) {
    message(sprintf("Simulating open-end critical values: q_max=%d, nrep=%d, train-grid=%d, open-grid=%d, ncores=%d",
                    q_max, nrep, n_train_grid, n_open_grid, ncores))
  }

  ncores <- as.integer(max(1L, ncores))

  if (ncores <= 1L) {
    values_list <- lapply(seeds, function(sd) {
      set.seed(sd)
      simulate_one_openend_rep(q_max = q_max,
                               n_train_grid = n_train_grid,
                               n_open_grid = n_open_grid,
                               gamma_vec = gamma_vec,
                               weight_names = weight_names,
                               ridge = ridge,
                               range_floor = range_floor)
    })
  } else {
    cl <- make_cluster(ncores)
    on.exit(stop_cluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c("simulate_one_openend_rep",
                  "build_openend_meta",
                  "trapz_unit_interval",
                  "safe_solve",
                  "openend_weight_x",
                  "q_max",
                  "n_train_grid",
                  "n_open_grid",
                  "gamma_vec",
                  "weight_names",
                  "ridge",
                  "range_floor"),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, NULL)

    values_list <- parallel::parLapply(cl, seeds, function(sd) {
      set.seed(sd)
      simulate_one_openend_rep(q_max = q_max,
                               n_train_grid = n_train_grid,
                               n_open_grid = n_open_grid,
                               gamma_vec = gamma_vec,
                               weight_names = weight_names,
                               ridge = ridge,
                               range_floor = range_floor)
    })
  }

  values_mat <- do.call(rbind, values_list)

  rows <- list()
  idx <- 1L
  for (j in seq_len(nrow(meta))) {
    colj <- values_mat[, j]
    for (a in alpha_levels) {
      rows[[idx]] <- data.frame(
        stat = meta$stat[j],
        type = meta$type[j],
        T = "Inf",
        gamma = meta$gamma[j],
        q = meta$q[j],
        alpha = a,
        critical_value = as.numeric(stats::quantile(colj, probs = 1 - a, names = FALSE)),
        weight = meta$weight[j],
        nrep = nrep,
        n_train_grid = n_train_grid,
        n_open_grid = n_open_grid,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, rows)
  ks <- out[out$type == "KS", c("stat", "type", "T", "gamma", "q", "alpha", "critical_value", "nrep", "n_train_grid", "n_open_grid")]
  cvm <- out[out$type == "CvM", c("stat", "type", "T", "gamma", "q", "alpha", "critical_value", "weight", "nrep", "n_train_grid", "n_open_grid")]

  list(meta = meta, values = values_mat, ks = ks, cvm = cvm)
}

write_openend_critical_values <- function(sim,
                                          out_dir = "outputs",
                                          prefix = "openend_critical_values") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ks_path <- file.path(out_dir, paste0(prefix, "_ks.csv"))
  cvm_path <- file.path(out_dir, paste0(prefix, "_cvm.csv"))

  utils::write.csv(sim$ks, ks_path, row.names = FALSE)
  utils::write.csv(sim$cvm, cvm_path, row.names = FALSE)

  invisible(list(ks_path = ks_path, cvm_path = cvm_path))
}

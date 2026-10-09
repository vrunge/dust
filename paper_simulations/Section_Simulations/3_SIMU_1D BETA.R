# Figure 5 in the DUST paper: sensitivity to the penalty multiplier.
this_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(this_file)) sub("^--file=", "", this_file[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
Sys.setenv(DUST_SIM_DIR = dirname(this_file))
source(file.path(dirname(this_file), "ENGINE.R"))
profile <- sim_profile()
models <- c("gauss", "poisson")
n <- sim_grid(profile, 1e7, 5000L)
factors <- sim_grid(profile, exp(seq(log(0.001), log(20), length.out = 100)),
                    exp(seq(log(0.001), log(20), length.out = 5)))
repetitions <- sim_grid(profile, 100L, 2L)
tasks <- lapply(seq_len(length(models) * repetitions), function(i) {
  g <- expand.grid(replicate = seq_len(repetitions), model = models,
                   stringsAsFactors = FALSE)[i, ]
  list(model = g$model, replicate = g$replicate)
})
rows <- sim_parallel(tasks, function(task) {
  y <- sim_data(n, task$model, changes = 0L)
  out <- lapply(factors, function(factor) {
    fit <- sim_dust(y, task$model, factor * log(n))
    sim_row("beta", task$model, "dust", n, task$replicate, 0L,
            factor = factor, candidates = tail(fit$nb, 1L))
  })
  message("Completed penalty sweep: model=", task$model, ", replicate=", task$replicate)
  do.call(rbind, out)
}, seed = as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
sim_write(do.call(rbind, rows), "beta", profile)

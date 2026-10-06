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
backend <- Sys.getenv("DUST_SIM_BACKEND", "highway")
rows <- vector("list", length(models) * length(factors) * repetitions)
k <- 0L
set.seed(as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
for (model in models) for (replicate in seq_len(repetitions)) {
  y <- sim_data(n, model, changes = 0L)
  for (factor in factors) {
    fit <- sim_dust(y, model, factor * log(n), backend = backend)
    k <- k + 1L
    rows[[k]] <- sim_row("beta", model, "dust", n, replicate, 0L,
                         factor = factor, candidates = tail(fit$nb, 1L))
  }
  message("Completed penalty sweep: model=", model, ", replicate=", replicate)
}
sim_write(do.call(rbind, rows), "beta", profile)

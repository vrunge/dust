# Figure 1 in the DUST paper: candidate counts through time.
this_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(this_file)) sub("^--file=", "", this_file[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
Sys.setenv(DUST_SIM_DIR = dirname(this_file))
source(file.path(dirname(this_file), "ENGINE.R"))
profile <- sim_profile()
models <- c("gauss", "negbin")
sizes <- sim_grid(profile, c(1e4, 1e8), c(500, 2000))
repetitions <- sim_grid(profile, 100L, 2L)
backend <- Sys.getenv("DUST_SIM_BACKEND", "highway")
rows <- list()
set.seed(as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
for (model in models) for (n in sizes) for (replicate in seq_len(repetitions)) {
  y <- sim_data(n, model, changes = 0L)
  fit <- sim_dust(y, model, sim_penalty(n, model), backend = backend)
  ticks <- unique(as.integer(round(seq(1, n, length.out = min(1000L, n)))))
  rows[[length(rows) + 1L]] <- data.frame(
    experiment = "trajectory", model = model, algorithm = "dust", n = n,
    replicate = replicate, changes = 0L, factor = 1, t = ticks,
    candidates = fit$nb[ticks], time_sec = NA_real_
  )
  message("Completed trajectory: model=", model, ", n=", n, ", replicate=", replicate)
}
sim_write(do.call(rbind, rows), "nb_plot", profile)

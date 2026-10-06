# Figure 2 (candidate count) and Figure 3 (runtime) in the DUST paper.
this_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(this_file)) sub("^--file=", "", this_file[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
Sys.setenv(DUST_SIM_DIR = dirname(this_file))
source(file.path(dirname(this_file), "ENGINE.R"))
profile <- sim_profile()
models <- c("gauss", "poisson")
sizes <- sim_grid(profile, exp(seq(log(100), log(1e6), length.out = 100)), c(100, 5000, 7500, 10000))
sizes <- unique(as.integer(round(sizes)))
repetitions <- sim_grid(profile, 100L, 2L)
backend <- Sys.getenv("DUST_SIM_BACKEND", "highway")
rows <- vector("list", length(models) * length(sizes) * repetitions * 2L)
k <- 0L
set.seed(as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
for (model in models) for (replicate in seq_len(repetitions)) for (n in sizes) {
  y <- sim_data(n, model, changes = 0L)
  penalty <- sim_penalty(n, model)
  for (algorithm in c("dust", "fpop")) {
    result <- sim_measure(y, model, penalty, algorithm, backend = backend)
    fit <- result$fit
    candidates <- if (algorithm == "dust") tail(fit$nb, 1L) else NA_integer_
    k <- k + 1L
    rows[[k]] <- sim_row("complexity", model, algorithm, n, replicate,
                         0L, candidates = candidates, time_sec = result$time_sec)
  }
  if (k %% 20L == 0L) message("Completed ", k, " of ", length(rows), " runs")
}
sim_write(do.call(rbind, rows[seq_len(k)]), "regressions", profile)

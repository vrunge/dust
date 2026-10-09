# Figure 4 in the DUST paper: timing against the number of true changes.
this_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(this_file)) sub("^--file=", "", this_file[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
Sys.setenv(DUST_SIM_DIR = dirname(this_file))
source(file.path(dirname(this_file), "ENGINE.R"))
profile <- sim_profile()
models <- c("gauss", "negbin")
sizes <- sim_grid(profile, c(1000, 10000), c(500, 1000))
# The archived study varies segment spacing over 0, 10, ..., 50.
spacing_grid <- sim_grid(profile, seq(0, 50, by = 10), c(0, 50, 100, 200, 250, 400))
repetitions <- sim_grid(profile, 100L, 2L)
rows <- list()
set.seed(as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
for (model in models) for (n in sizes) for (spacing in spacing_grid) {
  ends <- if (spacing == 0L) n else c(seq.int(spacing, n, by = spacing), if (n %% spacing) n)
  changes <- length(unique(ends)) - 1L
  for (replicate in seq_len(repetitions)) {
    y <- sim_data_spacing(n, model, spacing)
    penalty <- sim_penalty(n, model)
    for (algorithm in c("dust", "fpop")) {
      result <- sim_measure(y, model, penalty, algorithm)
      k <- length(rows) + 1L
      rows[[k]] <- sim_row("density", model, algorithm, n, replicate,
                           changes, time_sec = result$time_sec)
    }
  }
  message("Completed timing grid: model=", model, ", n=", n,
          ", segment spacing=", spacing, ", true changes=", changes)
}
sim_write(do.call(rbind, rows), "density", profile)

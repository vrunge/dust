# Paper Figure 3 (fig:nb_plot): remaining candidate indices over time for
# Gaussian and negative binomial data without change, n = 1e4 and 1e8, 100
# repetitions each. Runs the paper design by default (DUST_SIM_PROFILE=smoke
# gives a quick test), writes results/<profile>/nb_plot.csv and draws the 2 x 2
# figure results/<profile>/figures/figure_1_pruning_capacity.{png,pdf}.
# Repetitions run on 50 cores, or every core when fewer exist (DUST_SIM_CORES).
this_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- if (length(this_file)) sub("^--file=", "", this_file[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
Sys.setenv(DUST_SIM_DIR = dirname(this_file))
source(file.path(dirname(this_file), "ENGINE.R"))
if (!nzchar(Sys.getenv("DUST_SIM_PROFILE"))) Sys.setenv(DUST_SIM_PROFILE = "paper")
profile <- sim_profile()
models <- c("gauss", "negbin")
sizes <- sim_grid(profile, c(1e4, 1e8), c(500, 2000))
repetitions <- sim_grid(profile, 100L, 2L)
backend <- Sys.getenv("DUST_SIM_BACKEND", "highway")
tasks <- lapply(seq_len(length(models) * length(sizes) * repetitions), function(i) {
  g <- expand.grid(replicate = seq_len(repetitions), n = sizes, model = models,
                   stringsAsFactors = FALSE)[i, ]
  list(model = g$model, n = g$n, replicate = g$replicate)
})
rows <- sim_parallel(tasks, function(task) {
  y <- sim_data(task$n, task$model, changes = 0L)
  fit <- sim_dust(y, task$model, sim_penalty(task$n, task$model), backend = backend)
  ticks <- unique(as.integer(round(seq(1, task$n, length.out = min(2000L, task$n)))))
  message("Completed trajectory: model=", task$model, ", n=", task$n, ", replicate=", task$replicate)
  data.frame(
    experiment = "trajectory", model = task$model, algorithm = "dust", n = task$n,
    replicate = task$replicate, changes = 0L, factor = 1, t = ticks,
    candidates = fit$nb[ticks], time_sec = NA_real_
  )
}, seed = as.integer(Sys.getenv("DUST_SIM_SEED", "20261005")))
results <- do.call(rbind, rows)
sim_write(results, "nb_plot", profile)
sim_plot_figure3(results, file.path(sim_output_dir(profile), "figures"))

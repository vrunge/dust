# DUST paper, Figure 3 (fig:nb_plot): number of remaining candidate indices
# over time, for data without change (mean over the repetitions and 95% interval).
#
# Self-contained: needs only the packages dust and ggplot2.
# Run it from the folder that contains it, in R with source("FIGURE_3.R")
# (the figure is then displayed) or in a terminal:
#   Rscript FIGURE_3.R                   simulate, then draw (models chosen below)
#   Rscript FIGURE_3.R poisson exp       same for 2 other models, no editing needed
#   Rscript FIGURE_3.R plot              redraw from the saved results
#   Rscript FIGURE_3.R plot poisson exp  redraw for these 2 models
# Outputs (same folder), for the 2 models <m1> and <m2>:
#   figure3_<m1>_<m2>_results.rds, figure3_<m1>_<m2>.png and .pdf, and
#   figure3_<m1>_<m2>_changes.csv: number of change points detected in each
#   setting with the chosen penalty (should be 0 or close to 0).

library(dust)

## Settings ---------------------------------------------------------------
# The 2 models shown (one per row). Choose 2 among the 8 models:
#   "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance"
# (in a terminal they can also be given as arguments, see the top of the file).
models      <- c("gauss", "negbin")
sizes       <- c(1e4, 1e8)
repetitions <- 100
n_times     <- 2000      # time points kept per run for the figure
seed        <- 20261005
# 50 cores, or every core when fewer exist (Windows cannot fork: 1 core).
# One run with n = 1e8 needs about 4 GB of memory: lower `cores` if needed.
cores <- if (.Platform$OS.type == "windows") 1 else min(50, parallel::detectCores())

# Simulation framework of the paper (Table 2), data without change: first
# parameter value of each model and scale factor c0 of the penalty
# beta = 2 c0 log(n). Binomial and negative binomial use size 10.
# c0 values checked with CALIBRATION.R at n = 1e4 (no false change, exact
# detection of 1 to 9 changes); Poisson uses the calibrated 1/3 (Table 2: 2/3).
paper <- list(
  gauss    = list(parameter = 0,   c0 = 1,    label = "Gaussian"),
  poisson  = list(parameter = 3,   c0 = 1/3,  label = "Poisson"),
  exp      = list(parameter = 1,   c0 = 3/4,  label = "Exponential"),
  geom     = list(parameter = 0.5, c0 = 2/3,  label = "Geometric"),
  bern     = list(parameter = 0.5, c0 = 2/3,  label = "Bernoulli"),
  binom    = list(parameter = 0.5, c0 = 1/6,  label = "Binomial"),
  negbin   = list(parameter = 0.5, c0 = 1/10, label = "Negative binomial"),
  variance = list(parameter = 1,   c0 = 1,    label = "Variance")
)
size <- 10

# Command-line arguments: optional "plot", then optionally the 2 models.
args <- commandArgs(trailingOnly = TRUE)
plot_only <- "plot" %in% args
if (length(setdiff(args, "plot"))) models <- setdiff(args, "plot")
if (length(models) != 2 || !all(models %in% names(paper)))
  stop("choose 2 models among: ", paste(names(paper), collapse = ", "))
tag <- paste(models, collapse = "_")   # output files: figure3_<model1>_<model2>.*

# Data are normalized with data_normalization_1D() before running DUST, as
# in the paper (e.g. Gaussian data divided by the estimated noise standard
# deviation, binomial and negative binomial counts divided by the size).
# The paper calibrated c0 on these normalized data.
simulate_data <- function(n, model) {
  y <- dataGenerator_1D(chpts = n, parameters = paper[[model]]$parameter,
                        nbTrials = size, nbSuccess = size, type = model)
  if (model %in% c("binom", "negbin")) data_normalization_1D(y, type = model, size = size)
  else data_normalization_1D(y, type = model)
}
penalty <- function(n, model) 2 * paper[[model]]$c0 * log(n)

# Runs task(1), ..., task(n_tasks) on `cores` workers. Each worker gets a new
# task as soon as it finishes the previous one (no idle cores), and the
# progress is printed after every `cores` completed tasks.
run_with_progress <- function(n_tasks, task) {
  start <- Sys.time()
  out <- vector("list", n_tasks)
  report <- function(done) {
    seconds <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    duration <- function(s) if (s < 120) sprintf("%.0f s", s) else sprintf("%.1f min", s / 60)
    cat(sprintf("%s  %d / %d done (%3.0f%%)  elapsed %s, about %s left\n",
                format(Sys.time(), "%H:%M:%S"), done, n_tasks, 100 * done / n_tasks,
                duration(seconds), duration(seconds / done * (n_tasks - done))))
    flush.console()
  }
  if (cores <= 1) {
    for (i in seq_len(n_tasks)) {
      out[[i]] <- try(task(i), silent = TRUE)
      if (i %% 10 == 0 || i == n_tasks) report(i)
    }
    return(out)
  }
  jobs <- list()      # running jobs
  index <- integer()  # task number of each running job
  next_task <- 1
  done <- 0
  while (done < n_tasks) {
    while (length(jobs) < cores && next_task <= n_tasks) {
      i <- next_task
      jobs[[length(jobs) + 1]] <- parallel::mcparallel(try(task(i), silent = TRUE))
      index <- c(index, i)
      next_task <- next_task + 1
    }
    finished <- parallel::mccollect(jobs, wait = FALSE, timeout = 1)
    if (is.null(finished)) next
    pids <- vapply(jobs, function(j) j$pid, numeric(1))
    for (pid in names(finished)) {
      k <- which(pids == as.numeric(pid))
      out[[index[k]]] <- finished[[pid]]
      done <- done + 1
      if (done %% cores == 0 || done == n_tasks) report(done)
    }
    keep <- !(pids %in% as.numeric(names(finished)))
    jobs <- jobs[keep]
    index <- index[keep]
  }
  out
}

## Simulation ---------------------------------------------------------------
run_simulation <- function() {
  grid <- expand.grid(n = sizes, model = models, replicate = seq_len(repetitions),
                      stringsAsFactors = FALSE)
  one_run <- function(i) {
    set.seed(seed + i)   # same results whatever the number of cores
    model <- grid$model[i]
    n <- grid$n[i]
    fit <- dust.1D(simulate_data(n, model), penalty = penalty(n, model), model = model)
    t <- unique(round(seq(1, n, length.out = n_times)))
    # `changes`: number of change points detected with the chosen penalty
    # (the data have none, so it should be 0 or close to 0).
    data.frame(model = model, n = n, replicate = grid$replicate[i], t = t,
               candidates = fit$nb[t], changes = length(fit$changepoints) - 1)
  }
  cat(nrow(grid), "runs on", cores, "cores\n")
  runs <- run_with_progress(nrow(grid), one_run)
  failed <- vapply(runs, function(r) !is.data.frame(r), logical(1))
  if (any(failed)) stop(sum(failed), " runs failed, for example: ", runs[[which(failed)[1]]])
  results <- do.call(rbind, runs)
  saveRDS(results, paste0("figure3_", tag, "_results.rds"))
  results
}

# Axis labels 0, 2.5K, ..., 100M (written here so no extra package is needed).
short_label <- function(x) {
  out <- format(x, trim = TRUE, scientific = FALSE, drop0trailing = TRUE)
  big <- !is.na(x) & abs(x) >= 1e6
  mid <- !is.na(x) & abs(x) >= 1e3 & !big
  out[big] <- paste0(format(x[big] / 1e6, trim = TRUE, drop0trailing = TRUE), "M")
  out[mid] <- paste0(format(x[mid] / 1e3, trim = TRUE, drop0trailing = TRUE), "K")
  out[is.na(x)] <- NA
  out
}

## Figure -------------------------------------------------------------------
draw_figure <- function(results) {
  library(ggplot2)
  q <- function(p) function(x) quantile(x, p, names = FALSE)
  keys <- results[c("model", "n", "t")]
  s <- aggregate(results["candidates"], keys, mean)
  s$lo <- aggregate(results["candidates"], keys, q(0.025))$candidates
  s$hi <- aggregate(results["candidates"], keys, q(0.975))$candidates
  s$model <- factor(s$model, models, gsub(" ", "~", vapply(models, function(m) paper[[m]]$label, "")))
  s$size <- factor(paste0("n == 10^", log10(s$n)))
  fig <- ggplot(s, aes(t, candidates)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "95% interval"), alpha = 0.35) +
    geom_line(aes(colour = "Mean"), linewidth = 0.9) +
    facet_grid(model ~ size, scales = "free", labeller = label_parsed) +
    scale_x_continuous(labels = short_label) +
    scale_fill_manual(NULL, values = c("95% interval" = "#78a9cf")) +
    scale_colour_manual(NULL, values = c("Mean" = "#155b8a")) +
    labs(x = "Time t", y = "Remaining candidate indices") +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", axis.text = element_text(size = 14),
          legend.key.width = grid::unit(1.6, "cm"),
          strip.background = element_rect(fill = "grey92"),
          panel.spacing.x = grid::unit(1.6, "lines"))
  ggsave(paste0("figure3_", tag, ".png"), fig, width = 10, height = 8, dpi = 180)
  ggsave(paste0("figure3_", tag, ".pdf"), fig, width = 10, height = 8)
  cat("Figure written: figure3_", tag, ".png and .pdf\n", sep = "")
  if (interactive()) print(fig)   # shown in RStudio or an R console
  invisible(fig)
}

## Check of the penalty ----------------------------------------------------
# Number of change points detected in each run (data without change), by
# model and length: printed and saved to figure3_<m1>_<m2>_changes.csv.
report_changes <- function(results) {
  if (!"changes" %in% names(results)) {
    cat("These results were saved before the change count was added: rerun the simulation.\n")
    return(invisible(NULL))
  }
  runs <- unique(results[c("model", "n", "replicate", "changes")])
  check <- do.call(rbind, lapply(split(runs, list(runs$model, runs$n), drop = TRUE), function(d)
    data.frame(model = d$model[1], n = d$n[1], runs = nrow(d),
               runs_with_0_change = sum(d$changes == 0), mean_changes = mean(d$changes),
               max_changes = max(d$changes))))
  check <- check[order(match(check$model, models), check$n), ]
  rownames(check) <- NULL
  cat("\nChange points detected with the chosen penalty (data without change):\n")
  print(check, digits = 3)
  write.csv(check, paste0("figure3_", tag, "_changes.csv"), row.names = FALSE)
  invisible(check)
}

## Main ---------------------------------------------------------------------
results <- if (plot_only) readRDS(paste0("figure3_", tag, "_results.rds")) else run_simulation()
report_changes(results)
draw_figure(results)

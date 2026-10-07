# DUST paper, Figure 3 for multivariate data: number of remaining candidate
# indices over time with dust.MD (method "QN"), for data without change.
# One plot per model, one curve per dimension (mean over the repetitions).
#
# Self-contained: needs only the packages dust and ggplot2.
# Run it from the folder that contains it, in R with source("FIGURE_3_MD.R")
# (the figure is then displayed) or in a terminal:
#   Rscript FIGURE_3_MD.R          simulate, then draw the figure
#   Rscript FIGURE_3_MD.R plot     redraw the figure from the saved results
# Outputs (same folder): figure3_MD_results.rds, figure3_MD.png and .pdf,
#   figure3_MD_changes.csv (number of change points detected, should be 0).

library(dust)

## Settings ---------------------------------------------------------------
n           <- 1e4                     # data length
repetitions <- 100                     # simulations for each model and dimension
# models (one plot each) among "gauss", "poisson", "exp", "geom", "bern",
# "binom", "negbin", "variance"
models      <- c("gauss", "poisson", "bern", "negbin")
dimensions  <- c(2, 3, 4, 5, 10)       # one curve per dimension
constraints <- "max"                   # "max" (= dimension) or a number (at most the dimension)
nbIterations <- 10                     # quasi-Newton steps for each pruning test
n_times     <- 2000                    # time points kept per run for the figure
seed        <- 20261009
# 50 cores, or every core when fewer exist (Windows cannot fork: 1 core).
cores <- if (.Platform$OS.type == "windows") 1 else min(50, parallel::detectCores())

# Simulation framework of the paper (Table 2), data without change: first
# parameter value of each model and scale factor c0 of the penalty
# beta = 2 c0 d log(n) for d time series. Binomial and negative binomial
# use size 10. Poisson uses the calibrated c0 = 1/3 (CALIBRATION.R).
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
stopifnot(all(models %in% names(paper)))
plot_only <- "plot" %in% commandArgs(trailingOnly = TRUE)

# d independent time series (rows), each normalized with data_normalization_1D
simulate_data <- function(n, d, model) {
  t(vapply(seq_len(d), function(i) {
    y <- dataGenerator_1D(chpts = n, parameters = paper[[model]]$parameter,
                          nbTrials = size, nbSuccess = size, type = model)
    if (model %in% c("binom", "negbin")) data_normalization_1D(y, type = model, size = size)
    else data_normalization_1D(y, type = model)
  }, numeric(n)))
}
penalty <- function(n, d, model) 2 * paper[[model]]$c0 * d * log(n)
n_constraints <- function(d) if (identical(constraints, "max")) d else min(constraints, d)

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
  grid <- expand.grid(model = models, d = dimensions, replicate = seq_len(repetitions),
                      stringsAsFactors = FALSE)
  one_run <- function(i) {
    set.seed(seed + i)   # same results whatever the number of cores
    model <- grid$model[i]
    d <- grid$d[i]
    fit <- dust.MD(simulate_data(n, d, model), penalty = penalty(n, d, model), model = model,
                   method = "QN", constraints = n_constraints(d), nbIterations = nbIterations)
    t <- unique(round(seq(1, n, length.out = n_times)))
    data.frame(model = model, d = d, replicate = grid$replicate[i], t = t,
               candidates = fit$nb[t], changes = length(fit$changepoints) - 1)
  }
  cat(nrow(grid), "runs on", cores, "cores\n")
  runs <- run_with_progress(nrow(grid), one_run)
  failed <- vapply(runs, function(r) !is.data.frame(r), logical(1))
  if (any(failed)) stop(sum(failed), " runs failed, for example: ", runs[[which(failed)[1]]])
  results <- do.call(rbind, runs)
  saveRDS(results, "figure3_MD_results.rds")
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
  s <- aggregate(results["candidates"], results[c("model", "d", "t")], mean)
  s$model <- factor(s$model, models, vapply(models, function(m) paper[[m]]$label, ""))
  # same colour for a dimension in all the plots
  dims <- sort(unique(s$d))
  s$dimension <- factor(paste("d =", s$d), paste("d =", dims))
  colours <- setNames(c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
                        "#56B4E9", "#F0E442", "#000000")[seq_along(dims)], paste("d =", dims))
  fig <- ggplot(s, aes(t, candidates, colour = dimension)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ model, scales = "free_y") +
    scale_x_continuous(labels = short_label) +
    scale_colour_manual(NULL, values = colours) +
    labs(x = "Time t", y = "Remaining candidate indices (mean)",
         title = paste0("dust.MD, method QN, ",
                        if (identical(constraints, "max")) "constraints = d" else paste("constraints =", constraints),
                        ", n = ", format(n, big.mark = ","), ", ", max(results$replicate), " repetitions")) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", axis.text = element_text(size = 14),
          legend.key.width = grid::unit(1.4, "cm"),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 15),
          panel.spacing.x = grid::unit(1.6, "lines"))
  ggsave("figure3_MD.png", fig, width = 11, height = 9, dpi = 180)
  ggsave("figure3_MD.pdf", fig, width = 11, height = 9)
  cat("Figure written: figure3_MD.png and .pdf\n")
  if (interactive()) print(fig)   # shown in RStudio or an R console
  invisible(fig)
}

## Check of the penalty ----------------------------------------------------
# number of change points detected (data without change, should be 0)
report_changes <- function(results) {
  runs <- unique(results[c("model", "d", "replicate", "changes")])
  check <- do.call(rbind, lapply(split(runs, list(runs$model, runs$d), drop = TRUE), function(r)
    data.frame(model = r$model[1], d = r$d[1], runs = nrow(r),
               runs_with_0_change = sum(r$changes == 0), mean_changes = mean(r$changes),
               max_changes = max(r$changes))))
  check <- check[order(match(check$model, models), check$d), ]
  rownames(check) <- NULL
  cat("\nChange points detected with the chosen penalty (data without change):\n")
  print(check, digits = 3)
  write.csv(check, "figure3_MD_changes.csv", row.names = FALSE)
  invisible(check)
}

## Main ---------------------------------------------------------------------
results <- if (plot_only) readRDS("figure3_MD_results.rds") else run_simulation()
report_changes(results)
draw_figure(results)

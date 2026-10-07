# DUST paper, Figure 4 (fig:nb_reg): number of remaining candidate indices at
# the end of DUST as a function of the data length, for data without change,
# with a log-log linear regression.
#
# Self-contained: needs only the packages dust and ggplot2.
# Run it from the folder that contains it, in R with source("FIGURE_4.R")
# (the figure is then displayed) or in a terminal:
#   Rscript FIGURE_4.R         simulate, then draw the figure
#   Rscript FIGURE_4.R plot    redraw the figure from figure4_results.rds
# Outputs (in the same folder): figure4_results.rds, figure4.png, figure4.pdf

library(dust)

## Settings ---------------------------------------------------------------
models      <- c("gauss", "poisson")
sizes       <- unique(round(exp(seq(log(1e2), log(1e6), length.out = 100))))
repetitions <- 100
seed        <- 20261006
# 50 cores, or every core when fewer exist (Windows cannot fork: 1 core).
cores <- if (.Platform$OS.type == "windows") 1 else min(50, parallel::detectCores())

# Data and penalties of the paper (beta = 2 c0 log(n), Table 2), no change:
# Gaussian N(0, 1) with c0 = 1; Poisson(3) with c0 = 1/3, as calibrated with
# CALIBRATION.R (Table 2: 2/3). As in the paper,
# data are normalized with data_normalization_1D() before running DUST
# (Gaussian: divided by the estimated noise sd; Poisson: divided by the mean).
simulate_data <- function(n, model) {
  y <- dataGenerator_1D(chpts = n, parameters = c(gauss = 0, poisson = 3)[[model]], type = model)
  data_normalization_1D(y, type = model)
}
penalty <- function(n, model) 2 * c(gauss = 1, poisson = 1/3)[[model]] * log(n)

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
# One task = one model and one repetition, over all the data lengths.
run_simulation <- function() {
  grid <- expand.grid(model = models, replicate = seq_len(repetitions),
                      stringsAsFactors = FALSE)
  one_task <- function(i) {
    set.seed(seed + i)   # same results whatever the number of cores
    model <- grid$model[i]
    candidates <- vapply(sizes, function(n) {
      fit <- dust.1D(simulate_data(n, model), penalty = penalty(n, model), model = model)
      as.numeric(tail(fit$nb, 1))
    }, numeric(1))
    data.frame(model = model, n = sizes, replicate = grid$replicate[i], candidates = candidates)
  }
  cat(nrow(grid), "tasks of", length(sizes), "lengths on", cores, "cores\n")
  runs <- run_with_progress(nrow(grid), one_task)
  failed <- vapply(runs, function(r) !is.data.frame(r), logical(1))
  if (any(failed)) stop(sum(failed), " tasks failed, for example: ", runs[[which(failed)[1]]])
  results <- do.call(rbind, runs)
  saveRDS(results, "figure4_results.rds")
  results
}

## Figure -------------------------------------------------------------------
draw_figure <- function(results) {
  library(ggplot2)
  results$log_n <- log(results$n)
  results$log_c <- log(results$candidates)
  panels <- lapply(c("gauss", "poisson"), function(m) {
    d <- results[results$model == m, ]
    band <- aggregate(log_c ~ log_n, d, function(x) quantile(x, c(0.025, 0.975), names = FALSE))
    band <- data.frame(log_n = band$log_n, lo = band$log_c[, 1], hi = band$log_c[, 2])
    fit <- lm(log_c ~ log_n, data = d)
    grid <- data.frame(log_n = sort(unique(d$log_n)))
    pred <- predict(fit, grid, interval = "prediction", level = 0.95)
    label <- sprintf("%s (slope %.2f)", c(gauss = "Gaussian", poisson = "Poisson")[[m]], coef(fit)[[2]])
    list(band = cbind(band, panel = label),
         fit = data.frame(grid, fit = pred[, "fit"], lwr = pred[, "lwr"], upr = pred[, "upr"], panel = label))
  })
  band <- do.call(rbind, lapply(panels, `[[`, "band"))
  fit <- do.call(rbind, lapply(panels, `[[`, "fit"))
  fig <- ggplot() +
    geom_ribbon(data = band, aes(log_n, ymin = lo, ymax = hi, fill = "95% interval"), alpha = 0.5) +
    geom_line(data = fit, aes(log_n, fit, colour = "Linear regression", linetype = "Linear regression"),
              linewidth = 1.2) +
    geom_line(data = fit, aes(log_n, lwr, colour = "95% prediction interval",
                              linetype = "95% prediction interval"), linewidth = 0.9) +
    geom_line(data = fit, aes(log_n, upr, colour = "95% prediction interval",
                              linetype = "95% prediction interval"), linewidth = 0.9) +
    facet_wrap(~ panel, scales = "free_y") +
    scale_fill_manual(NULL, values = c("95% interval" = "#9ecae1")) +
    scale_colour_manual(NULL, values = c("Linear regression" = "#0072B2",
                                         "95% prediction interval" = "#0072B2")) +
    scale_linetype_manual(NULL, values = c("Linear regression" = "solid",
                                           "95% prediction interval" = "dotted")) +
    labs(x = "log(data length n)", y = "log(remaining candidates at the end)") +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", axis.text = element_text(size = 14),
          legend.key.width = grid::unit(1.6, "cm"),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(size = 16, face = "bold"))
  ggsave("figure4.png", fig, width = 12, height = 5.5, dpi = 180)
  ggsave("figure4.pdf", fig, width = 12, height = 5.5)
  cat("Figure written: figure4.png and figure4.pdf\n")
  if (interactive()) print(fig)   # shown in RStudio or an R console
  invisible(fig)
}

## Main ---------------------------------------------------------------------
if (identical(commandArgs(trailingOnly = TRUE)[1], "plot")) {
  draw_figure(readRDS("figure4_results.rds"))
} else {
  draw_figure(run_simulation())
}

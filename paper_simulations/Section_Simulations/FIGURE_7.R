# DUST paper, Figure 7 (fig:nb_beta): number of remaining candidate indices at
# the end of DUST as a function of the penalty beta = a log(n), for data of
# length 1e7 without change, with the mean number of change points found.
#
# Self-contained: needs only the packages dust and ggplot2.
# Run it from the folder that contains it, in R with source("FIGURE_7.R")
# (the figure is then displayed) or in a terminal:
#   Rscript FIGURE_7.R         simulate, then draw the figure
#   Rscript FIGURE_7.R plot    redraw the figure from figure7_results.rds
# Outputs (in the same folder): figure7_results.rds, figure7.png, figure7.pdf

library(dust)

## Settings ---------------------------------------------------------------
models      <- c("gauss", "poisson")
n           <- 1e7
a_values    <- exp(seq(log(0.001), log(20), length.out = 100))   # beta = a log(n)
repetitions <- 100
seed        <- 20261007
# 50 cores, or every core when fewer exist (Windows cannot fork: 1 core).
# One task needs about 0.5 GB of memory.
cores <- if (.Platform$OS.type == "windows") 1 else min(50, parallel::detectCores())

# Data of the paper (Table 2), no change: Gaussian N(0, 1); Poisson(3),
# normalized with data_normalization_1D() before running DUST as in the paper
# (Gaussian: divided by the estimated noise sd; Poisson: divided by the mean).
simulate_data <- function(n, model) {
  y <- dataGenerator_1D(chpts = n, parameters = c(gauss = 0, poisson = 3)[[model]], type = model)
  data_normalization_1D(y, type = model)
}

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
# One task = one model and one repetition: one series, all the penalties.
run_simulation <- function() {
  grid <- expand.grid(model = models, replicate = seq_len(repetitions),
                      stringsAsFactors = FALSE)
  one_task <- function(i) {
    set.seed(seed + i)   # same results whatever the number of cores
    model <- grid$model[i]
    y <- simulate_data(n, model)
    out <- t(vapply(a_values, function(a) {
      fit <- dust.1D(y, penalty = a * log(n), model = model)
      c(candidates = tail(fit$nb, 1), changes = length(fit$changepoints) - 1)
    }, numeric(2)))
    data.frame(model = model, a = a_values, replicate = grid$replicate[i], out)
  }
  cat(nrow(grid), "tasks of", length(a_values), "penalties on", cores, "cores\n")
  runs <- run_with_progress(nrow(grid), one_task)
  failed <- vapply(runs, function(r) !is.data.frame(r), logical(1))
  if (any(failed)) stop(sum(failed), " tasks failed, for example: ", runs[[which(failed)[1]]])
  results <- do.call(rbind, runs)
  saveRDS(results, "figure7_results.rds")
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
  keys <- results[c("model", "a")]
  s <- aggregate(results["candidates"], keys, median)
  s$lo <- aggregate(results["candidates"], keys, q(0.025))$candidates
  s$hi <- aggregate(results["candidates"], keys, q(0.975))$candidates
  s$changes <- aggregate(results["changes"], keys, mean)$changes
  s$model <- factor(s$model, c("gauss", "poisson"), c("Gaussian", "Poisson"))
  # The right axis (mean number of change points) is rescaled onto the left one.
  k <- max(s$hi) / max(s$changes, 1)
  fig <- ggplot(s, aes(a)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "95% interval"), alpha = 0.5) +
    geom_line(aes(y = candidates, colour = "Median remaining candidates",
                  linetype = "Median remaining candidates"), linewidth = 1.2) +
    geom_line(aes(y = changes * k, colour = "Mean change points found (right axis)",
                  linetype = "Mean change points found (right axis)"), linewidth = 1.1) +
    facet_wrap(~ model) +
    scale_x_log10(breaks = 10^(-3:1),
                  labels = expression(10^-3, 10^-2, 10^-1, 10^0, 10^1)) +
    scale_y_continuous(name = "Median remaining candidates",
                       sec.axis = sec_axis(~ . / k, name = "Mean change points found",
                                           labels = short_label)) +
    scale_fill_manual(NULL, values = c("95% interval" = "#9ecae1")) +
    scale_colour_manual(NULL, values = c("Median remaining candidates" = "#0072B2",
                                         "Mean change points found (right axis)" = "#D55E00")) +
    scale_linetype_manual(NULL, values = c("Median remaining candidates" = "solid",
                                           "Mean change points found (right axis)" = "dashed")) +
    labs(x = expression("Penalty coefficient " * a * " (" * beta == a ~ log(n) * ")")) +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", axis.text = element_text(size = 14),
          legend.key.width = grid::unit(1.6, "cm"),
          axis.title.y.right = element_text(colour = "#D55E00"),
          axis.text.y.right = element_text(colour = "#D55E00"),
          strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(size = 16, face = "bold"))
  ggsave("figure7.png", fig, width = 12, height = 5.5, dpi = 180)
  ggsave("figure7.pdf", fig, width = 12, height = 5.5)
  cat("Figure written: figure7.png and figure7.pdf\n")
  if (interactive()) print(fig)   # shown in RStudio or an R console
  invisible(fig)
}

## Main ---------------------------------------------------------------------
if (identical(commandArgs(trailingOnly = TRUE)[1], "plot")) {
  draw_figure(readRDS("figure7_results.rds"))
} else {
  draw_figure(run_simulation())
}

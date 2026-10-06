# Rebuild the five simulation figures from RUN_ALL.R output CSV files.
args <- commandArgs(FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else sys.frame(1)$ofile
this_file <- normalizePath(gsub("~\\+~", " ", this_file))
here <- dirname(this_file)
Sys.setenv(DUST_SIM_DIR = here)
source(file.path(here, "ENGINE.R"))
profile <- sim_profile()
input <- file.path(here, "results", profile)
figures <- file.path(input, "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("install ggplot2 to draw figures")
if (!requireNamespace("patchwork", quietly = TRUE)) stop("install patchwork to assemble figures")
suppressPackageStartupMessages(library(ggplot2))
read_result <- function(name) utils::read.csv(file.path(input, paste0(name, ".csv")))
save_plot <- function(plot, name, width = 7, height = 5) {
  ggsave(file.path(figures, name), plot, width = width, height = height, dpi = 180)
}
theme_paper <- theme_bw(base_size = 12) + theme(legend.position = "bottom")
interval <- function(x) c(median = median(x, na.rm = TRUE),
                          lo = unname(quantile(x, 0.025, na.rm = TRUE)),
                          hi = unname(quantile(x, 0.975, na.rm = TRUE)))
summarise_interval <- function(data, value, groups) {
  key <- data[groups]
  x <- data[[value]]
  med <- aggregate(x, key, median, na.rm = TRUE)
  lo <- aggregate(x, key, function(z) unname(quantile(z, 0.025, na.rm = TRUE)))
  hi <- aggregate(x, key, function(z) unname(quantile(z, 0.975, na.rm = TRUE)))
  names(med)[ncol(med)] <- "median"
  names(lo)[ncol(lo)] <- "lo"
  names(hi)[ncol(hi)] <- "hi"
  merge(merge(med, lo, by = groups), hi, by = groups)
}

# Figure 1: remaining candidate indices over time (4 paper panels).
d <- read_result("nb_plot")
summary <- summarise_interval(d, "candidates", c("model", "n", "t"))
panels <- list()
for (model in c("gauss", "negbin")) for (n in unique(summary$n[summary$model == model])) {
  z <- summary[summary$model == model & summary$n == n, ]
  p <- ggplot(z, aes(t, median)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#78a9cf", alpha = .35) +
    geom_line(colour = "#155b8a") + scale_x_log10() + scale_y_log10() + theme_paper +
    labs(title = paste(toupper(model), "n =", format(n, scientific = FALSE)),
         x = "Time", y = "Remaining candidates")
  file <- paste0("pruning_capacity_", model, "_size_", format(n, scientific = TRUE), ".png")
  if (n == 10000) file <- paste0("pruning_capacity_", model, "_size_10000.png")
  if (n == 1e8) file <- paste0("pruning_capacity_", model, "_size_1e+08.png")
  save_plot(p, file)
  panels[[paste(model, n)]] <- p
}
save_plot(panels[[1]] + panels[[2]] + panels[[3]] + panels[[4]] + patchwork::plot_layout(ncol = 2),
          "figure_1_pruning_capacity.png", 12, 8)

# Figure 2: candidate count against n, with log-log fit and 95% prediction band.
d <- read_result("regressions")
d <- d[d$algorithm == "dust" & d$candidates > 0, ]
summary <- summarise_interval(d, "candidates", c("model", "n"))
panels <- list()
for (model in c("gauss", "poisson")) {
  z <- d[d$model == model, ]
  s <- summary[summary$model == model, ]
  fit <- lm(log(candidates) ~ log(n), data = z)
  grid <- data.frame(n = exp(seq(log(min(z$n)), log(max(z$n)), length.out = 150)))
  pr <- predict(fit, newdata = grid, interval = "prediction")
  grid$predicted <- exp(pr[, "fit"]); grid$lo <- exp(pr[, "lwr"]); grid$hi <- exp(pr[, "upr"])
  p <- ggplot(s, aes(n, median)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#78a9cf", alpha = .3) +
    geom_point(colour = "#155b8a", size = 1.5) + geom_line(data = grid, aes(n, predicted), inherit.aes = FALSE) +
    geom_line(data = grid, aes(n, lo), linetype = 2, inherit.aes = FALSE) +
    geom_line(data = grid, aes(n, hi), linetype = 2, inherit.aes = FALSE) +
    scale_x_log10() + scale_y_log10() + theme_paper + labs(title = toupper(model), x = "Data length", y = "Candidates at exit")
  save_plot(p, paste0("nb_complexity_label_", model, ".png"))
  panels[[model]] <- p
}
save_plot(panels[[1]] + panels[[2]], "figure_2_candidate_complexity.png", 12, 5)

# Figure 3: DUST and FPOP-compatible comparison runtime against n.
d <- read_result("regressions")
d <- d[d$model %in% c("gauss", "poisson") & d$time_sec > 0, ]
summary <- summarise_interval(d, "time_sec", c("model", "algorithm", "n"))
panels <- list()
for (model in c("gauss", "poisson")) {
  z <- d[d$model == model & d$n >= 3125, ]
  s <- summary[summary$model == model, ]
  predictions <- list()
  for (algorithm in unique(z$algorithm)) {
    q <- z[z$algorithm == algorithm, ]
    fit <- lm(log(time_sec) ~ log(n), data = q)
    grid <- data.frame(n = exp(seq(log(min(q$n)), log(max(q$n)), length.out = 150)))
    pr <- predict(fit, newdata = grid, interval = "prediction")
    predictions[[algorithm]] <- data.frame(grid, algorithm = algorithm, predicted = exp(pr[, "fit"]),
                                           lo = exp(pr[, "lwr"]), hi = exp(pr[, "upr"]))
  }
  pred <- do.call(rbind, predictions)
  p <- ggplot(s, aes(n, median, colour = algorithm, fill = algorithm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .16, colour = NA) +
    geom_line() + geom_point(size = 1.2) +
    geom_line(data = pred, aes(n, predicted, colour = algorithm), inherit.aes = FALSE) +
    geom_line(data = pred, aes(n, lo, colour = algorithm), linetype = 2, inherit.aes = FALSE) +
    geom_line(data = pred, aes(n, hi, colour = algorithm), linetype = 2, inherit.aes = FALSE) +
    scale_x_log10() + scale_y_log10() + theme_paper + labs(title = toupper(model), x = "Data length", y = "Elapsed time (seconds)")
  save_plot(p, paste0("time_complexity_label_", model, ".png"))
  panels[[model]] <- p
}
save_plot(panels[[1]] + panels[[2]], "figure_3_runtime_complexity.png", 12, 5)

# Figure 4: runtime against the number of true changes.
d <- read_result("density")
d <- d[d$time_sec > 0, ]
panels <- list()
for (model in c("gauss", "negbin")) for (n in sort(unique(d$n))) {
  z <- d[d$model == model & d$n == n, ]
  s <- summarise_interval(z, "time_sec", c("algorithm", "changes"))
  s$x <- log10(s$changes + 1)
  p <- ggplot(s, aes(x, median, colour = algorithm, fill = algorithm)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .16, colour = NA) + geom_line() + geom_point() +
    scale_y_log10() + theme_paper + labs(title = paste(toupper(model), "n =", n),
                                         x = "log10(true changes + 1)", y = "Elapsed time (seconds)")
  file_n <- if (n == 1000) "1000" else if (n == 10000) "10000" else as.character(n)
  save_plot(p, paste0("cpt_", model, "_", file_n, ".png"))
  panels[[paste(model, n)]] <- p
}
save_plot(panels[[1]] + panels[[2]] + panels[[3]] + panels[[4]] + patchwork::plot_layout(ncol = 2),
          "figure_4_runtime_by_changes.png", 12, 8)

# Figure 5: candidate count against the penalty factor.
d <- read_result("beta")
summary <- summarise_interval(d, "candidates", c("model", "factor"))
panels <- list()
for (model in c("gauss", "poisson")) {
  z <- summary[summary$model == model, ]
  p <- ggplot(z, aes(factor, median)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#78a9cf", alpha = .3) +
    geom_line(colour = "#155b8a") + geom_point(colour = "#155b8a", size = 1.3) +
    scale_x_log10() + scale_y_log10() + theme_paper + labs(title = toupper(model),
                                                           x = "Penalty factor", y = "Candidates at exit")
  size_label <- if (nrow(d) && unique(d$n[d$model == model]) == 1e7) "1e+07" else format(unique(d$n[d$model == model]), scientific = TRUE, trim = TRUE)
  save_plot(p, paste0("beta_nb_", model, "_size_", size_label, ".png"))
  panels[[model]] <- p
}
save_plot(panels[[1]] + panels[[2]], "figure_5_penalty_sensitivity.png", 12, 5)

message("Figures saved under: ", figures)

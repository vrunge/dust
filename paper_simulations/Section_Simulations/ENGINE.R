# Shared simulation helpers for the DUST paper figures.
# Uses the public API of the current dust package (0.2.0 or later).

paper_models <- c("gauss", "poisson", "exp", "geom", "bern", "negbin", "variance")
paper_c0 <- c(gauss = 1, poisson = 2/3, exp = 3/4, geom = 2/3,
              bern = 2/3, binom = 1/6, negbin = 1/10, variance = 1)
paper_parameters <- list(
  gauss = c(0, 1), poisson = c(3, 4), exp = c(1, 0.5),
  geom = c(0.5, 0.7), bern = c(0.5, 0.7), negbin = c(0.5, 0.7),
  variance = c(1, 2)
)

sim_profile <- function() {
  profile <- tolower(Sys.getenv("DUST_SIM_PROFILE", "smoke"))
  if (!profile %in% c("smoke", "paper")) stop("DUST_SIM_PROFILE must be 'smoke' or 'paper'")
  profile
}

sim_output_dir <- function(profile = sim_profile()) {
  path <- file.path(sim_script_dir(), "results", profile)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || is.na(x[1])) y else x

sim_data <- function(n, model, changes = 0L) {
  n <- as.integer(n)
  changes <- as.integer(changes)
  if (n < 2L || changes < 0L || changes >= n) stop("invalid series length or change count")
  ends <- if (changes == 0L) n else as.integer(floor(seq(0, n, length.out = changes + 2L)[-1L]))
  params <- paper_parameters[[model]]
  if (is.null(params)) stop("no generator parameters for model: ", model)
  params <- rep(params, length.out = length(ends))
  y <- dust::dataGenerator_1D(chpts = ends, parameters = params,
                              nbSuccess = 10, type = model)
  if (model == "negbin") y <- y / 10
  y
}

sim_data_spacing <- function(n, model, spacing) {
  n <- as.integer(n)
  spacing <- as.integer(spacing)
  if (spacing == 0L) return(sim_data(n, model, changes = 0L))
  ends <- seq.int(spacing, n, by = spacing)
  if (tail(ends, 1L) < n) ends <- c(ends, n)
  params <- rep(paper_parameters[[model]], length.out = length(ends))
  y <- dust::dataGenerator_1D(chpts = ends, parameters = params,
                              nbSuccess = 10, type = model)
  if (model == "negbin") y <- y / 10
  y
}

sim_penalty <- function(n, model, factor = 1) {
  scale <- if (model == "negbin") 10 else 1
  as.numeric(2 * paper_c0[[model]] * factor * log(n) / scale)
}

sim_dust <- function(y, model, penalty, method = "DUST", backend = "highway") {
  dust::dust.1D(y, penalty = penalty, model = model,
                method = method, backend = backend)
}

sim_compare <- function(y, model, penalty) {
  if (model == "gauss") {
    if (!requireNamespace("fpopw", quietly = TRUE)) stop("install fpopw to run timing comparisons")
    fit <- fpopw::Fpop(y, lambda = penalty)
    list(changepoints = fit$t.est, nb = NA_integer_)
  } else if (model %in% c("poisson", "negbin")) {
    if (!requireNamespace("gfpop", quietly = TRUE)) stop("install gfpop to run timing comparisons")
    graph <- gfpop::graph(type = "std", penalty = penalty)
    fit <- gfpop::gfpop(y, graph, type = model)
    list(changepoints = fit$changepoints, nb = NA_integer_)
  } else stop("no comparison implementation for model: ", model)
}

sim_measure <- function(y, model, penalty, algorithm = "dust", method = "DUST", backend = "highway") {
  call_fit <- function() {
    fit <- if (algorithm == "dust") sim_dust(y, model, penalty, method, backend) else sim_compare(y, model, penalty)
    fit
  }
  benchmark <- microbenchmark::microbenchmark(fit <- call_fit(), times = 1L)
  list(fit = fit, time_sec = unname(benchmark$time[1] / 1e9))
}

sim_row <- function(experiment, model, algorithm, n, replicate, changes,
                    factor = 1, t = NA_integer_, candidates = NA_integer_, time_sec = NA_real_) {
  data.frame(experiment, model, algorithm, n, replicate, changes, factor, t,
             candidates, time_sec)
}

sim_write <- function(data, name, profile = sim_profile()) {
  path <- file.path(sim_output_dir(profile), paste0(name, ".csv"))
  utils::write.csv(data, path, row.names = FALSE)
  message("Wrote ", path)
  invisible(path)
}

sim_grid <- function(profile, paper, smoke) if (profile == "paper") paper else smoke

# Worker processes for the candidate-count runners: 50 by default, or every
# available core when fewer exist. DUST_SIM_CORES overrides the default.
# Forked workers (parallel::mclapply) are not available on Windows, which
# runs the tasks sequentially.
sim_cores <- function() {
  configured <- Sys.getenv("DUST_SIM_CORES", "")
  cores <- if (nzchar(configured)) suppressWarnings(as.integer(configured)) else
    min(50L, parallel::detectCores(), na.rm = TRUE)
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (.Platform$OS.type == "windows") 1L else cores
}

# Applies f to each task on sim_cores() workers. Task i always uses seed
# seed + i, so the results do not depend on the number of workers.
sim_parallel <- function(tasks, f, seed) {
  run <- function(i) { set.seed(seed + i); f(tasks[[i]]) }
  cores <- sim_cores()
  message("Running ", length(tasks), " tasks on ", cores, " cores")
  out <- if (cores > 1L)
    parallel::mclapply(seq_along(tasks), run, mc.cores = cores, mc.preschedule = FALSE) else
    lapply(seq_along(tasks), run)
  failed <- vapply(out, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(failed)) stop("simulation task ", which(failed)[1], " failed: ",
                        paste(as.character(out[[which(failed)[1]]]), collapse = " "))
  out
}

# Median and 2.5% / 97.5% quantiles of a column by group.
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

# Paper Figure 3 (fig:nb_plot): remaining candidate indices over time, as one
# 2 x 2 grid (models in rows, data lengths in columns) in PNG and PDF, plus
# the four single panels that DUST.tex currently includes.
sim_plot_figure3 <- function(d, figures) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("install ggplot2 to draw figures")
  if (!requireNamespace("scales", quietly = TRUE)) stop("install scales to draw figures")
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  summary <- summarise_interval(d, "candidates", c("model", "n", "t"))
  panels <- list()
  # The paper shows each panel at 40% of the line width, so its text and lines
  # are drawn large enough to stay readable after that reduction.
  theme_fig3 <- theme_bw(base_size = 26) +
    theme(plot.title = element_text(size = 28, face = "bold"),
          axis.title = element_text(size = 26), axis.text = element_text(size = 22))
  for (model in c("gauss", "negbin")) for (n in sort(unique(summary$n[summary$model == model]))) {
    z <- summary[summary$model == model & summary$n == n, ]
    p <- ggplot(z, aes(t, median)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#78a9cf", alpha = .35) +
      geom_line(colour = "#155b8a", linewidth = 1.2) + scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + theme_fig3 +
      labs(title = paste(toupper(model), "n =", format(n, scientific = FALSE)),
           x = "Time", y = "Remaining candidates")
    file <- paste0("pruning_capacity_", model, "_size_", format(n, scientific = TRUE), ".png")
    if (n == 10000) file <- paste0("pruning_capacity_", model, "_size_10000.png")
    if (n == 1e8) file <- paste0("pruning_capacity_", model, "_size_1e+08.png")
    ggplot2::ggsave(file.path(figures, file), p, width = 7, height = 5, dpi = 180)
    panels[[paste(model, n)]] <- p
  }
  # Paper Figure 3 as one 2 x 2 grid: models in rows, data lengths in columns.
  size_label <- function(n) {
    k <- log10(n)
    if (abs(k - round(k)) < 1e-9) sprintf("n == 10^%d", as.integer(round(k))) else sprintf("n == %d", as.integer(n))
  }
  grid <- summary
  grid$model_label <- factor(c(gauss = "Gaussian", negbin = "Negative~binomial")[grid$model],
                             levels = c("Gaussian", "Negative~binomial"))
  lengths <- sort(unique(grid$n))
  grid$n_label <- factor(vapply(grid$n, size_label, ""), levels = vapply(lengths, size_label, ""))
  fig3 <- ggplot(grid, aes(t, median)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = "95% interval"), alpha = .35) +
    geom_line(aes(colour = "Median"), linewidth = 0.9) +
    facet_grid(model_label ~ n_label, scales = "free", labeller = label_parsed) +
    scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
    scale_fill_manual(NULL, values = c("95% interval" = "#78a9cf")) +
    scale_colour_manual(NULL, values = c("Median" = "#155b8a")) +
    labs(x = "Time t", y = "Remaining candidate indices") +
    theme_bw(base_size = 16) +
    theme(legend.position = "bottom", legend.text = element_text(size = 15),
          axis.text = element_text(size = 14),
          legend.key.width = grid::unit(1.6, "cm"),
          strip.text = element_text(size = 16, face = "bold"),
          strip.background = element_rect(fill = "grey92"),
          panel.spacing.x = grid::unit(1.6, "lines"), panel.spacing.y = grid::unit(0.8, "lines"))
  ggplot2::ggsave(file.path(figures, "figure_1_pruning_capacity.png"), fig3, width = 10, height = 8, dpi = 180)
  ggplot2::ggsave(file.path(figures, "figure_1_pruning_capacity.pdf"), fig3, width = 10, height = 8,
         device = grDevices::pdf)
  message("Figure 3 saved under: ", figures)
  invisible(fig3)
}

sim_script_dir <- function() {
  configured <- Sys.getenv("DUST_SIM_DIR", "")
  if (nzchar(configured)) return(normalizePath(configured))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(f) f$ofile %||% "", character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) dirname(normalizePath(gsub("~\\+~", " ", tail(ofiles, 1)))) else getwd()
}

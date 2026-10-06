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

sim_script_dir <- function() {
  configured <- Sys.getenv("DUST_SIM_DIR", "")
  if (nzchar(configured)) return(normalizePath(configured))
  frames <- sys.frames()
  ofiles <- vapply(frames, function(f) f$ofile %||% "", character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) dirname(normalizePath(gsub("~\\+~", " ", tail(ofiles, 1)))) else getwd()
}

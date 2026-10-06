library(dust)

# Four-panel version (1 row, 4 columns) of the DUST pruning illustration.
# Uses n = 5, without formula titles and with a single legend for all panels.
# PDF and PNG are written to the same folder as this script. The folder is found
# with Source in RStudio, source() or Rscript; otherwise the working directory is used.
script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(script_file))
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
out_dir <- if (length(script_file)) dirname(normalizePath(script_file)) else getwd()
out_name <- "DUST_pruning_1x4"

n <- 5
pen <- 2*log(n)
x_vals <- seq(-2, 2, length.out = 1000)

# Builds the n parabolas q_n^i(theta) and the lower envelope for one random dataset
make_panel <- function() {
  data <- dataGenerator_1D(c(n), parameters = c(0), sdNoise = 1, type = "gauss")
  res_dust <- dust.1D(data, penalty = pen, model = "gauss", method = "DUST")

  v1 <- c(-pen, res_dust$costQ)[1:n] + sum(data[1:n]^2)/2 + pen
  v2 <- n:1
  v3 <- rev(cumsum(rev(data[1:n])))

  q_vals <- sapply(1:n, function(i) v1[i] + v2[i] * x_vals^2/2 - v3[i] * x_vals)
  q_min <- apply(q_vals, 1, min)

  # Point where q^0 is above q^3, taken at the minimum of q^3 on that set
  indices <- which(q_vals[, 1] > q_vals[, 4])
  point <- indices[which.min(q_vals[indices, 4])]

  list(q_vals = q_vals, q_min = q_min, point = point)
}

draw_panel <- function(p) {
  q_vals <- p$q_vals
  q_min <- p$q_min

  par(mar = c(1, 1, 1, 1))
  matplot(x_vals, q_vals, type = "l", lty = 1, col = 1:n, lwd = 1.5,
          ylab = NA, xlab = NA, main = "",
          xlim = range(x_vals),
          ylim = c(min(q_min) - 1, min(q_min) + 10))

  # Thick segments where each parabola is the lower envelope
  for (i in 1:n) {
    is_min <- abs(q_vals[, i] - q_min) < 1e-8
    rle_min <- rle(is_min)
    idx <- cumsum(rle_min$lengths)
    start <- c(1, head(idx + 1, -1))
    for (j in which(rle_min$values)) {
      seg <- start[j]:idx[j]
      lines(x_vals[seg], q_vals[seg, i], col = i, lwd = 5)
    }
  }

  abline(h = min(q_min), lty = 3, lwd = 2, col = "gray40")
  points(x_vals[p$point], q_vals[p$point, 4], col = "red", pch = 19, cex = 2.5)
}

draw_figure <- function(panels) {
  layout(matrix(c(1:length(panels), rep(length(panels) + 1, length(panels))),
                nrow = 2, byrow = TRUE),
         heights = c(4.4, 1.0))

  for (p in panels) draw_panel(p)

  # Single legend shared by the four panels
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = lapply(0:(n - 1), function(i) bquote(q[.(n)]^.(i))),
         col = 1:n, lty = 1, lwd = 4, horiz = TRUE, bty = "n", cex = 2.6)
}

panels <- replicate(4, make_panel(), simplify = FALSE)

pdf(file.path(out_dir, paste0(out_name, ".pdf")), width = 16, height = 5.4)
draw_figure(panels)
invisible(dev.off())

png(file.path(out_dir, paste0(out_name, ".png")), width = 16, height = 5.4,
    units = "in", res = 200)
draw_figure(panels)
invisible(dev.off())

# Show the figure in the RStudio Plots pane when run interactively
if (interactive()) draw_figure(panels)

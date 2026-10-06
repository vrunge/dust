# Usage: Rscript RUN_ALL.R [smoke|paper]
args <- commandArgs(FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
this_file <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg[1])))
here <- dirname(normalizePath(this_file))
passed_args <- commandArgs(TRUE)
profile <- if (length(passed_args) && passed_args[1] %in% c("smoke", "paper")) passed_args[1] else "smoke"
scripts <- c("1_SIMU_1D REGRESSIONS.R", "2_SIMU_1D NB PLOT.R",
             "3_SIMU_1D BETA.R", "4_SIMU_1D DENSITY.R")
rscript <- file.path(R.home("bin"), "Rscript")
for (script in scripts) {
  status <- system2(rscript, shQuote(file.path(here, script)),
                    env = c(paste0("DUST_SIM_PROFILE=", profile), paste0("DUST_SIM_DIR=", here)), stdout = "", stderr = "")
  if (status != 0L) stop("simulation failed: ", script, " (status ", status, ")")
}
status <- system2(rscript, shQuote(file.path(here, "PLOT_FIGURES.R")),
                  env = c(paste0("DUST_SIM_PROFILE=", profile), paste0("DUST_SIM_DIR=", here)), stdout = "", stderr = "")
if (status != 0L) stop("figure rendering failed (status ", status, ")")

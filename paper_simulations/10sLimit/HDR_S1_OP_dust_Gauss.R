
library(dust)
library(dplyr)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  means <- rep(c(0, 10), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = means, type = "gauss", sdNoise = 1)
}

################################################################################
################################################################################
##### method = "OP": a single anchor for all nb_seg (dust.1D)
################################################################################
################################################################################
### method="OP" is now DUST_1D_OP_T/DUST_1D_HW_OP_T, a standalone engine with
### no active-set bookkeeping -- update_partition() is a flat loop over all
### s < t with no data-dependent early exit, so its instruction count is a
### pure function of n, identical regardless of nb_seg. That is confirmed
### (nb == t at every step, for every model, on every dataset tested this
### session) and is the actual fix: the old method="OP" silently routed
### through the general pruning engine and had wildly (100-1000x)
### nb_seg-dependent runtime from genuinely skipped work.
###
### A first attempt at this file called dust.1D(backend = "scalar") (as
### every sibling Poisson file in this directory does) and found wall-clock
### time was NOT flat despite the identical instruction count: at n=86500,
### nb_seg=10 took 11.1s, nb_seg=100 took 9.9s, nb_seg=1000 took 2.6s. That
### gap is real but specific to the scalar engine: its inner
### `if (c < minCost_t)` comparison is a genuine CPU branch, and how often
### it mispredicts depends on how "contested" the running minimum is as s
### sweeps -- degenerate/undetectable segmentations (fewer real contested
### updates) run faster than decisively-detected ones, for the same n.
### The Highway backend's op_step_scan does the same reduction with a vectorized,
### branch-free SIMD compare-and-select instead -- there is no scalar
### branch to mispredict, so wall-clock time really is flat, measured at
### two different n: 1.927/1.927/1.931s at n=86500, and 10.461/10.515/
### 10.560s at the final anchor below, across nb_seg=10/100/1000 both times
### (<1% spread, pure measurement noise). HW is also just faster outright
### for this method (unlike DUST/PELT in the sibling files, where HW only
### wins for Gauss -- for OP's simple dense unpruned scan it wins for both
### models, Poisson included -- see HDR_S1_OP_dust_Poisson.R).
###
### So, unlike every other method in this directory, OP needs only ONE
### n-range, reused unchanged across all three nb_seg blocks below -- the
### blocks are kept (rather than collapsed into one) specifically so the
### merged summary demonstrates, empirically, that mean_elapsed does not
### vary with nb_seg. Anchor: 197000 (measured).
###
### The means alternation was originally 0/1 (1 SD), where nb_seg=1000 at
### this anchor (197 points/segment) sat in the same detectability
### transition documented for Poisson below: n_detected measured 824 out
### of 1000 true changepoints, not 1000. Re-verified after switching to
### the stronger 0/10 alternation: detection is now full (n_detected=1000
### at nb_seg=1000 too) -- confirmed the anchor itself still lands at
### ~10.0s unchanged across all three nb_seg under the new parameters, as
### expected (OP's timing is provably data-independent). mean_n_detected
### is reported as-measured either way, not assumed.
################################################################################

n_range <- c(185000, 191000, 197000, 203000, 209000, 215000)

run_op_block <- function(nb_seg) {
  results <- data.frame()

  for (n in n_range) {

    print(n)

    for (j in 1:10) {

      data <- build_data(n, nb_seg)

      timing <- system.time(
        res <- dust.1D(data, penalty = 2 * log(n), model = "gauss", method = "OP", threads = 1)
      )

      n_detected <- length(res$changepoints)

      results <- rbind(
        results,
        data.frame(
          n = n, run = j,
          user = timing[["user.self"]], system = timing[["sys.self"]], elapsed = timing[["elapsed"]],
          n_detected = n_detected
        )
      )

      print(j)
    }
  }

  results %>%
    group_by(n) %>%
    summarise(
      N_runs = n(),
      mean_elapsed = mean(elapsed),
      sd_elapsed = sd(elapsed),
      se_elapsed = sd_elapsed / sqrt(N_runs),
      CI_lower = mean_elapsed - qt(0.975, df = N_runs - 1) * se_elapsed,
      CI_upper = mean_elapsed + qt(0.975, df = N_runs - 1) * se_elapsed,
      mean_n_detected = mean(n_detected),
      .groups = "drop"
    ) %>%
    mutate(nb_seg = nb_seg)
}

summary_results10 <- run_op_block(10)
print(summary_results10)

summary_results100 <- run_op_block(100)
print(summary_results100)

summary_results1000 <- run_op_block(1000)
print(summary_results1000)

################################################################################
################################################################################
##### Merged summary across all nb_seg
################################################################################
################################################################################

summary_results_OP_dust_Gauss <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_OP_dust_Gauss)

write.csv(summary_results_OP_dust_Gauss,
          file = "paper_simulations/10sLimit/results/HDR_S1_OP_dust_Gauss_summary.csv",
          row.names = FALSE)

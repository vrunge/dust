
library(dust)
library(dplyr)

build_data <- function(n, nb_seg) {
  cpts <- floor(seq(from = 1/nb_seg, to = 1, by = 1/nb_seg) * n)
  rates <- rep(c(4, 20), length.out = nb_seg)
  dataGenerator_1D(chpts = cpts, parameters = rates, type = "poisson")
}

################################################################################
################################################################################
##### method = "OP": a single anchor for all nb_seg (dust.1D, backend = "highway")
################################################################################
################################################################################
### method="OP" is now DUST_1D_OP_T/DUST_1D_HW_OP_T, a standalone engine with
### no active-set bookkeeping -- update_partition() is a flat loop over all
### s < t with no data-dependent early exit, so its instruction count is a
### pure function of n, identical regardless of nb_seg. That is confirmed
### (nb == t at every step, for every model, on every dataset tested this
### session) and is the actual fix: the old method="OP" silently routed
### through the general pruning engine, and at nb_seg=1000 specifically
### that engine's active set genuinely never shrank either (mean(nb)/n =
### 0.5000 exactly, confirmed on 8+ seeds at n=84875, n=100000, and
### n=400000) -- but at nb_seg=10/100 it pruned aggressively, giving
### 100-1000x faster runtime purely from skipped work. That nb_seg-
### dependent skipped work is what's gone now.
###
### A first attempt at this file kept dust.1D(backend = "scalar") (the
### convention every sibling Poisson file in this directory uses for
### DUST/PELT) and found wall-clock time was NOT flat despite the
### identical instruction count: at n=60000, nb_seg=10 took 16.0s,
### nb_seg=100 took 9.6s, nb_seg=1000 took 4.3s. That gap is real but
### specific to the scalar engine: its inner `if (c < minCost_t)`
### comparison is a genuine CPU branch, and how often it mispredicts
### depends on how "contested" the running minimum is as s sweeps --
### degenerate/undetectable segmentations (fewer real contested updates,
### as at nb_seg=1000 here, see the known limitation below) run faster
### than decisively-detected ones (nb_seg=10/100), for the same n.
### The Highway backend's op_step_scan does the same reduction with a vectorized,
### branch-free SIMD compare-and-select instead -- there is no scalar
### branch to mispredict, so wall-clock time really is flat, measured at
### two different n: ~4.9s (4.891/4.968/4.984) at n=60000, and ~10.3s
### (10.311/10.303/10.248) at the final anchor below, across
### nb_seg=10/100/1000 both times (<1% spread, pure measurement noise).
### HW is also just faster outright here: at n=60000, scalar took
### 16.6/10.1/4.8s versus HW's flat ~4.9s. This is the opposite of the
### DUST/PELT sibling files' choice (they use scalar for Poisson because
### HW was NOT fastest there) -- but OP's workload is a simple dense
### unpruned scan, not a branchy pruning test, and for that, HW wins for
### both models, so this file switches to backend = "highway" despite its DUST/PELT
### siblings using scalar.
###
### So, unlike every other method in this directory, OP needs only ONE
### n-range, reused unchanged across all three nb_seg blocks below -- the
### blocks are kept (rather than collapsed into one) specifically so the
### merged summary demonstrates, empirically, that mean_elapsed does not
### vary with nb_seg. Anchor: 86000 (measured).
###
### KNOWN, CONFIRMED-NOT-A-BUG LIMITATION (unchanged from before the engine
### and backend changes, and still holds after switching rates from 4/8 to
### 4/20): at nb_seg=1000, <100 points/segment is genuinely too weak a
### signal to detect under the default penalty -- n_detected will read 1
### (or occasionally 2-3; re-verified at n_detected=2 under 4/20, still
### effectively collapsed), not 1000. PELT and DUST give the identical
### costQ/changepoints on this exact data (mathematically required: all
### three recursions are pruning-independent), so this is a property of
### the statistical problem, not of OP or of the chosen backend.
################################################################################

n_range <- c(80000, 83000, 86000, 89000, 92000, 95000)

run_op_block <- function(nb_seg) {
  results <- data.frame()

  for (n in n_range) {

    print(n)

    for (j in 1:10) {

      data <- data_normalization_1D(build_data(n, nb_seg), type = "poisson")

      timing <- system.time(
        res <- dust.1D(data, penalty = 2 * log(n), model = "poisson", method = "OP", backend = "highway")
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

summary_results_OP_dust_Poisson <- bind_rows(summary_results10, summary_results100, summary_results1000)

print(summary_results_OP_dust_Poisson)

write.csv(summary_results_OP_dust_Poisson,
          file = "paper_simulations/10sLimit/results/HDR_S1_OP_dust_Poisson_summary.csv",
          row.names = FALSE)

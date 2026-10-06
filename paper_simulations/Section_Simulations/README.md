# DUST paper simulations

This directory contains four simulation runners and one plotting script for the five figures in the univariate simulation study in `PAPER_DUST/DUST.tex`. The paper-scale configuration writes into `results/paper/`; the default smoke configuration writes into `results/smoke/`.

From this directory run:

```sh
Rscript RUN_ALL.R smoke
Rscript RUN_ALL.R paper
```

`smoke` uses small series, three lengths for the complexity figures, five penalty values, and two repetitions. It checks that the full data-to-CSV-to-figure workflow runs. `paper` uses the design in the paper: 100 repetitions; 100 log-spaced lengths from 100 to 1,000,000; trajectory lengths 10,000 and 100,000,000; 100 penalty factors; and 1,000 and 10,000 point timing experiments. The change-density simulation varies segment spacing over 0, 10, ..., 50 observations, as in the archived runner; its plotted horizontal coordinate is `log10(true changes + 1)`. The paper profile is computationally and memory intensive, especially the 100 trajectories of length 100,000,000. Use `DUST_SIM_BACKEND=scalar` to select the scalar implementation; the default is `highway`. The random DUST evaluator discussed in the older paper text is no longer part of the current 1D API.

The scripts use `dust::dataGenerator_1D()` and `dust::dust.1D()` with `method = "DUST"`. Gaussian and Poisson runtime comparisons use `fpopw` and `gfpop`, respectively; the latter provides the Poisson FPOP-style comparator. Install `dust`, `fpopw`, `gfpop`, `microbenchmark`, `ggplot2`, and `patchwork` before running. Results record elapsed time in seconds, so timing curves depend on the machine, compiler, package versions, and backend.

## Figure map

| Paper figure | Runner | Data file | Panel image files |
| --- | --- | --- | --- |
| `fig:nb_plot` | `2_SIMU_1D NB PLOT.R` | `nb_plot.csv` | `pruning_capacity_{gauss,negbin}_size_{10000,1e+08}.png` |
| `fig:nb_reg` | `1_SIMU_1D REGRESSIONS.R` | `regressions.csv` | `nb_complexity_label_{gauss,poisson}.png` |
| `fig:time_reg` | `1_SIMU_1D REGRESSIONS.R` | `regressions.csv` | `time_complexity_label_{gauss,poisson}.png` |
| `fig:timeExecution` | `4_SIMU_1D DENSITY.R` | `density.csv` | `cpt_{gauss,negbin}_{1000,10000}.png` |
| `fig:nb_beta` | `3_SIMU_1D BETA.R` | `beta.csv` | `beta_nb_{gauss,poisson}_size_1e+07.png` |

`PLOT_FIGURES.R` also writes one assembled `figure_1_...png` through `figure_5_...png` image for each paper figure. It can be rerun by itself after the four CSV files exist, for example:

```sh
DUST_SIM_PROFILE=paper Rscript PLOT_FIGURES.R
```

The five historic `cpt6_*.csv` files are retained for reference. They use the old column layout and do not contain the input needed to redraw the five paper figures. The archived paper images and paper text were used to identify the five figure targets; no complete original plotting pipeline was present in this directory. The scripts here regenerate all five figure designs with the current package. Exact pixel and numerical reproduction of the published plots is not guaranteed because the dust API and algorithms have changed since those archived runs.

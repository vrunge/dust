<a id="top"></a>

# dust Vignette

### Vincent Runge and Simon Querné

#### LaMME, Evry University

<center>
<img src="man/figures/dust.png" alt="" style="width:30%;"/>
</center>

## Introduction

The `dust` package detects multiple change points in a time series by
minimizing a penalized segment cost. It uses optimal partitioning with a
pruning rule that discards candidates for the last change point. The
DUST rule evaluates a dual decision function; see the [DUST
paper](https://doi.org/10.48550/arXiv.2507.02467).

This version covers eight one-parameter models for 1D or independent
multivariate data, and a Gaussian model with changes in both mean and
variance. They can be run on a complete series or through an Rcpp object
that receives successive batches. The optional Highway backend
accelerates supported calculations when available at build time.

> [Quick start](#start)

> [Rcpp Object Structure](#rcpp)

> [Models And Data Generators](#Models)

> [dust Algorithms](#dust1D)

> [Pruning Capacity](#pruning)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="start"></a>

## Quick start

### Installing the dust Rcpp package

After the repository has been updated, install it with:

    remotes::install_github("vrunge/dust")
    library(dust)

### A simple 1D example

The following Gaussian series has changes in its mean. `dust.1D()`
returns the estimated change points, including the final index of the
series.

    set.seed(11)
    y <- c(rnorm(80), rnorm(80, mean = 2), rnorm(80, mean = -1))
    fit <- dust.1D(y, model = "gauss", method = "DUST")
    fit$changepoints

The default penalty is `2 * log(length(y))`. The result also contains
`costQ`, the optimal penalized cost through each observation; `nb`, the
number of active candidates over time; `lastIndexSet`, the candidates
still active at the end; and the backend that ran.

### Independent multivariate series

`dust.MD()` takes a matrix with one component per row and one
observation per column. It adds the component costs and finds change
points shared across components. Its default penalty is
`2 * nrow(data) * log(ncol(data))`.

    set.seed(13)
    multi <- rbind(c(rnorm(80), rnorm(80, mean = 2)),
                   c(rnorm(80), rnorm(80, mean = -1)))
    fit_md <- dust.MD(multi, model = "gauss", constraints = 2)
    fit_md$changepoints

The `constraints` largest still active indices below each candidate are used for its pruning test. The available methods are `"coordinateDescent"` (default), `"iterative"` (projected gradient with backtracking), `"QN"` (Armijo/BFGS), `"randomEval"`, `"exact"`, `"PELT"`, and unpruned `"OP"`. For Gaussian data, `"exact"` maximizes the joint decision over all selected constraints by checking a short sequence of stationary faces, then enumerating faces if needed. The fallback can take exponential time in `constraints`. For other models, `"exact"` uses PELT.

`"iterative"` and `"QN"` jointly search the selected constraints, for at most `nbIterations` optimizer iterations per candidate (effective default 10). If `nbIterations` is omitted, `epsilon` can stop a numerical search when the gain in the decision function is at most the given threshold, with a cap of 1000 iterations. An explicit `nbIterations` takes priority. Every trial is checked for model-domain feasibility. Pruning requires a strictly positive decision value with a numerical tolerance; an exhausted budget or unsuccessful search retains the candidate. A finite search does not guarantee finding the maximum. Every search uses the available selected earlier indices, up to `constraints`, even when fewer are available.

One `nbIterations` unit depends on the method. The budget restarts for each candidate at each time point; it is not a limit for the full segmentation. Every pruned method first checks the decision at zero (the PELT test) and stops immediately when it finds a feasible positive value.

| Method | One `nbIterations` unit | `epsilon` stopping |
|:--|:--|:--|
| `coordinateDescent` | One full sweep through the selected multiplier coordinates. | After a sweep, stop if the normalized decision value increased by at most `epsilon`. |
| `iterative` | One projected-gradient step, with up to 60 backtracking trials. | After an accepted step, stop if the gain is at most `epsilon`. |
| `QN` | One quasi-Newton attempt, with up to 60 backtracking trials and, if needed, a gradient fallback with up to 60 more. | After an accepted step, stop if the gain is at most `epsilon`. |
| `randomEval` | One random direction and radius, with up to 40 feasibility reductions. | Unavailable: consecutive random draws have no convergence meaning. |
| `exact` | The budget is ignored. Gaussian data use the finite joint active-face solver; other models use PELT. | Ignored. |
| `PELT` / `OP` | The budget is ignored. PELT tests only zero; OP does no pruning. | Ignored. |

When `nbIterations` is omitted, its effective value is 10. If `epsilon` is supplied for one of the three numerical optimizers, the cap becomes 1000 instead. Supplying `nbIterations` explicitly disables `epsilon`, even when both arguments are present. The gain threshold is absolute and uses the normalized decision function; it does not certify that the maximum has been found.

### Changes in mean and variance

For a Gaussian series in which both parameters can change, use
`dust.meanVar()`. Its default penalty is `4 * log(length(y))`. The
`"1D"` and `"2D"` methods use one and two pruning constraints,
respectively, for the same segmentation objective.

    set.seed(12)
    z <- c(rnorm(100, sd = 0.6), rnorm(100, sd = 2),
           rnorm(100, mean = 1, sd = 0.8))
    fit_mv <- dust.meanVar(z, method = "1D")
    fit_mv$changepoints
    dust.meanVar(z, method = "2D")$changepoints

Each segment must contain at least two nonidentical observations to have
finite Gaussian mean-and-variance cost. The separate 1D `"variance"`
model keeps the mean fixed at zero.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="rcpp"></a>

## Rcpp Object Structure

The object forms accept data in batches. Append observations, call
`update_partition()`, then read the current result with
`get_partition()`. If the final series length is known, supply a penalty
on the first append: otherwise its default is calculated from the first
batch size.

    penalty <- 2 * log(length(y))
    obj <- dust.object.1D(model = "gauss", method = "DUST")
    obj$append_data(y[1:120], penalty)
    obj$update_partition()
    obj$append_data(y[121:240], NULL)
    obj$update_partition()
    obj$get_partition()$changepoints

`dust.object.meanVar()` has the same append/update/get pattern:

    obj_mv <- dust.object.meanVar(method = "2D")
    obj_mv$append_data(z[1:150], 4 * log(length(z)))
    obj_mv$update_partition()
    obj_mv$append_data(z[151:300], NULL)
    obj_mv$update_partition()
    obj_mv$get_partition()$changepoints

`dust.object.MD()` accepts matrix batches in the same way:

    obj_md <- dust.object.MD(model = "gauss", constraints = 2)
    obj_md$append_data(multi[, 1:100, drop = FALSE],
                       2 * nrow(multi) * log(ncol(multi)))
    obj_md$update_partition()
    obj_md$append_data(multi[, 101:160, drop = FALSE], NULL)
    obj_md$update_partition()
    obj_md$get_partition()$changepoints

These objects also provide `get_info()`, including the backend actually
used.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="Models"></a>

## Models And Data Generators

`dataGenerator_1D()` generates univariate examples. `dataGenerator_MD()`
uses the same eight models to generate independent components with shared
change points. Supply segment parameters as a matrix or data frame with
segments in rows and components in columns; the result has components in
rows and time in columns, as expected by `dust.MD()`. Noise levels and
count-model sizes can be shared or specified per component.
`data_normalization_1D()` applies the model-specific transformations
needed before segmentation.

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;"><code>model</code></th>
<th style="text-align: left;">Segment parameter</th>
<th style="text-align: left;">Data supplied to <code>dust.1D()</code> or
each row of <code>dust.MD()</code></th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align: left;"><code>gauss</code></td>
<td style="text-align: left;">Mean, known unit variance</td>
<td style="text-align: left;">Finite real values</td>
</tr>
<tr>
<td style="text-align: left;"><code>poisson</code></td>
<td style="text-align: left;">Rate</td>
<td style="text-align: left;">Nonnegative counts</td>
</tr>
<tr>
<td style="text-align: left;"><code>exp</code></td>
<td style="text-align: left;">Rate</td>
<td style="text-align: left;">Strictly positive values</td>
</tr>
<tr>
<td style="text-align: left;"><code>geom</code></td>
<td style="text-align: left;">Success probability</td>
<td style="text-align: left;">Trials through the first success</td>
</tr>
<tr>
<td style="text-align: left;"><code>bern</code></td>
<td style="text-align: left;">Success probability</td>
<td style="text-align: left;">Values in <code>[0, 1]</code></td>
</tr>
<tr>
<td style="text-align: left;"><code>binom</code></td>
<td style="text-align: left;">Success probability</td>
<td style="text-align: left;">Counts divided by the known number of
trials</td>
</tr>
<tr>
<td style="text-align: left;"><code>negbin</code></td>
<td style="text-align: left;">Success probability</td>
<td style="text-align: left;">Counts divided by the known size</td>
</tr>
<tr>
<td style="text-align: left;"><code>variance</code></td>
<td style="text-align: left;">Variance, zero mean</td>
<td style="text-align: left;">Nonzero residuals</td>
</tr>
</tbody>
</table>

For Binomial and Negative Binomial observations, pass the known trial
count or size to `data_normalization_1D(size = ...)` for each component.
When translating an objective based on raw counts, divide the penalty by
the same known value. For Gaussian data, divide by a noise estimate if
using the default penalty.

The mean-and-variance model is available through `dust.meanVar()` and
`dust.object.meanVar()`; it does not require a separate data generator.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="dust1D"></a>

## dust Algorithms

For one-parameter models, `method = "DUST"` evaluates the dual decision
test. `"DUSTib"` uses an inequality certificate, `"PELT"` uses the PELT
rule, and `"OP"` scans all candidate endpoints without pruning. All four
methods minimize the same segment objective.

For multivariate data, `"coordinateDescent"`, `"iterative"`, `"QN"`, and
`"randomEval"` seek a feasible positive value of the multivariate decision
function. A finite search may leave some candidates active. The `nbIterations`
argument controls the number of coordinate sweeps, optimizer iterations,
or random evaluations. Regression tests compare all prefix costs and
changepoints with `"PELT"` and unpruned `"OP"`. Additional checks against
analytic decision maxima can be run from the package source directory
with `Rscript tools/check-md-search.R` and `Rscript tools/check-md-exact.R`.

    set.seed(5)
    counts <- dataGenerator_1D(chpts = c(60, 120, 180),
                               parameters = c(2, 8, 3), type = "poisson")
    dust.1D(counts, model = "poisson", method = "DUST")$changepoints

The default `backend = "highway"` uses SIMD when Highway was available
at build time and otherwise selects the scalar implementation.
`backend = "scalar"` selects it explicitly. Installation detects Highway
through `pkg-config` on platforms where it is available. Set
`DUST_FORCE_SCALAR=1` before installation to request the scalar build on
any platform.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="pruning"></a>

## Pruning Capacity

The `nb` field records the number of active candidates after pruning at
each time point. It is useful for inspecting how much work the algorithm
saves on a particular series. The final `lastIndexSet` and the
trajectory `nb` can differ between the scalar and Highway backends
because they visit candidates in a different order; the optimal cost is
the quantity to compare across backends.

The standalone research scripts are kept in `paper_simulations/` and are
excluded from the CRAN package archive.

[(Back to Top)](#top)

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

This version covers eight one-parameter models and a Gaussian model with
changes in both mean and variance. Both can be run on a complete series
or through an Rcpp object that receives successive batches. The optional
Highway backend accelerates supported calculations when available at
build time.

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

Both objects also provide `get_info()`, including the backend actually
used.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="Models"></a>

## Models And Data Generators

`dataGenerator_1D()` generates one-parameter examples.
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
<th style="text-align: left;">Data supplied to
<code>dust.1D()</code></th>
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
count or size to `data_normalization_1D(size = ...)`. When translating
an objective based on raw counts, divide the penalty by the same known
value. For Gaussian data, divide by a noise estimate if using the
default penalty.

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

<a id="top"></a>

# dust Vignette

### Vincent Runge and Simon Querné

#### LaMME, Evry University, November 24, 2025

<center>
<img src="man/figures/dust.png" alt="" style="width:30%;"/>
</center>

> [Quick start](#start)

> [Models And Data Generators](#Models)

> [dust 1D and MD Algorithms](#dust1D)

> [dust 2D Algorithms](#dust2D)

> [Rcpp Structure and new cost integration](#rcpp)

> [Pruning Capacity](#pruning)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

## Introduction

The `dust` package contains methods **for detecting multiple
change-points within time-series** based on the optimal partitioning
algorithm, which is a dynamic programming (DP) algorithm. Our algorithms
optimize a penalized likelihood and the DP algorithm is encoded with
pruning rules for reducing execution time. The novelty of the `dust`
package consists in its pruning step. We use a **new pruning rule**,
different from the two standard ones: [PELT rule
(2012)](\doi%7B10.1080/01621459.2012.737745%7D) and [FPOP rule
(2017)](\doi%7B10.1007/s11222-016-9636-3%7D).

We called this method the **DUST** pruning rule, standing for
**Du**ality **S**imple **T**est. This method is based on considering
some optimization problems under inequality constraints and its dual for
discarding indices in the search for the last change-point index.

We propose change-point detection algorithms working with different
models derived from the **exponential family** (Gauss, Poisson,
Exponential…).

We provide a polymorphic structure in Rcpp that enables users **to
easily add custom cost functions of their choice**. Detailed
instructions are available [in this Section](#rcpp), and users only need
to define a few Rcpp functions: the minimal cost function, the dual
function, and its derivative, along with some domain information for the
dual.

Various tests and simulations are included in this **README file** and
the **simulations folder**, demonstrating that the DUST dual approach is
**highly efficient across all regimes** (whether detecting few or many
changes) and shows improved computational efficiency compared to PELT
and FPOP. Additionally, unlike these methods, DUST can also reduce
computation time for multivariate cost functions (See Section
[pruning](#pruning)).

<a id="start"></a>

## Quick start

### Installing the dust Rcpp package

**REQUIREMENTS:** - R &gt;= 3.4 - devtools :
`install.packages('devtools')`

The package can then be installed from the github repo with the
following command:

    devtools::install_github("vrunge/dust")

and imported with:

    library(dust)

### A simple example

We generate some 1D time series of length `400` from the Gaussian model
and two changes in the sequence.

    library(dust)
    set.seed(2)
    data <- dataGenerator_1D(chpts = c(110, 260, 400), c(0,0.8,0), type = "gauss")

We segment data using the dust 1D method coded in Rcpp.

    dust.1D(data)

    ## $changepoints
    ## [1] 105 261 400
    ## 
    ## $lastIndexSet
    ##  [1] 400 399 397 388 261 260 105 104 103  96
    ## 
    ## $nb
    ##   [1]  1  1  2  2  1  1  1  1  2  2  2  2  2  3  3  3  3  4  3  2  3  2  2  3  3
    ##  [26]  3  4  4  4  4  4  2  2  2  3  3  4  4  4  4  4  5  5  4  4  4  4  4  4  5
    ##  [51]  5  5  6  5  4  4  4  3  4  4  4  3  4  4  5  4  5  5  5  5  6  5  6  6  7
    ##  [76]  7  6  5  5  4  4  3  3  4  5  4  5  5  5  6  6  5  5  4  4  4  3  3  4  4
    ## [101]  5  5  6  7  7  6  5  6  6  6  6  7  7  8  8  8  9  9  9  9  9  9 10  9  8
    ## [126]  8  8  8  8  9 10  9  9  8  8  8  8  8  9  9  9 10 10  9 10 10 10 10  9 10
    ## [151]  9 10  9  9  8  8  9  8  8  9  9  8  8  9 10 10  9  9  9 10 11  9  9  9  9
    ## [176]  9  9 10 10 10 10 10 11 11 11 12 12 12 12 12 12 12 11 11 10  9 10 11 10 10
    ## [201] 10 10  8  7  8  9  9  9 10 10 10  9  9  9 10 10 10  8  7  7  6  6  7  7  7
    ## [226]  7  8  8  7  8  7  7  6  6  7  6  6  6  6  6  5  6  6  5  5  6  7  8  8  7
    ## [251]  8  8  8  9  9  8  9  9  9 10 10  9  7  8  8  9  9  9  9  9  9  9  9 10  9
    ## [276] 10 11 11  9  8  7  8  8  9 10  9  9  9 10 10 10  9  9 10  8  8  8  9  9  9
    ## [301]  9  9 10 10 10 10 11 11 12 12 12 13 12 12 12 12 12 11 10 11 11 11 11 11 11
    ## [326] 10 11 12 11 12 12 11 12 12 12 10 10 10 10  9 10  9 10 11 12 13 13 14 12 11
    ## [351] 12 12 11 11 12 12 12 13 13 14 15 14 13 13 12 11 11 11 12 11 10 11 11 12 13
    ## [376] 12 11 11 12 12 12 12 10 10 10  9  9  9  9  8  9  8  8  8  9  9 10 10 10  9
    ## 
    ## $costQ
    ##   [1] -4.022279e-01 -1.267593e-01 -1.278318e-01 -8.102372e-03 -1.121228e-02
    ##   [6] -3.414730e-03 -1.825415e-02 -4.416584e-03 -2.813259e-01 -2.229251e-01
    ##  [11] -2.907585e-01 -5.136066e-01 -3.739740e-01 -1.542998e-01 -4.968553e-01
    ##  [16] -7.505057e-02 -1.734335e-01 -1.686645e-01 -3.181374e-01 -3.820502e-01
    ##  [21] -8.571542e-01 -5.236612e-01 -8.875855e-01 -1.450606e+00 -1.394230e+00
    ##  [26] -6.688866e-01 -7.525741e-01 -5.962307e-01 -7.443397e-01 -7.843620e-01
    ##  [31] -9.313906e-01 -9.796183e-01 -1.225696e+00 -1.115664e+00 -8.991228e-01
    ##  [36] -7.478083e-01 -4.255632e-01 -2.917919e-01 -2.208117e-01 -1.904747e-01
    ##  [41] -1.511024e-01 -2.900472e-02 -6.014370e-03 -7.816770e-02 -1.170170e-01
    ##  [46] -2.980141e-01 -2.586330e-01 -2.439990e-01 -2.211755e-01 -1.195013e-01
    ##  [51] -6.722658e-02 -2.110421e-01 -1.603429e-01 -2.698379e-01 -1.720865e-01
    ##  [56] -5.078464e-02 -3.729530e-02 -7.747343e-02 -1.450438e-01 -2.811718e-01
    ##  [61] -1.324903e-01 -2.953442e-01 -2.270373e-01 -2.369032e-01 -2.781173e-01
    ##  [66] -2.042917e-01 -7.613547e-02 -5.419173e-02 -5.676973e-02 -2.588030e-02
    ##  [71] -6.793863e-03 -1.196573e-02 -9.391976e-03 -1.742385e-02 -1.606051e-02
    ##  [76] -2.737087e-03 -2.465413e-02 -4.743650e-02 -9.009050e-02 -3.489261e-02
    ##  [81] -6.963843e-02 -1.686360e-02 -7.687353e-03 -3.503929e-04 -3.532448e-02
    ##  [86] -2.296002e-03 -9.442798e-03 -1.394265e-02 -2.143790e-02 -1.363725e-02
    ##  [91] -6.218275e-06 -1.598112e-02 -1.517003e-03 -3.640274e-03 -2.881743e-02
    ##  [96] -6.723937e-02 -1.375737e-02 -1.348983e-02 -3.077909e-02 -4.711884e-02
    ## [101] -1.971013e-02 -1.475189e-02 -2.038119e-02 -3.765625e-02 -6.381851e-02
    ## [106] -1.226973e-02 -2.115852e-03 -8.260794e-03 -3.835574e-03 -1.443742e-03
    ## [111] -5.091466e-04 -3.501923e-03 -2.059726e-02 -8.171931e-02 -1.403034e-01
    ## [116] -2.074001e-01 -3.436647e-01 -5.047928e-01 -5.871578e-01 -5.839094e-01
    ## [121] -7.960009e-01 -8.998970e-01 -1.220665e+00 -1.262988e+00 -1.218543e+00
    ## [126] -1.401257e+00 -1.409587e+00 -1.619355e+00 -1.464274e+00 -1.315521e+00
    ## [131] -1.524276e+00 -1.688795e+00 -1.951892e+00 -1.736029e+00 -2.059734e+00
    ## [136] -2.406187e+00 -3.170026e+00 -4.244818e+00 -6.689517e+00 -5.499448e+00
    ## [141] -5.226633e+00 -5.931769e+00 -5.455332e+00 -5.844996e+00 -6.869700e+00
    ## [146] -6.550065e+00 -7.485719e+00 -8.323082e+00 -7.856652e+00 -7.234646e+00
    ## [151] -8.455441e+00 -8.607038e+00 -9.596529e+00 -9.116256e+00 -1.062502e+01
    ## [156] -9.680241e+00 -1.070138e+01 -1.144659e+01 -1.200873e+01 -1.294965e+01
    ## [161] -1.363772e+01 -1.334158e+01 -1.340098e+01 -1.404865e+01 -1.314488e+01
    ## [166] -1.265856e+01 -1.486995e+01 -1.326953e+01 -1.249146e+01 -1.368028e+01
    ## [171] -1.497024e+01 -1.604714e+01 -1.641266e+01 -1.702226e+01 -1.651440e+01
    ## [176] -1.624543e+01 -1.632841e+01 -1.582393e+01 -1.686106e+01 -1.576459e+01
    ## [181] -1.520532e+01 -1.446924e+01 -1.511439e+01 -1.450043e+01 -1.526703e+01
    ## [186] -1.656493e+01 -1.690534e+01 -1.765790e+01 -1.742761e+01 -1.816925e+01
    ## [191] -1.742806e+01 -1.768436e+01 -1.690489e+01 -1.700973e+01 -1.821308e+01
    ## [196] -1.910874e+01 -1.906397e+01 -1.874014e+01 -1.836839e+01 -1.809215e+01
    ## [201] -1.864947e+01 -1.816118e+01 -2.078445e+01 -2.128319e+01 -2.081626e+01
    ## [206] -2.144749e+01 -2.166842e+01 -2.170451e+01 -2.251053e+01 -2.301885e+01
    ## [211] -2.418073e+01 -2.406899e+01 -2.587658e+01 -2.500596e+01 -2.543991e+01
    ## [216] -2.613253e+01 -2.746428e+01 -2.668102e+01 -2.606467e+01 -2.776438e+01
    ## [221] -2.841525e+01 -2.933334e+01 -2.898941e+01 -3.092587e+01 -3.130447e+01
    ## [226] -3.202231e+01 -3.226446e+01 -3.334703e+01 -3.334549e+01 -3.561298e+01
    ## [231] -3.731119e+01 -3.591642e+01 -3.569843e+01 -3.599067e+01 -3.483718e+01
    ## [236] -3.426718e+01 -3.410860e+01 -3.483186e+01 -3.523156e+01 -3.505912e+01
    ## [241] -3.531036e+01 -3.806675e+01 -3.939168e+01 -4.018457e+01 -4.117881e+01
    ## [246] -4.100553e+01 -4.141253e+01 -4.197446e+01 -4.188126e+01 -4.180422e+01
    ## [251] -4.255477e+01 -4.252352e+01 -4.268920e+01 -4.273344e+01 -4.377632e+01
    ## [256] -4.451179e+01 -4.456762e+01 -4.464030e+01 -4.580429e+01 -4.597397e+01
    ## [261] -4.676196e+01 -4.630346e+01 -4.565974e+01 -4.368476e+01 -4.280080e+01
    ## [266] -4.131680e+01 -4.256045e+01 -4.260310e+01 -4.270026e+01 -4.233948e+01
    ## [271] -3.984221e+01 -3.949527e+01 -3.919293e+01 -3.794479e+01 -3.784328e+01
    ## [276] -3.878083e+01 -3.845887e+01 -3.895406e+01 -3.835712e+01 -3.812067e+01
    ## [281] -3.937891e+01 -3.958339e+01 -4.026913e+01 -3.950032e+01 -4.023754e+01
    ## [286] -4.122335e+01 -4.182886e+01 -4.153988e+01 -4.224644e+01 -4.188271e+01
    ## [291] -4.074066e+01 -4.089958e+01 -4.068531e+01 -4.023377e+01 -4.106454e+01
    ## [296] -4.052442e+01 -4.056035e+01 -3.875710e+01 -3.751270e+01 -3.677429e+01
    ## [301] -3.630245e+01 -3.583774e+01 -3.621117e+01 -3.478768e+01 -3.493712e+01
    ## [306] -3.524401e+01 -3.514345e+01 -3.553942e+01 -3.615529e+01 -3.545655e+01
    ## [311] -3.557587e+01 -3.503055e+01 -3.490408e+01 -3.480417e+01 -3.478949e+01
    ## [316] -3.480054e+01 -3.484682e+01 -3.485144e+01 -3.485464e+01 -3.484363e+01
    ## [321] -3.482168e+01 -3.487112e+01 -3.492632e+01 -3.494976e+01 -3.496916e+01
    ## [326] -3.500862e+01 -3.509244e+01 -3.529818e+01 -3.549147e+01 -3.543737e+01
    ## [331] -3.556347e+01 -3.553953e+01 -3.550127e+01 -3.559690e+01 -3.574485e+01
    ## [336] -3.568363e+01 -3.579296e+01 -3.577633e+01 -3.601496e+01 -3.607927e+01
    ## [341] -3.609456e+01 -3.614861e+01 -3.590791e+01 -3.593258e+01 -3.604057e+01
    ## [346] -3.598304e+01 -3.591656e+01 -3.606716e+01 -3.619732e+01 -3.625017e+01
    ## [351] -3.650732e+01 -3.661994e+01 -3.652706e+01 -3.681972e+01 -3.745469e+01
    ## [356] -3.813789e+01 -3.736211e+01 -3.767457e+01 -3.793871e+01 -3.782012e+01
    ## [361] -3.859822e+01 -3.929625e+01 -3.887237e+01 -3.836563e+01 -3.759426e+01
    ## [366] -3.731099e+01 -3.718873e+01 -3.742049e+01 -3.732602e+01 -3.731163e+01
    ## [371] -3.725095e+01 -3.773219e+01 -3.752499e+01 -3.777749e+01 -3.774333e+01
    ## [376] -3.756970e+01 -3.749140e+01 -3.803875e+01 -3.804824e+01 -3.819868e+01
    ## [381] -3.806926e+01 -3.787276e+01 -3.787970e+01 -3.769667e+01 -3.796457e+01
    ## [386] -3.759501e+01 -3.758239e+01 -3.755434e+01 -3.783358e+01 -3.799231e+01
    ## [391] -3.800834e+01 -3.828683e+01 -3.854187e+01 -3.872304e+01 -3.866551e+01
    ## [396] -3.871912e+01 -3.845738e+01 -3.869554e+01 -3.872190e+01 -3.891457e+01

Here the penalty value is by default set to `2 log(n)` for `n` data
points and the model to `gauss`. That is, we did
`dust.1D(data,  penalty = 2*log(length(data)), model = "gauss")`

The `penalty` value and the type of `model` can be explicitly given. It
can be one of the following: `"gauss"` (additional parameters `sdNoise`
and `gamma`), `"exp"`, `"poisson"`, `"geom"`, `"bern"`, `"binom"`
(additional parameter `nbTrials`), `"negbin"` (additional parameter
`nbSuccess`) and `variance`. See next [Section](#Models).

The result is a list whose elements are:

-   `changepoints`: A vector of change points (the index ending each of
    the segments)

-   `lastIndexSet`: The list of non-pruned indices at the end of the
    analysis

-   `nb`: The number of indices to consider at each time step (its
    length is equal to data length)

-   `costQ` The minimal (penalized) cost of the optimization problem at
    each time step.

Vector `nb` is a kind of complexity control vector, its values are
directly related to the algorithm’s time complexity.

And we can also use multivariate data with Poisson model.

    library(dust)
    set.seed(5)
    data <- dataGenerator_MD(chpts = c(50,100,150),
                     parameters = data.frame(ts1 = c(1,4,5),ts2 = c(10,10,5)),
                     type = "poisson")

    dust.MD(data, model = "poisson", method = "randIndex_Eval0")

    ## nb_l = 2; nb_r = 0; nb_max = 2

    ## $changepoints
    ## [1]  50 100 150
    ## 
    ## $lastIndexSet
    ##  [1] 150 149 148 147 146 145 144 143 142 141 140 125 124 118 117 115 109 106 102
    ## [20] 101 100  99  98  97  96  95  51  50  49  48  46  44  43  42  40  39  38   0
    ## 
    ## $nb
    ##   [1]  1  1  1  2  3  1  2  2  3  3  4  3  4  3  2  3  2  3  3  4  5  6  6  4  5
    ##  [26]  6  6  4  5  5  6  2  3  4  5  6  6  7  7  8  8  9  9 10 11 11 12 13 14 14
    ##  [51] 14 15 16 17 18 19 20 21 22 21 22 21 19 19 20 21 21 21 22 23 24 25 26 25 25
    ##  [76] 23 24 24 24 24 23 23 23 23 24 25 25 26 26 26 24 23 24 24 24 25 26 27 27 24
    ## [101] 22 22 23 24 25 26 27 28 29 30 30 31 32 33 34 33 33 34 35 36 37 37 38 39 40
    ## [126] 41 41 41 41 41 42 43 43 43 43 44 45 46 44 41 42 43 43 43 41 40 40 37 36 37
    ## 
    ## $costQ
    ##   [1]   -10.77502   -29.06055   -40.74796   -52.89471   -60.59678   -74.84412
    ##   [7]   -73.16865   -87.16346   -92.44774   -97.75021   -96.71445  -106.47458
    ##  [13]  -114.21088  -130.44482  -133.89960  -135.33146  -153.62380  -167.47645
    ##  [19]  -175.06350  -182.52172  -192.21256  -201.98204  -198.90313  -213.05924
    ##  [25]  -225.12154  -230.56859  -246.79305  -256.32434  -266.18251  -275.80940
    ##  [31]  -281.25111  -293.23161  -294.44123  -295.72423  -309.90580  -321.81106
    ##  [37]  -331.56732  -349.97061  -366.25069  -382.56066  -400.97767  -419.67107
    ##  [43]  -422.71653  -439.09756  -453.14865  -465.13452  -468.24446  -478.05064
    ##  [49]  -492.21029  -497.56775  -522.83606  -541.64737  -552.15700  -561.24443
    ##  [55]  -576.02167  -586.19067  -587.26390  -609.47764  -612.03583  -619.90318
    ##  [61]  -637.85130  -645.40156  -669.29298  -686.57240  -715.69918  -730.31811
    ##  [67]  -747.99431  -738.79518  -753.53255  -757.81337  -771.03823  -780.45523
    ##  [73]  -804.95970  -813.70613  -818.97365  -827.14159  -834.67578  -841.92947
    ##  [79]  -854.18800  -869.27854  -897.18349  -905.95271  -930.56296  -950.89430
    ##  [85]  -976.62823  -986.60771 -1000.43764 -1011.98755 -1029.47442 -1035.01262
    ##  [91] -1044.33126 -1052.86580 -1066.61260 -1077.41043 -1099.41666 -1121.42077
    ##  [97] -1132.99242 -1140.78386 -1146.34879 -1158.64544 -1172.39084 -1174.30599
    ## [103] -1174.81205 -1179.91474 -1187.98617 -1189.35682 -1195.11493 -1202.40702
    ## [109] -1209.71448 -1206.90897 -1213.37390 -1220.73339 -1226.87809 -1234.87926
    ## [115] -1244.40627 -1250.29197 -1261.94863 -1260.58079 -1268.01394 -1266.11972
    ## [121] -1270.56251 -1279.33399 -1285.78048 -1297.93767 -1306.00214 -1302.76391
    ## [127] -1308.48811 -1319.48942 -1322.29233 -1328.37791 -1337.75907 -1343.92269
    ## [133] -1345.28200 -1356.21503 -1362.37363 -1357.23921 -1358.65262 -1364.84369
    ## [139] -1369.37821 -1373.79921 -1386.65669 -1384.65586 -1387.54212 -1391.64379
    ## [145] -1394.38591 -1397.17713 -1403.26461 -1417.38362 -1418.56742 -1424.62087

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="Models"></a>

## Models And Data Generators

### Data Generators in 1D and MultiD

**dataGenerator\_1D** is used to generate data with a given vector of
change-point (e.g. `chpts = c(50,100)` for one change at position `50`
and data length `100`), parameter vector (e.g. `parameters = c(0,1)`)
and a type of probability distribution (from the exponential family) in
`type`. The following types are available in the current package
version:

-   `type = "gauss"` (additional parameters `sdNoise` and `gamma`)

-   `type = "exp"`

-   `type = "poisson"`

-   `type = "geom"`

-   `type = "bern"`

-   `type = "binom"` (additional parameter `nbTrials`)

-   `type = "negbin"` (additional parameter `nbSuccess`)

-   `type = "variance"`

We show two data examples with Gaussian and Exponential models
(`"gauss"` and `"exp"`)

<center>
<img src="man/figures/cost1.png" alt="" style="width:80%;"/>
</center>
and some other examples with integer-valued cost (`"poisson"`,`"geom"`,
`"binom"`, `"negbin"`):
<center>
<img src="man/figures/cost2.png" alt="" style="width:80%;"/>
</center>

**dataGenerator\_MD** concatenates `p` copies of `dataGenerator_1D`
function.

Additional information and examples are easily accessible in the help of
these functions (e.g. run `?dataGenerator_MD`).

### Data Generators in 2D

**dataGenerator\_meanVar** is used for change in mean and variance for
the Gaussian problem

**dataGenerator\_Reg** generates 2-dimensional data frame `(x,y)`
following a simple linear regression link (`y = Ax + B + noise`) with
`A` and `B` changing over time (after each change-point)

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="dust1D"></a>

## dust 1D and MD Algorithms

We start with one simple example with the `exp` model:

    data <- dataGenerator_1D(chpts = c(5*1e5,1e6), parameters = c(2,1), type = "exp")
    system.time(res <- dust.1D(data = data, model = "exp"))[[1]]

    ## [1] 0.821

    res$changepoints

    ## [1]  499991 1000000

A fundamental information relies in the number of indices to consider at
each data step. It is saved into the field `nb`.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="dust2D"></a>

## dust 2D Algorithms

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="rcpp"></a>

## Rcpp Structure and new cost integration

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="pruning"></a>

## Pruning Capacity

Analysis of the pruning capacity (return field `nb`) for some of our
algorithms. We explore in particular the impact of the different choices
for the dual max evaluation.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

## Hidden Functions and Parameters (Package Development)

### OP in R

The base function `OP_R` is used to compute the change-point vector with
the simplest dynamic programming algorithm with no pruning. This method
is of quadratic time complexity. We propose 3 such functions:`OP_R_1D`,
`OP_R_MD`, `OP_R_2param`.

`OP_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss")`

Example:
`OP_R_1D(dataGenerator_1D(chpts = c(200,400), c(0,1), type = "gauss"))`

`OP_R_2param` is used for:

-   `type = "meanVar"` change in Gaussian data in mean and variance

-   `type = regression` change in simple regression model

### Dual Functions

`dual_1D` returns the value of the dual at a point `mu` when comparing
index `s1` with the constraint from index `s2` at time `t`. With option
`OP = TRUE` the optimal partitioning algorithm is used to have the true
constants in the cost functions with penalty `penalty` and a pruning
option `pruningOpt`.

`dual_1D <- function(mu, data, s1, s2, t, type = "gauss", OP = FALSE, penalty = 2*length(data), pruningOpt = 3)`

-   `data` is raw data

-   If `OP` is `true`, we run the OP algorithm to have the optimal cost
    vector to use in cost functions. See the function `OP_R`.

-   at time `t`, we evaluate the dual function at point `mu` when trying
    to remove index `s1` using function linked to index `s2` (we have a
    unique constraint, which means that the dual is a one-parametric
    function)

-   Depending on the `type`, different functions `A`, `B`, `statistic`,
    `mu_max` and `evalDual` are used (see the code in file
    `functions_by_type.R`)

Function `dual_1D` allows us to study the shape of the dual.

### dust\_R

We propose a few R functions computing the change-point location with
dust method: `dust_R_1D`, `dust_R_MD`, `dust_R_2Dquad`.

The function `dust_R_1D` has the following parameters:

`dust_R_1D <- function(data, penalty = 2*log(length(data)), type = "gauss", pruningOpt = 2)`

We have different type of possible pruning:

-   `pruningOpt == 0`: nothing

-   `pruningOpt == 1`: PELT

-   `pruningOpt == 2`: dust

-   `pruningOpt == 3`: dust + PELT

and returns a list of two elements:

-   `changepoints`: the change-points found by the algo

-   `nb`: the number of indices to consider in the minimization at each
    time step

-   `lastIndexSet`: the vecotr of indices saved by the algo in the
    dynamic programming algorithm at the last iteration

-   `costQ`: the vector of optimal cost (of size `length(data)`)

### Plot functions

`plot_dual_1D` is a function displaying the result of `dual_1D` for a
vector of mu values.

`plot_dual_1D <- function(mu =  (1:99)/100,` `data, s1, s2,`
`type = "gauss",` `OP = FALSE,` `penalty = 2*length(data))`

We use the function `plot_dual_1D` with `OP = TRUE` to plot the true
dual function seen by the dynamic programming algorithm.

What we called the “pruning interval” is the interval of values between
the vertical green lines for which the dual function takes a value
higher than the pruning threshold (horizontal line in red), so that the
index considered `s1` is pruned by `s2` at time `n`.

Using function `barplot_dual_1D` we can repeat the generation of the
pruning interval `nb` and count the number of time each value mu is in
this interval.

We add the values in the bar plot only if at the final time step `n`,
the index `s1` has not been removed by the algorithm (the pruning option
is given by option `pruningOpt`).

`barplot_dual_1D <- function(nb = 1000, s1 = 18, s2 = 15,` `n = 20,`
`oneParam = 0,` `type = "gauss",` `penalty = 2*log(n),`
`pruningOpt = 3)`

[(Back to Top)](#top)

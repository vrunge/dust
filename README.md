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
    data <- dataGenerator_1D(chpts = c(110, 260, 400), c(0,0.4,0), type = "gauss")

We segment data using the dust 1D method coded in Rcpp.

    dust.1D(data)

    ## $changepoints
    ## [1] 400
    ## 
    ## $lastIndexSet
    ##  [1] 400 399 397 394 388 362 105 104 103  96   0
    ## 
    ## $nb
    ##   [1]  1  1  2  2  1  1  1  1  2  2  2  2  2  3  3  3  3  4  3  2  3  2  2  3  3
    ##  [26]  3  4  4  4  4  4  2  2  2  3  3  4  4  4  4  4  5  5  4  4  4  4  4  4  5
    ##  [51]  5  5  6  5  4  4  4  3  4  4  4  3  4  4  5  4  5  5  5  5  6  5  6  6  7
    ##  [76]  7  6  5  5  4  4  3  3  4  5  4  5  5  5  6  6  5  5  4  4  4  3  3  4  4
    ## [101]  5  5  6  7  7  6  5  6  6  6  7  7  8  8  8  8  9  9  9  9  9  9 10 10  9
    ## [126]  9 10  9  9  9  9  9  9  9  9  7  7  7  8  8  8  9 10 11 10  9  8  8  8  8
    ## [151]  9  8  8  7  7  7  8  7  7  8  8  8  8  8  7  7  6  6  6  7  8  6  6  6  7
    ## [176]  7  7  7  8  7  7  7  8  8  8  9  8  8  8  8  9 10 10 10 11  9  9 10 11 10
    ## [201] 10 10  9  7  8  9  9 10  9  8  9  8  9  9  9 10 10  8  8  9  8  8  9  9  9
    ## [226]  9  9 10  9 10 10 10  9 10 10 10  9  9  9 10 11 12 12 11 11 12 13 13 12 12
    ## [251] 13 13 13 14 13 12 12 11 11 12 12 12 12 12 11 10 11 10 10 11  9  9  9 10 10
    ## [276] 11 12 11 10 10  9  9  9 10 10 11 11 11 11 11 11 12 13 13 13 13 14 13 13 13
    ## [301] 13 13 14 13 12 13 12 11 12 11 12 12 13 13 13 14 13 14 14 15 15 12 11 11 11
    ## [326] 11 12 12 12 13 13 12 13 14 13 11 11 12 11 11 11 11 12 12 12 12 13 13 13 13
    ## [351] 12 12 12 11 11 11 11 12 12 13 13 12 12 12 12 11 10 11 12 10 10 10 10 11 12
    ## [376] 12 10 10 10  9 10 11 10 11 12 10  9  9  9  8  9  8  8  8  9  9 10 10 10 10
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
    ## [111] -1.833523e-05 -3.277338e-05 -4.057008e-03 -3.236520e-02 -5.890077e-02
    ## [116] -8.871122e-02 -1.625602e-01 -2.521902e-01 -2.839919e-01 -2.559761e-01
    ## [121] -3.713027e-01 -4.113177e-01 -5.979884e-01 -5.901726e-01 -5.247600e-01
    ## [126] -6.093126e-01 -5.784975e-01 -6.765699e-01 -5.430457e-01 -4.235701e-01
    ## [131] -5.121723e-01 -5.744637e-01 -6.939110e-01 -5.346114e-01 -6.832619e-01
    ## [136] -8.474975e-01 -1.007503e+00 -1.148144e+00 -1.484006e+00 -1.324483e+00
    ## [141] -1.290201e+00 -1.392192e+00 -1.326444e+00 -1.383726e+00 -1.536522e+00
    ## [146] -1.489542e+00 -1.633185e+00 -1.765839e+00 -1.692135e+00 -1.592613e+00
    ## [151] -1.787504e+00 -1.811164e+00 -1.975546e+00 -1.892405e+00 -2.150524e+00
    ## [156] -1.983493e+00 -2.159833e+00 -2.290967e+00 -2.391001e+00 -2.563594e+00
    ## [161] -2.691697e+00 -2.631992e+00 -2.639599e+00 -2.761182e+00 -2.580476e+00
    ## [166] -2.480558e+00 -2.909743e+00 -2.585703e+00 -2.426239e+00 -2.651875e+00
    ## [171] -2.904253e+00 -3.119017e+00 -3.188093e+00 -3.309448e+00 -3.193398e+00
    ## [176] -3.127531e+00 -3.135698e+00 -3.018809e+00 -3.229757e+00 -2.984703e+00
    ## [181] -2.855143e+00 -2.689243e+00 -2.811700e+00 -2.670351e+00 -2.818309e+00
    ## [186] -3.083559e+00 -3.144568e+00 -3.297002e+00 -3.232390e+00 -3.383772e+00
    ## [191] -3.204315e+00 -3.246797e+00 -3.059342e+00 -3.067298e+00 -3.320411e+00
    ## [196] -3.508624e+00 -3.482769e+00 -3.392873e+00 -3.292211e+00 -3.213352e+00
    ## [201] -3.323074e+00 -3.195472e+00 -3.785017e+00 -3.886293e+00 -3.758070e+00
    ## [206] -3.890678e+00 -3.925858e+00 -3.916575e+00 -4.093253e+00 -4.199053e+00
    ## [211] -4.467838e+00 -4.421311e+00 -4.858304e+00 -4.618246e+00 -4.708833e+00
    ## [216] -4.866091e+00 -5.190855e+00 -4.967722e+00 -4.789109e+00 -5.209243e+00
    ## [221] -5.360017e+00 -5.583254e+00 -5.470775e+00 -5.970857e+00 -6.052950e+00
    ## [226] -6.228582e+00 -6.273850e+00 -6.552938e+00 -6.530599e+00 -7.149655e+00
    ## [231] -7.616976e+00 -7.191917e+00 -7.106777e+00 -7.167707e+00 -6.815187e+00
    ## [236] -6.630505e+00 -6.562306e+00 -6.743223e+00 -6.833081e+00 -6.759954e+00
    ## [241] -6.807329e+00 -7.579907e+00 -7.945843e+00 -8.156983e+00 -8.430168e+00
    ## [246] -8.352811e+00 -8.450043e+00 -8.594456e+00 -8.540435e+00 -8.491149e+00
    ## [251] -8.692974e+00 -8.657133e+00 -8.681134e+00 -8.667963e+00 -8.960580e+00
    ## [256] -9.160447e+00 -9.150491e+00 -9.145593e+00 -9.480101e+00 -9.505375e+00
    ## [261] -9.834176e+00 -9.769761e+00 -9.645212e+00 -9.099022e+00 -8.896237e+00
    ## [266] -8.505219e+00 -8.966877e+00 -9.052994e+00 -9.156441e+00 -9.112993e+00
    ## [271] -8.388886e+00 -8.346366e+00 -8.316962e+00 -7.987997e+00 -8.019435e+00
    ## [276] -8.381301e+00 -8.342324e+00 -8.565957e+00 -8.436589e+00 -8.422944e+00
    ## [281] -8.897916e+00 -9.029986e+00 -9.323587e+00 -9.130766e+00 -9.443138e+00
    ## [286] -9.844072e+00 -1.011879e+01 -1.008529e+01 -1.039764e+01 -1.033695e+01
    ## [291] -1.000327e+01 -1.012308e+01 -1.011194e+01 -1.001664e+01 -1.037349e+01
    ## [296] -1.024507e+01 -1.032065e+01 -9.742309e+00 -9.363601e+00 -9.162444e+00
    ## [301] -9.053810e+00 -8.946870e+00 -9.133180e+00 -8.647579e+00 -8.794172e+00
    ## [306] -8.955674e+00 -8.973553e+00 -9.167799e+00 -9.442546e+00 -9.244980e+00
    ## [311] -9.341132e+00 -9.197155e+00 -9.046659e+00 -8.550960e+00 -8.391855e+00
    ## [316] -8.475028e+00 -8.731651e+00 -8.731751e+00 -8.725034e+00 -8.651056e+00
    ## [321] -8.508920e+00 -8.734455e+00 -8.922729e+00 -8.980819e+00 -9.022599e+00
    ## [326] -9.119844e+00 -9.321844e+00 -9.753318e+00 -1.009247e+01 -9.985034e+00
    ## [331] -1.018873e+01 -1.013637e+01 -1.005826e+01 -1.020764e+01 -1.043224e+01
    ## [336] -1.032490e+01 -1.048440e+01 -1.044805e+01 -1.079052e+01 -1.087170e+01
    ## [341] -1.088398e+01 -1.094993e+01 -1.059522e+01 -1.062188e+01 -1.077028e+01
    ## [346] -1.067660e+01 -1.056696e+01 -1.077999e+01 -1.095565e+01 -1.102020e+01
    ## [351] -1.135607e+01 -1.149349e+01 -1.136683e+01 -1.172991e+01 -1.230244e+01
    ## [356] -1.263802e+01 -1.235376e+01 -1.253883e+01 -1.269862e+01 -1.269605e+01
    ## [361] -1.307616e+01 -1.342485e+01 -1.329061e+01 -1.312012e+01 -1.258432e+01
    ## [366] -1.225327e+01 -1.210231e+01 -1.237180e+01 -1.225561e+01 -1.223320e+01
    ## [371] -1.215488e+01 -1.271249e+01 -1.246975e+01 -1.275580e+01 -1.271281e+01
    ## [376] -1.250755e+01 -1.241025e+01 -1.303415e+01 -1.304120e+01 -1.320439e+01
    ## [381] -1.305752e+01 -1.283114e+01 -1.283480e+01 -1.261742e+01 -1.292350e+01
    ## [386] -1.248639e+01 -1.246584e+01 -1.242648e+01 -1.275383e+01 -1.293380e+01
    ## [391] -1.294787e+01 -1.326098e+01 -1.354102e+01 -1.373590e+01 -1.367048e+01
    ## [396] -1.372607e+01 -1.343391e+01 -1.369439e+01 -1.372027e+01 -1.392731e+01

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
    ##  [1] 150 149 148 147 146 144 143 142 141 139 138 137 135 133 132 130 125 118 117
    ## [20] 113 109 108 105 104 101 100  98  97  96  52  51  50  49  48  47  46  45  44
    ## [39]  43  42  40  38  37  36  34   0
    ## 
    ## $nb
    ##   [1]  1  1  2  2  3  1  2  2  3  4  4  4  5  6  5  6  4  3  4  4  5  6  3  4  5
    ##  [26]  6  5  4  4  4  5  3  3  4  5  6  5  6  6  6  7  8  9  9 10 11 12 13 14 15
    ##  [51] 16 17 18 19 20 21 22 21 22 22 22 23 21 21 22 23 23 24 25 26 27 28 27 27 28
    ##  [76] 28 29 29 30 29 27 24 24 25 26 27 27 27 27 27 26 27 27 27 27 28 29 29 29 28
    ## [101] 26 24 24 25 26 27 28 29 30 31 32 33 34 35 35 35 36 37 38 39 40 41 42 43 44
    ## [126] 44 45 46 46 46 45 45 46 45 46 47 48 48 48 49 50 50 50 50 48 48 49 45 44 45
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

    ## [1] 0.841

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

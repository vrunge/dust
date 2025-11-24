<a id="top"></a>

# dust Vignette

### Vincent Runge and Simon Querné

#### LaMME, Evry University, November 24, 2025

<center>
<img src="man/figures/dust.png" alt="" style="width:30%;"/>
</center>

> [Quick start](#start)

> [Rcpp Object Structure](#rcpp)

> [dust Algorithms](#dust1D)

> [Models And Data Generators](#Models)

> [New cost integration](#cost)

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
    set.seed(11)
    data <- dataGenerator_1D(chpts = c(40, 160, 200), c(0,1.3,0), type = "gauss")

We segment data using the dust 1D method coded in Rcpp.

    dust.1D(data)

    ## $changepoints
    ## [1]  41 158 200
    ## 
    ## $lastIndexSet
    ## [1] 200 199 197 195 194 160 158
    ## 
    ## $nb
    ##   [1] 1 1 1 1 2 2 2 2 2 2 1 1 2 2 2 2 1 2 2 2 2 2 2 2 2 2 1 2 2 2 3 3 2 2 2 2 2
    ##  [38] 3 4 4 3 3 4 3 4 4 4 4 4 4 4 4 5 5 4 4 4 4 5 6 6 5 5 5 6 5 4 4 4 4 3 3 4 4
    ##  [75] 4 4 3 3 4 4 4 3 3 4 4 4 4 4 4 5 4 5 6 7 7 7 5 4 4 4 3 4 4 3 4 5 5 6 5 4 4
    ## [112] 5 4 5 4 4 3 4 5 5 5 6 6 6 6 6 6 6 7 7 6 7 7 7 6 5 4 5 5 6 6 6 6 6 5 6 5 6
    ## [149] 7 8 7 7 7 7 8 8 7 6 6 6 6 6 4 4 5 5 6 6 5 5 6 6 7 7 6 5 5 5 5 6 6 6 6 6 7
    ## [186] 6 7 5 4 4 4 3 3 4 3 4 4 5 5 6
    ## 
    ## $costQ
    ##   [1]   -0.17465888   -0.07964721   -0.72175311   -1.48233479   -0.51309228
    ##   [6]   -0.85296289   -0.25130353   -0.09777846   -0.09338471   -0.26464387
    ##  [11]   -0.44504583   -0.50384930   -0.96758787   -0.99236705   -1.37439967
    ##  [16]   -1.28355731   -1.29357419   -0.91651424   -1.05652476   -1.22218706
    ##  [21]   -1.40232339   -1.34411910   -1.43792819   -1.26113201   -1.18802200
    ##  [26]   -1.14020774   -1.15213070   -1.33715207   -1.35784737   -1.61966884
    ##  [31]   -1.93825216   -2.21281710   -1.90801487   -2.40915056   -2.66926069
    ##  [36]   -2.41434475   -2.41585792   -1.83989705   -1.98288693   -2.04294780
    ##  [41]   -2.09942991   -1.65771373   -1.07354750   -0.95792134   -0.54309089
    ##  [46]   -0.39590486   -0.50971677   -0.23422424   -0.07538479   -0.01435639
    ##  [51]   -1.32017343   -1.84275887   -1.59593395   -3.09880205   -4.15633867
    ##  [56]   -5.73817538   -6.36606355   -7.81132599  -10.76828580  -11.68829310
    ##  [61]  -11.30862198  -15.55923953  -16.20051075  -14.20396735  -15.80426188
    ##  [66]  -19.10100166  -18.77178732  -18.07619799  -18.07207083  -19.96122688
    ##  [71]  -21.43849869  -22.10367325  -21.88271719  -23.03956617  -25.44463704
    ##  [76]  -26.50451772  -26.37785265  -27.85086682  -27.21230480  -28.60670837
    ##  [81]  -29.35946615  -30.63635446  -30.64680592  -30.26156952  -34.19818126
    ##  [86]  -34.20388668  -33.27848332  -34.80713531  -34.50074525  -33.35352378
    ##  [91]  -35.78420616  -35.27174772  -37.35423129  -37.48062161  -39.00359890
    ##  [96]  -39.38188636  -41.99229902  -41.31762575  -44.04730362  -44.08423472
    ## [101]  -44.34085353  -45.50920963  -47.15044669  -47.83402563  -45.94520735
    ## [106]  -47.55665693  -49.04599933  -48.56944633  -51.08330949  -53.42460377
    ## [111]  -55.40210602  -58.45288164  -59.91094326  -57.62256425  -58.36006770
    ## [116]  -58.50311637  -58.79063716  -60.04135191  -60.84393218  -61.82858076
    ## [121]  -65.99562618  -65.96822750  -68.46240795  -70.28176591  -70.84504518
    ## [126]  -73.63745906  -76.27841141  -77.94139201  -80.59105686  -80.19058919
    ## [131]  -80.83402567  -79.99846036  -81.29868493  -80.63410088  -80.76123546
    ## [136]  -81.49561442  -81.85789768  -79.62989199  -80.08807096  -81.30964616
    ## [141]  -80.07036195  -81.64044839  -80.44150161  -80.90055295  -81.59897528
    ## [146]  -82.78446220  -84.14290916  -87.67096558  -88.08180099  -89.43304095
    ## [151]  -88.15010160  -89.88622801  -91.29824453  -93.06509976  -92.58564731
    ## [156]  -94.55392416  -97.68270697 -100.68279475 -100.48154709 -100.94213224
    ## [161] -100.03822185  -98.38764312  -98.28454299  -97.55628227  -95.00366986
    ## [166]  -94.05052730  -91.44183716  -90.88781895  -90.81742064  -90.78049875
    ## [171]  -90.61473438  -90.52948916  -90.39709647  -90.37109600  -90.47043669
    ## [176]  -90.36934139  -90.44770671  -90.35968250  -90.34623420  -90.38967439
    ## [181]  -90.38162422  -90.37495952  -90.41657289  -90.36995176  -90.57122809
    ## [186]  -90.52955716  -90.48124056  -90.46221764  -90.41741677  -90.44035115
    ## [191]  -90.43497986  -90.39941291  -90.45834265  -90.46106642  -90.49347374
    ## [196]  -90.84380211  -91.18447544  -90.89627946  -90.91145973  -91.08919884

Here the penalty value is by default set to `2 log(n)` for `n` data
points and the model to `gauss`. That is, we did
`dust.1D(data,  penalty = 2*log(length(data)), model = "gauss")`

The `penalty` value and the type of `model` can be explicitly given. It
can be one of the following: `"gauss"` (additional parameters `sdNoise`
and `gamma`), `"exp"`, `"poisson"`, `"geom"`, `"bern"`, `"binom"`
(additional parameter `nbTrials`), `"negbin"` (additional parameter
`nbSuccess`) and `variance`. See next [Section](#Models).

*The result is a list whose elements are:*

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
    ##  [1] 150 149 148 147 146 145 144 143 142 141 140 135 134 119 118 117 116 115 110
    ## [20] 107 105 101 100  98  97  96  52  51  50  49  48  46  44  43  42  40  39  37
    ## [39]   0
    ## 
    ## $nb
    ##   [1]  1  1  1  2  3  1  2  3  3  2  3  3  4  3  4  5  4  3  3  2  3  1  2  3  4
    ##  [26]  5  6  4  4  4  5  5  4  5  6  7  8  9  9  9  9 10  9 10 11 12 13 14 15 16
    ##  [51] 15 16 17 18 19 19 20 21 22 23 24 22 21 22 23 24 24 25 25 26 27 27 24 25 25
    ##  [76] 24 24 24 24 22 22 23 22 22 23 24 25 25 24 24 23 24 25 26 26 25 26 27 27 26
    ## [101] 25 23 23 24 25 26 27 28 29 30 31 32 33 34 33 34 34 35 36 37 38 39 40 41 41
    ## [126] 42 42 42 42 39 38 38 39 39 39 39 40 41 42 41 41 40 40 41 41 41 42 42 42 38
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

<a id="rcpp"></a>

## Rcpp Object Structure

    obj_dust <- new(DUST_1D, "gauss", "randIndex_Eval3", 5)
    obj_dust$get_info()

    ## $data_statistic
    ## numeric(0)
    ## 
    ## $data_length
    ## [1] 0
    ## 
    ## $current_penalty
    ## [1] 0
    ## 
    ## $model
    ## [1] "gauss"
    ## 
    ## $pruning_algo
    ## [1] 3
    ## 
    ## $pruning_constraints_type
    ## [1] 0
    ## 
    ## $pruning_nb_loops
    ## [1] 5

    data <- dataGenerator_1D(chpts = c(500,1000), parameters = c(0,1), type = "gauss")
    obj_dust$append_c(data, 2*log(5000))

    obj_dust$update_partition()
    obj_dust$get_info()

    ## $data_statistic
    ##    [1]  0.000000e+00 -4.755915e-01 -1.359694e+00 -4.857753e+00 -5.239736e+00
    ##    [6] -4.262048e+00 -4.820089e+00 -5.446544e+00 -5.976995e+00 -4.079374e+00
    ##   [11] -2.683833e+00 -3.429859e+00 -3.735432e+00 -2.565754e+00 -2.261367e+00
    ##   [16] -2.378865e+00 -2.438951e+00 -9.680116e-01 -2.446159e+00 -3.129772e+00
    ##   [21] -2.669232e+00 -2.850733e+00 -4.009550e+00 -3.600531e+00 -3.858738e+00
    ##   [26] -4.125637e+00 -3.961481e+00 -4.354940e+00 -6.198678e+00 -7.740966e+00
    ##   [31] -8.327206e+00 -9.179345e+00 -8.401021e+00 -8.431338e+00 -9.886995e+00
    ##   [36] -9.793210e+00 -8.810861e+00 -9.407571e+00 -9.332766e+00 -7.135337e+00
    ##   [41] -6.340314e+00 -6.879256e+00 -8.480539e+00 -9.211913e+00 -9.567653e+00
    ##   [46] -1.055307e+01 -1.128424e+01 -9.818914e+00 -7.960299e+00 -7.956802e+00
    ##   [51] -9.300577e+00 -9.149283e+00 -8.859274e+00 -8.981752e+00 -8.856598e+00
    ##   [56] -9.629032e+00 -1.064200e+01 -9.675079e+00 -1.009841e+01 -1.093001e+01
    ##   [61] -9.530048e+00 -9.513034e+00 -7.665537e+00 -8.351877e+00 -8.570601e+00
    ##   [66] -7.888008e+00 -7.363878e+00 -7.283328e+00 -7.229542e+00 -7.966801e+00
    ##   [71] -7.000363e+00 -6.016056e+00 -5.828765e+00 -5.555780e+00 -4.345634e+00
    ##   [76] -4.156976e+00 -2.194478e+00 -2.055766e+00 -3.634393e+00 -4.431414e+00
    ##   [81] -3.207060e+00 -3.572394e+00 -3.734984e+00 -3.174505e+00 -4.035231e+00
    ##   [86] -2.796596e+00 -2.029161e+00 -3.116570e+00 -3.049066e+00 -1.443925e+00
    ##   [91] -2.116959e-01 -5.908097e-01 -1.940676e+00 -1.575758e+00 -1.939353e+00
    ##   [96] -5.646998e-01 -2.728040e-01  4.377877e-01 -4.999732e-01 -1.614036e+00
    ##  [101] -9.796644e-01 -1.210857e+00 -2.579052e+00 -3.333959e+00 -4.459556e+00
    ##  [106] -4.678915e+00 -4.813223e+00 -5.631243e+00 -5.158909e+00 -6.028165e+00
    ##  [111] -7.360454e+00 -7.289891e+00 -6.825798e+00 -6.536639e+00 -9.421580e+00
    ##  [116] -1.175627e+01 -1.348716e+01 -1.266215e+01 -1.370719e+01 -1.458439e+01
    ##  [121] -1.498478e+01 -1.625297e+01 -1.611438e+01 -1.493081e+01 -1.704136e+01
    ##  [126] -1.678069e+01 -1.583512e+01 -1.645508e+01 -1.646418e+01 -1.594396e+01
    ##  [131] -1.413770e+01 -1.605022e+01 -1.585093e+01 -1.557445e+01 -1.641073e+01
    ##  [136] -1.438113e+01 -1.395202e+01 -1.288807e+01 -1.349390e+01 -1.241952e+01
    ##  [141] -1.313187e+01 -1.331073e+01 -1.281097e+01 -1.301054e+01 -1.292065e+01
    ##  [146] -1.191574e+01 -1.378868e+01 -1.326390e+01 -1.377852e+01 -1.256759e+01
    ##  [151] -1.308039e+01 -1.198976e+01 -1.149792e+01 -1.174054e+01 -9.624497e+00
    ##  [156] -8.435858e+00 -7.540596e+00 -6.672056e+00 -5.015878e+00 -3.570242e+00
    ##  [161] -2.825641e+00 -3.515832e+00 -4.307243e+00 -4.569324e+00 -4.977242e+00
    ##  [166] -4.775931e+00 -5.507878e+00 -4.778742e+00 -4.452094e+00 -7.159972e+00
    ##  [171] -7.752110e+00 -7.247567e+00 -8.770264e+00 -9.129139e+00 -9.599266e+00
    ##  [176] -9.129943e+00 -1.060528e+01 -1.152124e+01 -1.106326e+01 -1.177577e+01
    ##  [181] -1.215420e+01 -1.014734e+01 -1.069751e+01 -1.267830e+01 -1.164370e+01
    ##  [186] -1.248854e+01 -1.350549e+01 -1.283073e+01 -1.286396e+01 -1.061847e+01
    ##  [191] -1.107261e+01 -1.087774e+01 -1.084624e+01 -1.263557e+01 -1.270606e+01
    ##  [196] -1.090192e+01 -1.045728e+01 -1.164822e+01 -1.140000e+01 -1.151897e+01
    ##  [201] -1.083199e+01 -1.223976e+01 -1.166444e+01 -1.102674e+01 -1.197486e+01
    ##  [206] -1.349039e+01 -1.322386e+01 -1.293232e+01 -1.122340e+01 -1.074670e+01
    ##  [211] -9.439148e+00 -7.633136e+00 -6.597750e+00 -6.081231e+00 -6.475335e+00
    ##  [216] -7.985575e+00 -7.979585e+00 -7.124578e+00 -6.394727e+00 -6.045968e+00
    ##  [221] -4.684601e+00 -4.973263e+00 -3.084859e+00 -3.102235e+00 -2.561152e+00
    ##  [226] -3.101193e+00 -3.182042e+00 -3.211676e+00 -3.356927e+00 -3.367789e+00
    ##  [231] -4.746508e+00 -4.074972e+00 -5.229178e+00 -4.806580e+00 -5.193570e+00
    ##  [236] -6.434096e+00 -5.977134e+00 -7.387927e+00 -5.454020e+00 -4.884459e+00
    ##  [241] -3.951679e+00 -5.497026e+00 -4.895408e+00 -5.846174e+00 -4.485808e+00
    ##  [246] -4.752332e+00 -7.787278e+00 -8.026098e+00 -8.298284e+00 -7.107887e+00
    ##  [251] -7.056691e+00 -6.725643e+00 -6.194024e+00 -5.579993e+00 -5.503277e+00
    ##  [256] -5.764156e+00 -6.686027e+00 -6.568420e+00 -8.706831e+00 -9.026250e+00
    ##  [261] -9.766945e+00 -1.069274e+01 -9.934802e+00 -9.575242e+00 -1.040448e+01
    ##  [266] -1.085023e+01 -1.152114e+01 -1.232206e+01 -1.071111e+01 -1.069259e+01
    ##  [271] -1.177599e+01 -1.048654e+01 -1.063046e+01 -1.291846e+01 -1.300318e+01
    ##  [276] -1.280412e+01 -1.157430e+01 -1.156974e+01 -1.286435e+01 -1.326151e+01
    ##  [281] -1.152203e+01 -1.049901e+01 -1.102658e+01 -1.018058e+01 -9.359031e+00
    ##  [286] -1.069480e+01 -9.651683e+00 -1.077467e+01 -1.008375e+01 -9.715145e+00
    ##  [291] -9.066631e+00 -9.506462e+00 -1.037288e+01 -9.507582e+00 -1.001740e+01
    ##  [296] -1.077617e+01 -7.374301e+00 -7.824717e+00 -6.095519e+00 -6.814007e+00
    ##  [301] -6.280388e+00 -6.350575e+00 -7.901138e+00 -7.627711e+00 -5.742476e+00
    ##  [306] -4.899277e+00 -4.564743e+00 -4.544357e+00 -5.551649e+00 -5.138955e+00
    ##  [311] -5.905987e+00 -4.980996e+00 -4.671846e+00 -2.705001e+00 -3.642538e+00
    ##  [316] -3.222065e+00 -3.761267e+00 -2.892673e+00 -1.742515e+00 -1.724819e+00
    ##  [321] -2.652664e+00 -3.158396e+00 -1.932662e+00 -7.222847e-01 -2.420558e+00
    ##  [326] -2.325862e+00 -1.563723e+00  8.034809e-01  7.597781e-01  8.268761e-01
    ##  [331]  9.821730e-01  4.493001e-01 -1.379673e+00 -7.791188e-01 -6.937787e-01
    ##  [336] -1.090352e+00 -2.707767e+00 -2.730672e+00 -2.324391e+00 -1.815345e+00
    ##  [341] -2.286739e-01  8.176379e-01 -2.303004e-01 -1.020901e+00 -5.536027e-01
    ##  [346]  2.241328e-04 -6.807543e-01 -3.475714e+00 -4.560300e+00 -5.217398e+00
    ##  [351] -6.186689e+00 -5.773925e+00 -6.689336e+00 -7.137652e+00 -7.936154e+00
    ##  [356] -8.953441e+00 -7.938453e+00 -8.589960e+00 -8.551744e+00 -8.367666e+00
    ##  [361] -7.915943e+00 -7.608063e+00 -7.323205e+00 -6.658654e+00 -6.989629e+00
    ##  [366] -5.498030e+00 -5.325940e+00 -5.134780e+00 -5.414395e+00 -4.639375e+00
    ##  [371] -4.975221e+00 -4.985750e+00 -5.124885e+00 -5.178351e+00 -3.848855e+00
    ##  [376] -2.859086e+00 -1.194520e+00 -3.314403e+00 -3.172923e+00 -3.350729e+00
    ##  [381] -3.052845e+00 -1.913412e+00 -1.786674e+00 -2.272624e+00 -3.276456e+00
    ##  [386] -1.492253e+00 -2.359397e+00 -3.329972e+00 -4.532847e+00 -4.900160e+00
    ##  [391] -4.886387e+00 -3.566298e+00 -3.077024e+00 -4.279405e+00 -5.152962e+00
    ##  [396] -5.041008e+00 -3.526776e+00 -3.150217e+00 -1.852935e+00  7.094214e-01
    ##  [401] -2.926033e-01  5.080641e-01 -5.279241e-01  8.379534e-01  6.087323e-01
    ##  [406]  3.409192e-01  7.146103e-01  1.695146e+00  3.976437e+00  4.508762e+00
    ##  [411]  5.007593e+00  5.161694e+00  5.244467e+00  3.658489e+00  2.128934e+00
    ##  [416]  3.853090e+00  4.057659e+00  5.295468e+00  5.538432e+00  6.551975e+00
    ##  [421]  7.834465e+00  8.123674e+00  5.719988e+00  7.988582e+00  7.377542e+00
    ##  [426]  6.914062e+00  5.432350e+00  5.128210e+00  5.441475e+00  5.176456e+00
    ##  [431]  3.816118e+00  4.210523e+00  4.749239e+00  6.886531e+00  7.093829e+00
    ##  [436]  5.817190e+00  5.301520e+00  5.222624e+00  5.860973e+00  6.673898e+00
    ##  [441]  5.780713e+00  5.379118e+00  6.203246e+00  4.899235e+00  5.772501e+00
    ##  [446]  6.462206e+00  5.608812e+00  6.262978e+00  7.662987e+00  8.665419e+00
    ##  [451]  5.623309e+00  5.381035e+00  7.232845e+00  5.556063e+00  4.928046e+00
    ##  [456]  5.181379e+00  5.156190e+00  4.969570e+00  4.736190e+00  7.460397e+00
    ##  [461]  4.667804e+00  6.319221e+00  8.928768e+00  8.325521e+00  8.735952e+00
    ##  [466]  9.287820e+00  7.165673e+00  7.051774e+00  6.778102e+00  6.986479e+00
    ##  [471]  9.351347e+00  1.095691e+01  1.261080e+01  1.482507e+01  1.438395e+01
    ##  [476]  1.494729e+01  1.612738e+01  1.519361e+01  1.517160e+01  1.493836e+01
    ##  [481]  1.371421e+01  1.427852e+01  1.499665e+01  1.422465e+01  1.465889e+01
    ##  [486]  1.495505e+01  1.398340e+01  1.334250e+01  1.426542e+01  1.317913e+01
    ##  [491]  1.363922e+01  1.342107e+01  1.474408e+01  1.502938e+01  1.413491e+01
    ##  [496]  1.411715e+01  1.453946e+01  1.466283e+01  1.300435e+01  1.287648e+01
    ##  [501]  1.198590e+01  1.388289e+01  1.623382e+01  1.620177e+01  1.771542e+01
    ##  [506]  1.905513e+01  1.979102e+01  2.123612e+01  2.360414e+01  2.378263e+01
    ##  [511]  2.497395e+01  2.742077e+01  2.741734e+01  2.903271e+01  3.027730e+01
    ##  [516]  3.057195e+01  3.001073e+01  3.025301e+01  3.279511e+01  3.373044e+01
    ##  [521]  3.376386e+01  3.513973e+01  3.411892e+01  3.536371e+01  3.638807e+01
    ##  [526]  3.645066e+01  3.746053e+01  3.724159e+01  3.893222e+01  4.039323e+01
    ##  [531]  4.195909e+01  4.291673e+01  4.441122e+01  4.463804e+01  4.721890e+01
    ##  [536]  4.639374e+01  4.772061e+01  4.843992e+01  5.033802e+01  5.069649e+01
    ##  [541]  5.183769e+01  5.118318e+01  5.019845e+01  5.068413e+01  4.972085e+01
    ##  [546]  5.102402e+01  5.263778e+01  5.438646e+01  5.536072e+01  5.707767e+01
    ##  [551]  5.915679e+01  6.003318e+01  6.176282e+01  6.297277e+01  6.496008e+01
    ##  [556]  6.444364e+01  6.542528e+01  6.583590e+01  6.704216e+01  6.849430e+01
    ##  [561]  7.016035e+01  7.063338e+01  7.165851e+01  7.321050e+01  7.462551e+01
    ##  [566]  7.597953e+01  7.621493e+01  7.670583e+01  7.590840e+01  7.660843e+01
    ##  [571]  7.715821e+01  7.908457e+01  7.868322e+01  8.094930e+01  8.171174e+01
    ##  [576]  8.179743e+01  8.424080e+01  8.522710e+01  8.569670e+01  8.750650e+01
    ##  [581]  8.817127e+01  8.773213e+01  9.059342e+01  9.164660e+01  9.269497e+01
    ##  [586]  9.239866e+01  9.232666e+01  9.330613e+01  9.359764e+01  9.562159e+01
    ##  [591]  9.770394e+01  9.925871e+01  1.018455e+02  1.031044e+02  1.053695e+02
    ##  [596]  1.075248e+02  1.101982e+02  1.118279e+02  1.116534e+02  1.134523e+02
    ##  [601]  1.156623e+02  1.161175e+02  1.183003e+02  1.187441e+02  1.204542e+02
    ##  [606]  1.235285e+02  1.244459e+02  1.242544e+02  1.226496e+02  1.240675e+02
    ##  [611]  1.251995e+02  1.266975e+02  1.279146e+02  1.289683e+02  1.299047e+02
    ##  [616]  1.301724e+02  1.326725e+02  1.333696e+02  1.352233e+02  1.358215e+02
    ##  [621]  1.361079e+02  1.349765e+02  1.365218e+02  1.377393e+02  1.388945e+02
    ##  [626]  1.411440e+02  1.426409e+02  1.433318e+02  1.437826e+02  1.451401e+02
    ##  [631]  1.472235e+02  1.490647e+02  1.498002e+02  1.499554e+02  1.505317e+02
    ##  [636]  1.517993e+02  1.531131e+02  1.550815e+02  1.555607e+02  1.565549e+02
    ##  [641]  1.585078e+02  1.602236e+02  1.607021e+02  1.613549e+02  1.611613e+02
    ##  [646]  1.626742e+02  1.637144e+02  1.654091e+02  1.672494e+02  1.675943e+02
    ##  [651]  1.693570e+02  1.715545e+02  1.703765e+02  1.712274e+02  1.724724e+02
    ##  [656]  1.732614e+02  1.745680e+02  1.757848e+02  1.761181e+02  1.761543e+02
    ##  [661]  1.758618e+02  1.775519e+02  1.807492e+02  1.816499e+02  1.820769e+02
    ##  [666]  1.822000e+02  1.843076e+02  1.843942e+02  1.866198e+02  1.864666e+02
    ##  [671]  1.866645e+02  1.863632e+02  1.868190e+02  1.886815e+02  1.904538e+02
    ##  [676]  1.904593e+02  1.917996e+02  1.923630e+02  1.935580e+02  1.944768e+02
    ##  [681]  1.968811e+02  1.985425e+02  2.000724e+02  2.027344e+02  2.033718e+02
    ##  [686]  2.029307e+02  2.055693e+02  2.046459e+02  2.059557e+02  2.077845e+02
    ##  [691]  2.073225e+02  2.082190e+02  2.103510e+02  2.107454e+02  2.112914e+02
    ##  [696]  2.121938e+02  2.128567e+02  2.141069e+02  2.158320e+02  2.151481e+02
    ##  [701]  2.153907e+02  2.164149e+02  2.154838e+02  2.168153e+02  2.174392e+02
    ##  [706]  2.191752e+02  2.193419e+02  2.198618e+02  2.184551e+02  2.198680e+02
    ##  [711]  2.196936e+02  2.208344e+02  2.208720e+02  2.213679e+02  2.228691e+02
    ##  [716]  2.248991e+02  2.251209e+02  2.267263e+02  2.265271e+02  2.259684e+02
    ##  [721]  2.272937e+02  2.288366e+02  2.301410e+02  2.304572e+02  2.315579e+02
    ##  [726]  2.334122e+02  2.340309e+02  2.344557e+02  2.345397e+02  2.348370e+02
    ##  [731]  2.361683e+02  2.386936e+02  2.389645e+02  2.422002e+02  2.436741e+02
    ##  [736]  2.462881e+02  2.475133e+02  2.477410e+02  2.489448e+02  2.510422e+02
    ##  [741]  2.517738e+02  2.513205e+02  2.525585e+02  2.542016e+02  2.544932e+02
    ##  [746]  2.548801e+02  2.560184e+02  2.576897e+02  2.592502e+02  2.615511e+02
    ##  [751]  2.606216e+02  2.593409e+02  2.591042e+02  2.603087e+02  2.622008e+02
    ##  [756]  2.637474e+02  2.626219e+02  2.638882e+02  2.647717e+02  2.652274e+02
    ##  [761]  2.659256e+02  2.668835e+02  2.685694e+02  2.684625e+02  2.716370e+02
    ##  [766]  2.723247e+02  2.737181e+02  2.743237e+02  2.771632e+02  2.761338e+02
    ##  [771]  2.755759e+02  2.780915e+02  2.780814e+02  2.801711e+02  2.797161e+02
    ##  [776]  2.819606e+02  2.825287e+02  2.835356e+02  2.846601e+02  2.852505e+02
    ##  [781]  2.868139e+02  2.894208e+02  2.894189e+02  2.909282e+02  2.922819e+02
    ##  [786]  2.934788e+02  2.955029e+02  2.955489e+02  2.970182e+02  2.979341e+02
    ##  [791]  2.978036e+02  2.985731e+02  2.995399e+02  2.995368e+02  3.006084e+02
    ##  [796]  3.031958e+02  3.046235e+02  3.059318e+02  3.070534e+02  3.067364e+02
    ##  [801]  3.068268e+02  3.072789e+02  3.074928e+02  3.088546e+02  3.095353e+02
    ##  [806]  3.109886e+02  3.121760e+02  3.139894e+02  3.138402e+02  3.138379e+02
    ##  [811]  3.161082e+02  3.173335e+02  3.169046e+02  3.172680e+02  3.182923e+02
    ##  [816]  3.209323e+02  3.244694e+02  3.255003e+02  3.266632e+02  3.269664e+02
    ##  [821]  3.298289e+02  3.312237e+02  3.313924e+02  3.308157e+02  3.323260e+02
    ##  [826]  3.318796e+02  3.317303e+02  3.312781e+02  3.318616e+02  3.345908e+02
    ##  [831]  3.350540e+02  3.371265e+02  3.376981e+02  3.397062e+02  3.396548e+02
    ##  [836]  3.410713e+02  3.428667e+02  3.453351e+02  3.464917e+02  3.474856e+02
    ##  [841]  3.490035e+02  3.491915e+02  3.507096e+02  3.511719e+02  3.528818e+02
    ##  [846]  3.560210e+02  3.577270e+02  3.611298e+02  3.616706e+02  3.623940e+02
    ##  [851]  3.624821e+02  3.632289e+02  3.628840e+02  3.631263e+02  3.640667e+02
    ##  [856]  3.651593e+02  3.657499e+02  3.668139e+02  3.670035e+02  3.679174e+02
    ##  [861]  3.691813e+02  3.688735e+02  3.709597e+02  3.729848e+02  3.733767e+02
    ##  [866]  3.744640e+02  3.769398e+02  3.802811e+02  3.792731e+02  3.796499e+02
    ##  [871]  3.811048e+02  3.817130e+02  3.827247e+02  3.835014e+02  3.833132e+02
    ##  [876]  3.842003e+02  3.846321e+02  3.829562e+02  3.837630e+02  3.854509e+02
    ##  [881]  3.880550e+02  3.902615e+02  3.920023e+02  3.932436e+02  3.952963e+02
    ##  [886]  3.980887e+02  4.004699e+02  4.021247e+02  4.050204e+02  4.062809e+02
    ##  [891]  4.061538e+02  4.060484e+02  4.064276e+02  4.058978e+02  4.079868e+02
    ##  [896]  4.088339e+02  4.100957e+02  4.119119e+02  4.106460e+02  4.106055e+02
    ##  [901]  4.106634e+02  4.127763e+02  4.143994e+02  4.174109e+02  4.165586e+02
    ##  [906]  4.177375e+02  4.184022e+02  4.195771e+02  4.202751e+02  4.216760e+02
    ##  [911]  4.224745e+02  4.228241e+02  4.237926e+02  4.260460e+02  4.279495e+02
    ##  [916]  4.282667e+02  4.277246e+02  4.292951e+02  4.291071e+02  4.297215e+02
    ##  [921]  4.293964e+02  4.300957e+02  4.308392e+02  4.309065e+02  4.330166e+02
    ##  [926]  4.332477e+02  4.347585e+02  4.377736e+02  4.408113e+02  4.428481e+02
    ##  [931]  4.439613e+02  4.453593e+02  4.464030e+02  4.466109e+02  4.470803e+02
    ##  [936]  4.495817e+02  4.511082e+02  4.523731e+02  4.552804e+02  4.569770e+02
    ##  [941]  4.576690e+02  4.582010e+02  4.593898e+02  4.610232e+02  4.610438e+02
    ##  [946]  4.595944e+02  4.606766e+02  4.627212e+02  4.637006e+02  4.646540e+02
    ##  [951]  4.654184e+02  4.674699e+02  4.684417e+02  4.700941e+02  4.723410e+02
    ##  [956]  4.742874e+02  4.753199e+02  4.770970e+02  4.787167e+02  4.792941e+02
    ##  [961]  4.807537e+02  4.817955e+02  4.831900e+02  4.843792e+02  4.858820e+02
    ##  [966]  4.867271e+02  4.877912e+02  4.884088e+02  4.896297e+02  4.912485e+02
    ##  [971]  4.921240e+02  4.918114e+02  4.933513e+02  4.938547e+02  4.953539e+02
    ##  [976]  4.978953e+02  4.992285e+02  5.008406e+02  5.032723e+02  5.046974e+02
    ##  [981]  5.053530e+02  5.048928e+02  5.055438e+02  5.071037e+02  5.099142e+02
    ##  [986]  5.113574e+02  5.126485e+02  5.137747e+02  5.161571e+02  5.179962e+02
    ##  [991]  5.186798e+02  5.211657e+02  5.223699e+02  5.218401e+02  5.245719e+02
    ##  [996]  5.265174e+02  5.264127e+02  5.267383e+02  5.272738e+02  5.278706e+02
    ## [1001]  5.296899e+02
    ## 
    ## $data_length
    ## [1] 1000
    ## 
    ## $current_penalty
    ## [1] 17.03439
    ## 
    ## $model
    ## [1] "gauss"
    ## 
    ## $pruning_algo
    ## [1] 3
    ## 
    ## $pruning_constraints_type
    ## [1] 0
    ## 
    ## $pruning_nb_loops
    ## [1] 5

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

## dust Algorithms

We start with one simple example with the `exp` model:

    data <- dataGenerator_1D(chpts = c(5*1e5,1e6), parameters = c(2,1), type = "exp")
    system.time(res <- dust.1D(data = data, model = "exp"))[[1]]

    ## [1] 0.799

    res$changepoints

    ## [1]  499991 1000000

A fundamental information relies in the number of indices to consider at
each data step. It is saved into the field `nb`.

[(Back to Top)](#top)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

<a id="cost"></a>

## Cost Integration

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

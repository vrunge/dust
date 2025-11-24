<a id="top"></a>

# dust Vignette

### Vincent Runge and Simon Querné

#### LaMME, Evry University, November 24, 2025

<center>
<img src="man/figures/dust.png" alt="" style="width:30%;"/>
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

> [Quick start](#start)

> [Rcpp Object Structure](#rcpp)

> [dust Algorithms](#dust1D)

> [Models And Data Generators](#Models)

> [New cost integration](#cost)

> [Pruning Capacity](#pruning)

<center>
<img src="man/figures/sep.png" alt="" style="width:100%;"/>
</center>

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
    ##  [1] 150 149 148 147 145 144 143 142 141 140 139 136 135 132 130 128 124 117 108
    ## [20] 105 102 101 100  99  98  97  96  95  52  51  50  49  48  47  46  44  43  42
    ## [39]  40  39  37  36   0
    ## 
    ## $nb
    ##   [1]  1  1  1  2  3  4  2  2  2  2  3  4  5  2  2  3  2  3  2  3  4  4  3  4  5
    ##  [26]  5  4  3  3  4  4  4  4  5  6  7  7  7  7  8  8  9  9  9 10 11 12 13 14 15
    ##  [51] 15 16 17 18 19 20 21 21 21 22 23 23 20 21 22 23 22 22 23 24 25 25 24 25 24
    ##  [76] 22 23 24 24 25 23 23 23 23 24 25 25 26 25 24 25 24 24 25 26 27 28 27 27 26
    ## [101] 25 26 25 26 27 28 29 30 31 32 33 34 35 36 36 37 37 38 39 40 40 40 41 42 42
    ## [126] 43 42 43 43 42 43 44 44 45 42 41 42 43 44 43 43 43 43 44 44 43 43 44 43 42
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

We can create a `DUST_1D` object with the `new` operator as follows:

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

The object is empty. We append data to this object using `append_c` and
run the dynamic programming algorithm with `update_partition()`. The
current partition is returned using a backtracking algorithm with
`get_partition()`.

    set.seed(1)
    data <- dataGenerator_1D(chpts = c(50,100), parameters = c(0,2), type = "gauss")
    obj_dust$append_c(data, 2*log(100))

    obj_dust$update_partition()
    obj_dust$get_partition()

    ## $changepoints
    ## [1]  50 100
    ## 
    ## $lastIndexSet
    ## [1] 100  98  97  96  50  49
    ## 
    ## $nb
    ##   [1] 1 1 1 2 1 1 1 1 1 1 2 2 2 3 3 3 2 2 3 2 3 3 3 3 4 4 3 4 4 5 3 4 3 2 3 3 3
    ##  [38] 3 2 2 3 2 2 2 2 3 3 3 1 2 3 3 4 4 5 6 6 6 4 4 4 5 4 4 4 3 4 3 4 5 5 5 4 4
    ##  [75] 5 5 6 4 4 5 6 6 4 4 5 4 4 3 2 2 3 3 4 5 4 4 4 5 5 5
    ## 
    ## $costQ
    ##   [1] -1.962222e-01 -4.902028e-02 -2.724011e-01 -1.254858e-02 -4.177676e-02
    ##   [6] -2.526450e-03 -7.011660e-03 -6.912099e-02 -1.471380e-01 -8.738788e-02
    ##  [11] -3.650215e-01 -4.329972e-01 -2.604826e-01 -5.368586e-03 -7.626959e-02
    ##  [16] -6.731781e-02 -6.196785e-02 -1.593813e-01 -2.722727e-01 -3.629935e-01
    ##  [21] -5.325653e-01 -6.904009e-01 -6.783726e-01 -2.695209e-01 -3.555994e-01
    ##  [26] -3.328802e-01 -2.969939e-01 -1.146593e-01 -7.286773e-02 -1.019902e-01
    ##  [31] -2.368948e-01 -2.173467e-01 -2.568520e-01 -2.428244e-01 -1.030997e-01
    ##  [36] -7.165948e-02 -4.761794e-02 -4.348111e-02 -1.091535e-01 -1.693764e-01
    ##  [41] -1.508041e-01 -1.267646e-01 -1.823557e-01 -2.318338e-01 -1.628204e-01
    ##  [46] -1.058452e-01 -1.292134e-01 -1.884746e-01 -1.750043e-01 -2.522464e-01
    ##  [51] -5.398443e-01 -7.460534e-01 -1.172772e+00 -1.337837e+00 -2.170942e+00
    ##  [56] -8.348940e+00 -9.428418e+00 -9.104201e+00 -1.231803e+01 -1.401394e+01
    ##  [61] -2.138320e+01 -2.323650e+01 -2.678952e+01 -2.879922e+01 -2.905885e+01
    ##  [66] -3.145242e+01 -2.948902e+01 -3.464769e+01 -3.696485e+01 -4.381931e+01
    ##  [71] -4.686827e+01 -4.720795e+01 -5.055695e+01 -5.042462e+01 -4.965436e+01
    ##  [76] -5.227209e+01 -5.330251e+01 -5.529438e+01 -5.744360e+01 -5.818112e+01
    ##  [81] -5.897889e+01 -6.069232e+01 -6.516463e+01 -6.396799e+01 -6.719953e+01
    ##  [86] -6.989124e+01 -7.412467e+01 -7.547530e+01 -7.825036e+01 -8.080892e+01
    ##  [91] -8.165785e+01 -8.620663e+01 -9.068094e+01 -9.418007e+01 -9.962469e+01
    ##  [96] -1.028336e+02 -1.020216e+02 -1.027666e+02 -1.021315e+02 -1.031187e+02

We can add data to the current time series and run again the algorithm:

    set.seed(13)
    data <- rnorm(40, mean = 1)
    obj_dust$append_c(data, 2*log(100))

    obj_dust$update_partition()
    obj_dust$get_partition()

    ## $changepoints
    ## [1]  50  96 140
    ## 
    ## $lastIndexSet
    ##  [1] 140 139 138 126 123 110 108 107  96  50  49
    ## 
    ## $nb
    ##   [1] 1 1 1 2 1 1 1 1 1 1 2 2 2 3 3 3 2 2 3 2 3 3 3 3 4 4 3 4 4 5 3 4 3 2 3 3 3
    ##  [38] 3 2 2 3 2 2 2 2 3 3 3 1 2 3 3 4 4 5 6 6 6 4 4 4 5 4 4 4 3 4 3 4 5 5 5 4 4
    ##  [75] 5 5 6 4 4 5 6 6 4 4 5 4 4 3 2 2 3 3 4 5 4 4 4 5 5 5 4 5 5 4 4 4 3 4 4 5 6
    ## [112] 6 7 8 9 8 9 7 7 6 6 6 7 7 6 6 7 8 9 8 8 9 9 7 7 8 7 7 8 9
    ## 
    ## $costQ
    ##   [1] -1.962222e-01 -4.902028e-02 -2.724011e-01 -1.254858e-02 -4.177676e-02
    ##   [6] -2.526450e-03 -7.011660e-03 -6.912099e-02 -1.471380e-01 -8.738788e-02
    ##  [11] -3.650215e-01 -4.329972e-01 -2.604826e-01 -5.368586e-03 -7.626959e-02
    ##  [16] -6.731781e-02 -6.196785e-02 -1.593813e-01 -2.722727e-01 -3.629935e-01
    ##  [21] -5.325653e-01 -6.904009e-01 -6.783726e-01 -2.695209e-01 -3.555994e-01
    ##  [26] -3.328802e-01 -2.969939e-01 -1.146593e-01 -7.286773e-02 -1.019902e-01
    ##  [31] -2.368948e-01 -2.173467e-01 -2.568520e-01 -2.428244e-01 -1.030997e-01
    ##  [36] -7.165948e-02 -4.761794e-02 -4.348111e-02 -1.091535e-01 -1.693764e-01
    ##  [41] -1.508041e-01 -1.267646e-01 -1.823557e-01 -2.318338e-01 -1.628204e-01
    ##  [46] -1.058452e-01 -1.292134e-01 -1.884746e-01 -1.750043e-01 -2.522464e-01
    ##  [51] -5.398443e-01 -7.460534e-01 -1.172772e+00 -1.337837e+00 -2.170942e+00
    ##  [56] -8.348940e+00 -9.428418e+00 -9.104201e+00 -1.231803e+01 -1.401394e+01
    ##  [61] -2.138320e+01 -2.323650e+01 -2.678952e+01 -2.879922e+01 -2.905885e+01
    ##  [66] -3.145242e+01 -2.948902e+01 -3.464769e+01 -3.696485e+01 -4.381931e+01
    ##  [71] -4.686827e+01 -4.720795e+01 -5.055695e+01 -5.042462e+01 -4.965436e+01
    ##  [76] -5.227209e+01 -5.330251e+01 -5.529438e+01 -5.744360e+01 -5.818112e+01
    ##  [81] -5.897889e+01 -6.069232e+01 -6.516463e+01 -6.396799e+01 -6.719953e+01
    ##  [86] -6.989124e+01 -7.412467e+01 -7.547530e+01 -7.825036e+01 -8.080892e+01
    ##  [91] -8.165785e+01 -8.620663e+01 -9.068094e+01 -9.418007e+01 -9.962469e+01
    ##  [96] -1.028336e+02 -1.020216e+02 -1.027666e+02 -1.021315e+02 -1.031187e+02
    ## [101] -1.041713e+02 -1.034875e+02 -1.071009e+02 -1.074035e+02 -1.096965e+02
    ## [106] -1.104834e+02 -1.129555e+02 -1.133804e+02 -1.125916e+02 -1.148047e+02
    ## [111] -1.125889e+02 -1.135169e+02 -1.108678e+02 -1.073577e+02 -1.066222e+02
    ## [116] -1.063729e+02 -1.091076e+02 -1.094161e+02 -1.093347e+02 -1.107724e+02
    ## [121] -1.113996e+02 -1.149335e+02 -1.157326e+02 -1.139648e+02 -1.152571e+02
    ## [126] -1.156956e+02 -1.142671e+02 -1.124902e+02 -1.118586e+02 -1.116007e+02
    ## [131] -1.120919e+02 -1.134718e+02 -1.135704e+02 -1.135272e+02 -1.137484e+02
    ## [136] -1.136377e+02 -1.138061e+02 -1.140670e+02 -1.137241e+02 -1.134121e+02

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

    ## [1] 0.85

    res$changepoints

    ## [1]  500014 1000000

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

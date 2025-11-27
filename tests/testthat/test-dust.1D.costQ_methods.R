
library(testthat)
library(dust)


##############################################################
############## all methods  - all costs as OP ################
##############################################################

test_that("costQ all the same, for all types and methods",
          {
            types <- c( "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance")
            for(i in types)
            {
            data <- dataGenerator_1D(parameters = 0.5, type = i)
            data <- data_normalization_1D(data, type = i)
            res0 <- dust.1D(data, model = i, method = "rand_DUSTr")
            res1 <- dust.1D(data, model = i, method = "rand_DUST")
            res2 <- dust.1D(data, model = i, method = "rand_DUSTgs")
            res3 <- dust.1D(data, model = i, method = "rand_DUSTbs")
            res4 <- dust.1D(data, model = i, method = "rand_DUSTqn")
            res5 <- dust.1D(data, model = i, method = "rand_PELT")
            expect_equal(res0$costQ, res1$costQ)
            expect_equal(res0$costQ, res2$costQ)
            expect_equal(res0$costQ, res3$costQ)
            expect_equal(res0$costQ, res4$costQ)
            expect_equal(res0$costQ, res5$costQ)
            res0 <- dust.1D(data, model = i, method = "det_DUSTr")
            res1 <- dust.1D(data, model = i, method = "det_DUST")
            res2 <- dust.1D(data, model = i, method = "det_DUSTgs")
            res3 <- dust.1D(data, model = i, method = "det_DUSTbs")
            res4 <- dust.1D(data, model = i, method = "det_DUSTqn")
            res5 <- dust.1D(data, model = i, method = "det_PELT")
            expect_equal(res0$costQ, res1$costQ)
            expect_equal(res0$costQ, res2$costQ)
            expect_equal(res0$costQ, res3$costQ)
            expect_equal(res0$costQ, res4$costQ)
            expect_equal(res0$costQ, res5$costQ)
            }
          })


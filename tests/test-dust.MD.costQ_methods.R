
library(testthat)
library(dust)


##############################################################
############## all methods  - all costs as OP ################
##############################################################

test_that("costQ all the same, for all types and methods",
          {
            types <- c( "gauss", "exp", "variance")
            for(i in types)
            {
              d <- 5
              data <- dataGenerator_MD(chpts = c(40,60,100),
                                       parameters =data.frame(matrix(abs(rnorm(3*d)), ncol = d)),
                                       type = i)
              data <- data_normalization_MD(data, type = i)
              res0 <- dust.MD(data, model = i, method = "rand_DUSTr")
              res1 <- dust.MD(data, model = i, method = "rand_DUST")
              res2 <- dust.MD(data, model = i, method = "rand_DUSTgs")
              res3 <- dust.MD(data, model = i, method = "rand_DUSTbs")
              res4 <- dust.MD(data, model = i, method = "rand_DUSTqn")
              res5 <- dust.MD(data, model = i, method = "rand_PELT")
              expect_equal(res0$costQ, res1$costQ)
              expect_equal(res0$costQ, res2$costQ)
              expect_equal(res0$costQ, res3$costQ)
              expect_equal(res0$costQ, res4$costQ)
              expect_equal(res0$costQ, res5$costQ)
              res0 <- dust.MD(data, model = i, method = "det_DUSTr")
              res1 <- dust.MD(data, model = i, method = "det_DUST")
              res2 <- dust.MD(data, model = i, method = "det_DUSTgs")
              res3 <- dust.MD(data, model = i, method = "det_DUSTbs")
              res4 <- dust.MD(data, model = i, method = "det_DUSTqn")
              res5 <- dust.MD(data, model = i, method = "det_PELT")
              expect_equal(res0$costQ, res1$costQ)
              expect_equal(res0$costQ, res2$costQ)
              expect_equal(res0$costQ, res3$costQ)
              expect_equal(res0$costQ, res4$costQ)
              expect_equal(res0$costQ, res5$costQ)
            }
          })





test_that("costQ all the same, for all types and methods",
          {
            types <- c( "bern", "binom")
            for(i in types)
            {
              d <- 5
              param <- data.frame(matrix(abs(rnorm(3*d)), ncol = d))
              param <- param/max(param)
              data <- dataGenerator_MD(chpts = c(40,60,100),
                                       parameters =param,
                                       type = i)
              data <- data_normalization_MD(data, type = i)
              res0 <- dust.MD(data, model = i, method = "rand_DUSTr")
              res1 <- dust.MD(data, model = i, method = "rand_DUST")
              res2 <- dust.MD(data, model = i, method = "rand_DUSTgs")
              res3 <- dust.MD(data, model = i, method = "rand_DUSTbs")
              res4 <- dust.MD(data, model = i, method = "rand_DUSTqn")
              res5 <- dust.MD(data, model = i, method = "rand_PELT")
              expect_equal(res0$costQ, res1$costQ)
              expect_equal(res0$costQ, res2$costQ)
              expect_equal(res0$costQ, res3$costQ)
              expect_equal(res0$costQ, res4$costQ)
              expect_equal(res0$costQ, res5$costQ)
              res0 <- dust.MD(data, model = i, method = "det_DUSTr")
              res1 <- dust.MD(data, model = i, method = "det_DUST")
              res2 <- dust.MD(data, model = i, method = "det_DUSTgs")
              res3 <- dust.MD(data, model = i, method = "det_DUSTbs")
              res4 <- dust.MD(data, model = i, method = "det_DUSTqn")
              res5 <- dust.MD(data, model = i, method = "det_PELT")
              expect_equal(res0$costQ, res1$costQ)
              expect_equal(res0$costQ, res2$costQ)
              expect_equal(res0$costQ, res3$costQ)
              expect_equal(res0$costQ, res4$costQ)
              expect_equal(res0$costQ, res5$costQ)
            }
          })


context("plot3D(): test plot3D function")

test_that("levelplot in plot3D works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("copula")
  skip_if_not_installed("plotly")
  skip('Skip levelplot in plot3D')

  library(copula)
  library(plotly)
  # Load data
  data("ANES2016")

  PAsso_4v <- PAsso(responses = c("PID", "selfLR", "TrumpLR", "ClinLR"),
                    adjustments = c("age", "edu.year", "income.num"),
                    data = ANES2016,
                    method = "kendall")

  testPlot <- plot3D(PAsso_4v, y1="PID", y2="selfLR", type = "contour")
  testPlots <- plot3D(PAsso_4v, type = "contour")

  # Expectations
  expect_s3_class(testPlot, "plotly")
  expect_is(testPlots, "list")

})

test_that("plot3D works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("copula")
  skip_if_not_installed("plotly")

  skip("Skip plot3D, it takes time and no display")

  library(copula)
  library(plotly)
  # Load data
  data("ANES2016")

  # multivariate analysis (2 variables) --------------------------------------------------------------------
  PAsso_2v <- PAsso(responses = c("PreVote.num", "PID"),
                    adjustments = c("income.num", "age", "edu.year"),
                    data = ANES2016, uni.model = "logit",
                    method = c("kendall"),
                    resids.type = "surrogate", jitter = "latent")

  testPlots <- plot3D(PAsso_2v, y1="PreVote.num", y2="PID")

  # Expectations
  expect_s3_class(testPlots, "plotly")

})


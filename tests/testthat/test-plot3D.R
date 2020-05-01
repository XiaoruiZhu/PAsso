context("PAsso: plot3D")

test_that("plot3D works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("copula")
  skip_if_not_installed("plotly")

  library(copula)
  library(plotly)
  # Load data
  data("nes2016")

  # multivariate analysis (5 variables) --------------------------------------------------------------------
  PAsso_2v <- PAsso(responses = c("Prevote.num", "PID"),
                    adjustments = c("income.num", "age", "edu.year"),
                    data = nes2016, uni.model = "logit",
                    method = c("kendall"),
                    resids.type = "surrogate", jitter = "latent")

  testPlots <- plot3D(PAsso_2v)

  # Expectations
  expect_s3_class(testPlots$plot_1, "plotly")

})


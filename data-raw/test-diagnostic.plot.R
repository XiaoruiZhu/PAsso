context("PAsso: 'diagnostic.plot()' Residual plots work for objects")


test_that("diagnostic.plot works for \"clm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "logit")

  # Construct residual plots
  p1 <- diagnostic.plot(fit, type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})


test_that("diagnostic.plot works for \"glm\" objects", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- glm(y ~ x + I(x ^ 2), data = df1, family = binomial)

  # Construct residual plots
  p1 <- diagnostic.plot(fit, jitter.scale = "probability", type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})


test_that("diagnostic.plot works for \"lrm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::lrm(y ~ x, data = df1)

  # Construct residual plots
  p1 <- diagnostic.plot(fit, type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})


test_that("diagnostic.plot works for \"orm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::orm(y ~ x, data = df1, family = logistic)

  # Construct residual plots
  p1 <- diagnostic.plot(fit, type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})


test_that("diagnostic.plot works for \"polr\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "logistic")

  # Construct residual plots
  p1 <- diagnostic.plot(fit, type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})


test_that("diagnostic.plot works for \"vglm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  suppressWarnings(
    fit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                      family = VGAM::cumulative(link = "logit",
                                                parallel = TRUE))
  )

  # Construct residual plots
  p1 <- diagnostic.plot(fit, type = "qq")
  p2 <- diagnostic.plot(fit, type = "covariate", x = df1$x)
  p3 <- diagnostic.plot(fit, type = "fitted")
  # p4 <- diagnostic.plot(fit, type = "fitted", nsim = 10)


  # Expectations
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
  expect_is(p3, "ggplot")
  # expect_is(p4, "ggplot")

})

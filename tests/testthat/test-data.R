context("PAsso: Data section")

test_that("df_AdjCat", {
  #
  # Adjacent Categories Regression Model Example to compare different residuals
  #

  data("df_AdjCat")

  # Expectations
  expect_is(df_AdjCat, "list")
})

test_that("nes2016", {
  #
  # Adjacent Categories Regression Model Example to compare different residuals
  #

  data("nes2016")

  # Expectations
  expect_is(nes2016, "data.frame")
})

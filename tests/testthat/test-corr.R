context("corr: marginal and partial association")


test_that("marginal correlation for two responses", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("parcor")
  skip_if_not_installed("tidyverse")
  library("MASS")
  library("parcor")
  library("tidyverse")

  # Load data
  data("nes96")

  # marginal association (Kendall's tau)
  tau <- cor(nes96[,c("vote.num", "PID")],method = "kendall")
  # standard error from bootstrap
  tau.sd.boot <- sd(sapply(1:200, function(b){
    index<- sample(nrow(nes96), replace=T)
    cor(nes96[index, c("vote.num", "PID")], method = "kendall")
  }))
  # tau; tau.sd.boot

  MAsso_1 <- corr(responses = c("vote.num", "PID"),
                  adjustments = c("income.num", "age", "edu.year"),
                  data = nes96,
                  association = c("marginal")
                  # models = c("probit", "probit"),
                  # method = c("kendall"),
                  # resids.method = "latent", fitted.models = NULL,
                  # rep_num = 100
  )


  # Expectations
  expect_equal(MAsso_1[1,2], tau[1,2])
  expect_is(MAsso_1, "MAsso")
})


test_that("partial correlation for two responses",{

  # "corr" advanced using of the function: The First way (Advanced), input a few models directly ------------------------------

  y1 <- nes96$vote.num
  y2 <- nes96$PID
  X <- as.matrix(nes96[c("income.num", "age", "edu.year")])
  fit.vote<- glm(y1 ~ X, family = binomial(link = "probit"))
  fit.PID<- polr(as.factor(y2)~ X, method="probit")
  fitted.temp <- list(fit.vote, fit.PID)
  PAsso_adv1 <- corr(fitted.models=fitted.temp,
                     association = c("partial"),
                     method = c("kendall"),
                     resids.method = "latent", rep_num=100)

  # Partial association coefficients (Parts of Table 7 in paper)
  expect_equal(dim(PAsso_adv1$corr)[1], length(fitted.temp))
  expect_is(PAsso_adv1, "PAsso")

  # "corr" function: The simple way, input response and confounders only ----------------------------
  PAsso_1 <- corr(responses = c("vote.num", "PID"),
                  adjustments = c("income.num", "age", "edu.year"),
                  data = nes96
                  # association = c("partial"),
                  # models = c("probit", "probit"),
                  # method = c("kendall"),
                  # resids.method = "latent", fitted.models = NULL,
                  # rep_num = 100
  )

  # Compare marginal correlation and partial correlation.
  expect_equal(dim(PAsso_adv1$corr)[1], dim(PAsso_1$corr)[1])

})


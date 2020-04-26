
library(MASS)
library(parcor)
library(parasol)
library(tidyverse)
library(matrixStats)
library(progress)

# import data -------------------------------------------------------------
data(nes96)
# parasol::nes96
summary(nes96)

# Bivariate analysis of VOTE and PID (Kendall's tao) ----------------------

# marginal association (Kendall's tau)
tau <- cor(nes96[,c("vote.num", "PID")],method = "kendall")
# standard error from bootstrap
tau.sd.boot <- sd(sapply(1:200, function(b){
  index<- sample(nrow(nes96), replace=T)
  cor(nes96[index, c("vote.num", "PID")], method = "kendall")
}))
tau; tau.sd.boot

# "PAsso" advanced using of the function: The First way (Advanced), input a few models directly ------------------------------
# Test "PAsso" function: Partial Association by surrogate residuals regression models
fit.vote<- glm(vote.num ~ income.num+ age + edu.year, data = nes96,
               family = binomial(link = "probit"))
fit.PID<- polr(as.factor(PID) ~ income.num+age+edu.year, data = nes96,
               method="probit", Hess = TRUE)

PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
                   association = c("partial"),
                   method = c("kendall"),
                   resids.method = "latent", rep_num=100)

# Partial association coefficients (Parts of Table 7 in paper)
PAsso_adv1$corr
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)

# "PAsso" function: The simple way, input response and confounders only ----------------------------
PAsso_1 <- PAsso(responses = c("vote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = nes96
                   # association = c("partial"),
                   # models = c("probit", "probit"),
                   # method = c("kendall"),
                   # resids.method = "latent", fitted.models = NULL,
                   # rep_num = 100
                )

# Compare marginal correlation and partial correlation.
tau; tau.sd.boot
PAsso_adv1$corr;
PAsso_1$corr;
print(PAsso_1, 5)
summary(PAsso_1)
plot(PAsso_1)

MAsso_1 <- PAsso(responses = c("vote.num", "PID"),
                data = nes96,
                association = c("marginal"))
summary(MAsso_1, digits=3)

# "PAsso" function: input three responses ----------------------------
PAsso_2 <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes96,
                models = c("probit", "probit", "probit"),
                association = c("partial"), method = c("kendall"),
                resids.method = "latent",
                rep_num = 100)

# Compare marginal correlation and partial correlation.
tau; tau.sd.boot
PAsso_adv1$corr;
PAsso_1$corr;
PAsso_2
print(PAsso_2, digits=3)
summary(PAsso_2, digits=3)
plot(PAsso_2)

# "PAsso" function: input three responses ----------------------------
PAsso_2_marg <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                data = nes96, association = c("marginal"))

# Compare marginal correlation and partial correlation.
PAsso_2
PAsso_2_marg
summary(PAsso_2_marg, digits=5)

# Pcor_SR.test function: Conduct inference based on object of "PartialCor" class ----------------------------
library(progress)

system.time(Pcor_SR_test1 <- test(object = PAsso_1, boot_SE=20, H0=0, parallel=FALSE))
Pcor_SR_test1$corr; Pcor_SR_test1$corr_stat; Pcor_SR_test1$corr_p.value; Pcor_SR_test1$CI_95

print(Pcor_SR_test1)
summary(PAsso_1, digits=4)
summary(Pcor_SR_test1, digits=4)

# Pcor_SR.test function: Test by parallel ----------------------------

library(doParallel); library(progress)

numCores <- detectCores()  # Not too aggressive!
cl <- makeCluster(numCores)
# registerDoSNOW(cl) # on Mac or Linux
registerDoParallel(cl) # Win

system.time(Pcor_SR_test2 <- test(object = PAsso_2, boot_SE=300, H0=0, parallel=TRUE))
Pcor_SR_test2$sd_MatCor; Pcor_SR_test2$corr; Pcor_SR_test2$corr_stat;
Pcor_SR_test2$corr_p.value;
Pcor_SR_test2$CI_95

print(Pcor_SR_test2)
summary(Pcor_SR_test2)

# test1 <- summary(Pcor_SR_test2$fitted.models[[1]])
# test2 <- Pcor_SR_test2$fitted.models[[2]]
# test2$zeta
#
# test3 <- summary(Pcor_SR_test2$fitted.models[[3]])

stopCluster(cl)

# multivariate analysis (5 variables) --------------------------------------------------------------------
# marginal correlation
dta <- dplyr::select(nes96, vote.num, PID, selfLR, ClinLR, DoleLR,
                     income.num, age, edu.year) %>%
  dplyr::mutate(vote.num=as.factor(vote.num), PID=as.factor(PID),
                selfLR=as.factor(selfLR), ClinLR=as.factor(ClinLR), DoleLR=as.factor(DoleLR))
head(dta)

# Marginal association matrix (Left part of Table 7 in paper)
Mcor_tau <- cor(
  dplyr::select(nes96, vote.num,
                PID, selfLR, ClinLR, DoleLR),
  method = "kendall")
Mcor_tau[lower.tri(Mcor_tau)] <- NA
print(Mcor_tau, digits=3)

# regression model
fit.selfLR<- polr(selfLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.ClinLR<- polr(ClinLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.DoleLR<- polr(DoleLR ~ income.num + age + edu.year,
                  data = dta, method="probit")

# surrogate residual
PAsso_adv_5v <- PAsso(
  fitted.models = list(fit.vote, fit.PID, fit.selfLR, fit.ClinLR, fit.DoleLR),
  association = c("partial"),
  method = c("kendall"), resids.method = "latent", rep_num=30)

PAsso_adv_5v # (Right part of Table 7 in paper)
print(PAsso_adv_5v, digits=4)
summary(PAsso_adv_5v, digits=4)

PAsso_5v <- PAsso(responses = c("vote.num", "PID", "selfLR", "ClinLR", "DoleLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes96,
                 association = c("partial"),
                 method = c("kendall"), resids.method = "latent", rep_num=30)

PAsso_5v # (Right part of Table 7 in paper)
summary(PAsso_5v)

system.time(PAsso_5v_test1 <- test(object = PAsso_1, boot_SE=20, H0=0, parallel=FALSE))
PAsso_5v_test1$sd_MatCor; PAsso_5v_test1$corr; PAsso_5v_test1$corr_p.value
print(PAsso_5v_test1, 5)

system.time(PAsso_5v_test <- test(object = PAsso_5v, boot_SE=20, H0=0, parallel=FALSE))
PAsso_5v_test$sd_MatCor; PAsso_5v_test$corr;
PAsso_5v_test$corr_stat; PAsso_5v_test$corr_p.value; PAsso_5v_test$CI_95
print(PAsso_5v_test, 3)

# Partial Regression plot matrix ------------------------------------------
ggpairs(object = PAsso_1, colour="blue")
plot(object = PAsso_1, colour="blue")
plot(object = PAsso_adv1, colour="blue")
plot(object = PAsso_2, colour="blue")
plot(object = PAsso_5v, colour="blue")
# plot(object = PAsso_5v_test, colour="blue")


# check.model function -----------------------------------------------------
library(parasol); library(ggplot2)
data("nes96")
PAsso_2 <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes96,
                models = c("probit", "probit", "probit"),
                association = c("partial"), method = c("kendall"),
                resids.method = "latent",
                rep_num = 100)

PAsso_2
print(PAsso_2, digits=3)

check_qq <- check.model(object = PAsso_2, color="blue",what = "qq")
do.call("grid.arrange", c(check_qq, ncol=1))

check_fitted <- check.model(object = PAsso_2, what = "fitted")
do.call("grid.arrange", c(check_fitted, ncol=1))

check_covar <- check.model(object = PAsso_2, what = "covariate")
do.call("grid.arrange", c(check_covar, ncol=3))

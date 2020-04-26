
library(MASS)
library(parcor)
library(parasol)
library(tidyverse)
library(matrixStats)
library(progress)

# import data -------------------------------------------------------------
data("nes2016_pre")
nes2016 <- nes2016_pre
summary(nes2016)

nes2016 <- nes2016 %>% rename(vote.num=Prevote.num)

# Bivariate analysis of VOTE and PID (Kendall's tao) ----------------------

# marginal association (Kendall's tau)
tau <- cor(nes2016[,c("vote.num", "PID")],method = "kendall")
# standard error from bootstrap
tau.sd.boot <- sd(sapply(1:200, function(b){
  index<- sample(nrow(nes2016), replace=T)
  cor(nes2016[index, c("vote.num", "PID")], method = "kendall")
}))
tau; tau.sd.boot

# "PAsso" advanced using of the function: The First way (Advanced), input a few models directly ------------------------------
# Test "PAsso" function: Partial Association by surrogate residuals regression models
fit.vote<- glm(vote.num ~ income.num+ age + edu.year, data = nes2016,
               family = binomial(link = "probit"))
fit.PID<- polr(as.factor(PID) ~ income.num+age+edu.year, data = nes2016,
               method="probit", Hess = TRUE)

system.time(PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
                   association = c("partial"),
                   method = c("kendall"),
                   resids.method = "latent")
)

system.time(PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                association = c("marginal"),
                                method = c("kendall"),
                                resids.method = "latent")
)
# Partial association coefficients (Parts of Table 7 in paper)
PAsso_adv1$corr
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)

# "PAsso" function: The simple way, input response and confounders only ----------------------------
PAsso_1 <- PAsso(responses = c("vote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = nes2016
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
                data = nes2016,
                association = c("marginal"))
summary(MAsso_1, digits=3)

# "PAsso" function: input three responses ----------------------------
PAsso_2 <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016,
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
                data = nes2016, association = c("marginal"))

# Compare marginal correlation and partial correlation.
PAsso_2
PAsso_2_marg
summary(PAsso_2_marg, digits=5)

# Pcor_SR.test function: Conduct inference based on object of "PartialCor" class ----------------------------
library(progress)

system.time(Pcor_SR_test1 <- test(object = PAsso_1, boot_SE=20, H0=0, parallel=FALSE))
Pcor_SR_test1$corr; Pcor_SR_test1$corr_stat; Pcor_SR_test1$corr_p.value; Pcor_SR_test1$CI_95

print(Pcor_SR_test1, digits=3)
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
dta <- dplyr::select(nes2016, vote.num, PID, selfLR, ClinLR, TrumpLR,
                     income.num, age, edu.year) %>%
  dplyr::mutate(vote.num=as.factor(vote.num), PID=as.factor(PID),
                selfLR=as.factor(selfLR), ClinLR=as.factor(ClinLR), TrumpLR=as.factor(TrumpLR))
head(dta)

# Marginal association matrix (Left part of Table 7 in paper)
Mcor_tau <- cor(
  dplyr::select(nes2016, vote.num,
                PID, selfLR, ClinLR, TrumpLR),
  method = "kendall")
Mcor_tau[lower.tri(Mcor_tau)] <- NA
print(Mcor_tau, digits=3)

# regression model
fit.selfLR<- polr(selfLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.ClinLR<- polr(ClinLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.TrumpLR<- polr(TrumpLR ~ income.num + age + edu.year,
                  data = dta, method="probit")

# surrogate residual
PAsso_adv_5v <- PAsso(
  fitted.models = list(fit.vote, fit.PID, fit.selfLR, fit.ClinLR, fit.TrumpLR),
  association = c("partial"),
  method = c("kendall"), resids.method = "latent", rep_num=30)

PAsso_adv_5v # (Right part of Table 7 in paper)
print(PAsso_adv_5v, digits=4)
summary(PAsso_adv_5v, digits=4)

PAsso_5v <- PAsso(responses = c("vote.num", "PID", "selfLR", "ClinLR", "TrumpLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
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
data("nes2016_pre")
PAsso_2 <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016_pre,
                models = c("probit", "probit", "probit"),
                association = c("partial"), method = c("kendall"),
                resids.method = "latent")

PAsso_2
print(PAsso_2, digits=3)
summary(PAsso_2, digits = 5)

check_qq <- check.model(object = PAsso_2, color="blue", output = "qq")

check_fitted <- check.model(object = PAsso_2, output = "fitted")

check_covar <- check.model(object = PAsso_2, output = "covariate")


## general association measure and 3-D plot for VOTE and PID
library("copula")
library("plotly")
# library("copBasic")

dim(PAsso_2$rep_SRs)
testPlots <- plot3D(PAsso_2)
testPlots$plot_1

WolfCop <- wolfCOP(para = data.frame(r_vote, r_PID), as.sample = TRUE)
WolfCop
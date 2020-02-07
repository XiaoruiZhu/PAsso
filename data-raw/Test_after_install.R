
library(MASS)
library(parcor)
library(parasol)
library(tidyverse)
library(matrixStats)

# import data -------------------------------------------------------------
data(nes96)
# parasol::nes96
summary(nes96)

# Bivariate analysis of VOTE and PID (Kendall's tao) ----------------------
y1 <- nes96$vote.num
y2 <- nes96$PID
X <- as.matrix(nes96[c("income.num", "age", "edu.year")])

# marginal association (Kendall's tau)
tau <- cor(nes96[,c("vote.num", "PID")],method = "kendall")
# standard error from bootstrap
tau.sd.boot <- sd(sapply(1:200, function(b){
  index<- sample(nrow(nes96), replace=T)
  cor(nes96[index, c("vote.num", "PID")], method = "kendall")
}))
tau; tau.sd.boot

# "corr" advanced using of the function: The First way (Advanced), input a few models directly ------------------------------
# Test "corr" function: Partial Association by surrogate residuals regression models
fit.vote<- glm(y1 ~ X, family = binomial(link = "probit"))
fit.PID<- polr(as.factor(y2)~ X, method="probit")

PAsso_adv1 <- corr(fitted.models=list(fit.vote, fit.PID),
                   association = c("partial"),
                   method = c("kendall"),
                   resids.method = "latent", rep_num=100)

# Partial association coefficients (Parts of Table 7 in paper)
PAsso_adv1$corr
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1)

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
tau; tau.sd.boot
PAsso_adv1$corr;
PAsso_1$corr;
PAsso_1

# "corr" function: input three responses ----------------------------
PAsso_2 <- corr(responses = c("vote.num", "PID", "selfLR"),
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

# Pcor_SR.test function: Conduct inference based on object of "PartialCor" class ----------------------------

system.time(Pcor_SR_test1 <- corr.test(object = PAsso_1, boot_SE=20, H0=0, parallel=FALSE))
Pcor_SR_test1$corr; Pcor_SR_test1$corr_stat; Pcor_SR_test1$corr_p.value; Pcor_SR_test1$CI_95

# Pcor_SR.test function: Test by parallel ----------------------------

library(parallel); library(doSNOW); library(progress)

system.time(Pcor_SR_test2 <- corr.test(object = PAsso_1, boot_SE=20, H0=0, parallel=TRUE))
Pcor_SR_test2$corr; Pcor_SR_test2$corr_stat; Pcor_SR_test2$corr_p.value; Pcor_SR_test2$CI_95

# multivariate analysis (5 variables) --------------------------------------------------------------------
y3<- nes96clean$selfLR
y4<- nes96clean$ClinLR
y5<- nes96clean$DoleLR
# marginal correlation
dta <- dplyr::select(nes96clean, vote.num, PID, selfLR, ClinLR, DoleLR,
                     income.num, age, edu.year) %>%
  dplyr::mutate(vote.num=as.factor(vote.num), PID=as.factor(PID),
                selfLR=as.factor(selfLR), ClinLR=as.factor(ClinLR), DoleLR=as.factor(DoleLR))
head(dta)

# Marginal association matrix (Left part of Table 7 in paper)
Mcor_tau <- cor(
  dplyr::select(nes96clean, vote.num,
                PID, selfLR, ClinLR, DoleLR),
  method = "kendall")
Mcor_tau

# regression model
fit.selfLR<- polr(selfLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.ClinLR<- polr(ClinLR ~ income.num + age + edu.year,
                  data = dta, method="probit")
fit.DoleLR<- polr(DoleLR ~ income.num + age + edu.year,
                  data = dta, method="probit")

# surrogate residual
PAsso_5v <- corr(
  fitted.models = list(fit.vote, fit.PID, fit.selfLR, fit.ClinLR, fit.DoleLR),
  association = c("partial"),
  method = c("kendall"), resids.method = "latent", rep_num=30)

PAsso_5v$corr # (Right part of Table 7 in paper)
summary(PAsso_5v)

# Partial Regression plot matrix ------------------------------------------
dim(PAsso_5v$rep_SRs[,1,])

ggpairs.PAsso(object = PAsso_1, colour="blue")
ggpairs.PAsso(object = PAsso_adv1, colour="blue")
ggpairs.PAsso(object = PAsso_2, colour="blue")
ggpairs.PAsso(object = PAsso_5v, colour="blue")


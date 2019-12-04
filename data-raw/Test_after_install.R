
library(MASS)
library(faraway)
library(parcor)
library(parasol); library(tidyverse)
# import data -------------------------------------------------------------
data(nes96clean)

# parasol::nes96clean
# Bivariate analysis of VOTE and PID (Kendall's tao) ----------------------
y1<- nes96clean$vote.num
y2<- nes96clean$PID
X<- as.matrix(nes96clean[c("income.num", "age", "edu.year")])

# marginal association (Kendall's tau)
tau <- cor(y1, y2, method = "kendall")
# standard error from bootstrap
tau.sd.boot <- sd(sapply(1:2000, function(b){
  index<- sample(nrow(nes96), replace=T)
  cor(y1[index], y2[index], method = "kendall")
}))
tau; tau.sd.boot

# Test "Pcor_SR" function: Partial Association by surrogate residuals ------------------------------
# regression models
fit.vote<- glm(y1 ~ X, family = binomial(link ="probit"))
fit.PID<- polr(as.factor(y2)~ X, method="probit")

Pcor_test1 <- Pcor_SR(models=list(fit.vote, fit.PID),
                       asso_type = c("PartialAsso"),
                       cor_method = c("kendall"), rep_num=10)

# Partial association coefficients (Parts of Table 7 in paper)
Pcor_test1$Mcor_phi

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
fit.DoleLR<- polr(ClinLR ~ income.num + age + edu.year,
                  data = dta, method="probit")

# surrogate residual
Pcor_5v <- Pcor_SR(
  models=list(fit.vote, fit.PID, fit.selfLR, fit.ClinLR, fit.DoleLR),
  asso_type = c("PartialAsso"),
  cor_method = c("kendall"), rep_num=30)

Pcor_5v$Mcor_phi # (Right part of Table 7 in paper)
# Pcor_5v$Mboot_phi


# Partial Regression plot matrix ------------------------------------------

colnames(Pcor_5v$SRs) <- c("VOTE", "PID", "selfLR", "ClinLR", "DoleLR")
ggpairs.resid(resid_Mat = as.data.frame(Pcor_5v$SRs), colour="blue")


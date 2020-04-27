
library(MASS)
library(PAsso)
library(tidyverse)
# library(matrixStats)
library(progress)
library(profvis)

library(PAsso)
# import data -------------------------------------------------------------
data("nes2016_pre")
nes2016 <- nes2016_pre
summary(nes2016)

nes2016 <- nes2016 %>% rename(vote.num=Prevote.num)

system.time(
  cor_matrix1 <- cor(nes2016[, c("vote.num","PID")], use = "pairwise.complete.obs", method = "pearson")
)

system.time(
  cor_matrix2 <- cor(nes2016[, c("vote.num","PID")], use = "pairwise.complete.obs", method = "kendall")
)
library(pcaPP)
system.time(
  cor_matrix3 <- cor.fk(nes2016[, c("vote.num","PID")])
)
cor_matrix2; cor_matrix3

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

# Partial association coefficients (Parts of Table 7 in paper)
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)

# Test jittering
system.time(PAsso_adv1_jit <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                association = c("partial"),
                                method = c("kendall"),
                                resids.method = "jitter")
)
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1_jit, digits = 3)

# "PAsso" function: The simple way, input response and confounders only ----------------------------
PAsso_1 <- PAsso(responses = c("vote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = nes2016
                   # association = c("partial"),
                   # models = c("probit", "probit"),
                   # method = c("kendall"),
                   # resids.method = "latent", fitted.models = NULL,
                   # rep_num = 30
                )

print(PAsso_1, 5)
summary(PAsso_1, 4)
plot(PAsso_1)

# "PAsso" function: input three responses ----------------------------
PAsso_2 <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016,
                models = c("probit", "probit", "probit"),
                association = c("partial"), method = c("kendall"),
                resids.method = "latent")

profvis({
  PAsso_2 <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = nes2016,
                   models = c("probit", "probit", "probit"),
                   association = c("partial"), method = c("kendall"),
                   resids.method = "latent",rep_num = 2)
})
# Compare marginal correlation and partial correlation.
print(PAsso_2, digits=4)
summary(PAsso_2, digits=4)
plot(PAsso_2)

PAsso_2_jit <- PAsso(responses = c("vote.num", "PID", "selfLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 uni.model = "probit",
                 association = c("partial"), method = c("kendall"),
                 resids.method = "jitter")
print(PAsso_2, digits=4)
summary(PAsso_2_jit, digits=4)

# "PAsso" function: input three responses ----------------------------

# Pcor_SR.test function: Conduct inference based on object of "PartialCor" class ----------------------------
library(progress); #library(doParallel)

system.time(Pcor_SR_test1 <- test(object = PAsso_2, boot_SE=100, H0=0, parallel=F))

print(Pcor_SR_test1, digits=3)
summary(PAsso_1, digits=4)

# Pcor_SR.test function: Test by parallel ----------------------------

library(doParallel); library(progress)

numCores <- detectCores()  # Not too aggressive!
cl <- makeCluster(numCores)
# registerDoSNOW(cl) # on Mac or Linux
registerDoParallel(cl) # Win

system.time(Pcor_SR_test2 <- test(object = PAsso_2, boot_SE=200, H0=0, parallel=TRUE))

profvis({Pcor_SR_test2 <- test(object = PAsso_2, boot_SE=20, H0=0, parallel=TRUE)
})
summaryRprof(Pcor_SR_test2 <- test(object = PAsso_2, boot_SE=20, H0=0, parallel=TRUE))

print(Pcor_SR_test2, digits=3)
summary(Pcor_SR_test2)
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
print(PAsso_adv_5v, digits=3)
summary(PAsso_adv_5v, digits=4)

PAsso_5v <- PAsso(responses = c("vote.num", "PID", "selfLR", "ClinLR", "TrumpLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 association = c("partial"),
                 method = c("kendall"), resids.method = "latent", rep_num=30)

PAsso_5v # (Right part of Table 7 in paper)
summary(PAsso_5v,3)

system.time(PAsso_5v_test <- test(object = PAsso_5v, boot_SE=200, H0=0, parallel=FALSE))

print(PAsso_5v_test, 3)

# Partial Regression plot matrix ------------------------------------------
plot(PAsso_1)
plot(x = PAsso_adv1)
plot(x = PAsso_2)
plot(x = PAsso_5v)
# plot(object = PAsso_5v_test, colour="blue")


# check.model function -----------------------------------------------------
library(PAsso); library(ggplot2)
data("nes2016_pre")
PAsso_2 <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016_pre,
                models = c("probit", "probit", "probit"),
                association = c("partial"), method = c("kendall"),
                resids.method = "latent")

summary(PAsso_2, digits = 4)

check_qq <- check.model(object = PAsso_2, color="blue", output = "qq")

check_fitted <- check.model(object = PAsso_2, output = "fitted")

check_covar <- check.model(object = PAsso_2, output = "covariate")

check_qq <- check.model(object = PAsso_5v, color="blue", output = "qq")

check_fitted <- check.model(object = PAsso_5v, output = "fitted")

check_covar <- check.model(object = PAsso_5v, output = "covariate")


## general association measure and 3-D plot for VOTE and PID
library("copula")
library("plotly")
# library("copBasic")

dim(PAsso_2$rep_SRs)
testPlots <- plot3D(PAsso_2)
testPlots$plot_1


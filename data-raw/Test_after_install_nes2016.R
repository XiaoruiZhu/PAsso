
library(MASS)
# library(tidyverse)
# library(progress)

library(PAsso)
# import data -------------------------------------------------------------
data("nes2016")
summary(nes2016)

# Partial association between "Prevote.num" and "PID" after adjusting ----------------------------
# "income.num", "age", "edu.year"
# "PAsso" function: Only need input responses, adjustments, data
# Other default arguments are displayed below as well
system.time(PAsso_1 <- PAsso(responses = c("Prevote.num", "PID"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 uni.model = "probit",
                 # association = c("partial"),
                 # models = c("probit", "probit"),
                 method = c("kendall"),
                 # resids.type = "surrogate", fitted.models = NULL,
                 n_draws = 30
                )
)
# Print the partial association matrix only
print(PAsso_1, 5)

# Provide partial association matrix, marginal association matrix, and summary of models' coefficients
summary(PAsso_1, 4)

# Plot partial association regression plot: residuals
plot(PAsso_1)

# Retrieve residuals
test_resids <- residuals(PAsso_1, draw = 1)
head(test_resids)
dim(test_resids)
# test_wolf <- t(copBasic::wolfCOP(para = data.frame(test_resids), as.sample = TRUE))

head(residuals(object = PAsso_1, draw_id=2))
# "PAsso" function: input three responses ----------------------------
PAsso_2 <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016,
                uni.model <- "probit",
                method = c("kendall"),
                # models = c("probit", "probit", "probit"),
                # association = c("partial")
                resids.type = "surrogate")

# Compare marginal correlation and partial correlation.
summary(PAsso_2, digits=4)
plot(PAsso_2)

PAsso_2_jit <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 uni.model = "logit",
                 association = c("partial"), method = c("kendall"),
                 resids.type = "surrogate", jitter = "uniform")
summary(PAsso_2_jit, digits=4)

# test function: Conduct inference based on object of "PAsso.test" class ----------------------------
library(progress); #library(doParallel)

system.time(Pcor_SR_test1 <- test(object = PAsso_2, boot_SE=100, H0=0, parallel=F))
print(Pcor_SR_test1, digits=3)

# test function: Test parallel ----------------------------
library(doParallel); library(progress)

numCores <- detectCores()  # Not too aggressive!
cl <- makeCluster(numCores)
# registerDoSNOW(cl) # on Mac or Linux
registerDoParallel(cl) # Win

system.time(Pcor_SR_test2 <- test(object = PAsso_2, boot_SE=200, H0=0, parallel=TRUE))

print(Pcor_SR_test2, digits=5)
stopCluster(cl)

# multivariate analysis (5 variables) --------------------------------------------------------------------
PAsso_5v <- PAsso(responses = c("Prevote.num", "PID", "selfLR", "ClinLR", "TrumpLR"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 association = c("partial"),
                 method = c("kendall"), resids.type = "surrogate")

print(PAsso_5v, 4) # (Right part of Table 7 in paper)
summary(PAsso_5v,4)

# Partial Regression plot matrix ------------------------------------------
plot(PAsso_1)
plot(x = PAsso_2)
plot(x = PAsso_5v)

# diagnostic.plot function -----------------------------------------------------

check_qq <- diagnostic.plot(object = PAsso_2, output = "qq")

check_fitted <- diagnostic.plot(object = PAsso_2, output = "fitted")

check_covar <- diagnostic.plot(object = PAsso_2, output = "covariate")

check_qq <- diagnostic.plot(object = PAsso_5v, color="blue", output = "qq")

check_fitted <- diagnostic.plot(object = PAsso_5v, output = "fitted")

check_covar <- diagnostic.plot(object = PAsso_5v, output = "covariate")

## general association measure and 3-D plot for VOTE and PID
library("copula")
library("plotly")

testPlots <- plot3D(PAsso_2)
testPlots$plot_1
# Above result need to be opened in browser through "Viewer" tap!

###########################################################
# Advanced user's guide
###########################################################

# "PAsso" advanced using of the function: Input a few models directly ------------------------------

fit.vote<- glm(Prevote.num ~ income.num+ age + edu.year, data = nes2016,
               family = binomial(link = "probit"))
fit.PID<- polr(as.factor(PID) ~ income.num+age+edu.year, data = nes2016,
               method="probit", Hess = TRUE)

system.time(PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                association = c("partial"),
                                method = c("kendall"),
                                resids.type = "surrogate")
)

# Partial association coefficients (Parts of Table 7 in paper)
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)

# Test jittering
system.time(PAsso_adv1_jit <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                    association = c("partial"),
                                    method = c("kendall"),
                                    resids.type = "surrogate", jitter = "uniform",
                                    jitter.uniform.scale = "response")
)
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1_jit, digits = 3)


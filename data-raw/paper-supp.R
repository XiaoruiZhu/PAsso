# setwd("C:/Users/Shaobo/Dropbox/Research/With Dungang/Framework/Rcode")

library(MASS)
library(faraway)
library(parcor)
library(doParallel)
source("data-raw/paper_supp_functions.R")
source("data-raw/paper_supp_data_preproc.R")
## import data
data("nes96")

# Frequency table (table 3)
table(nes96$PID, nes96$vote)

# data pre-processing
nes96clean <- nes.clean()

### Bivariate analysis of VOTE and PID ###
y1 <- nes96clean$vote.num
y2 <- nes96clean$PID
X <- as.matrix(nes96clean[c("income.num", "age", "edu.year")])

# marginal association (Kendall's tau)
tau <- cor(y1, y2, method = "kendall")
# standard error from bootstrap
tau.sd.boot <- sd(sapply(1:2000, function(b) {
  index <- sample(nrow(nes96), replace = T)
  cor(y1[index], y2[index], method = "kendall")
}))

# regression models
fit.vote <- glm(y1 ~ X, family = binomial(link = "probit"))
fit.PID <- polr(as.factor(y2) ~ X, method = "probit")

Pcor_test1 <- Pcor_SR2(
  models = list(fit.vote, fit.PID),
  asso_type = c("PartialAsso"),
  cor_method = c("kendall"), rep_num = 10
)
Pcor_test1$Mcor_phi
Pcor_test1$Mboot_phi
dim(Pcor_test1$boot_SRs)

# partial association (Kendall's tau)
set.seed(2020)
pcor.result <- Pcor_SR(y1, y2, X, model1 = fit.vote, model2 = fit.PID, link1 = "probit", link2 = "probit", n.avg = 30, method = "kendall")
pcor_tau <- pcor.result$phi
pcor_tau
# standard error from bootstrap (this takes 10-20 minutes)
pcor_tau.sd.boot <- Pcor_SR_boot(y1, y2, X, B = 2000, link1 = "probit", link2 = "probit", n.avg = 30, method = "kendall")

# visualization (Figure 7)
r_vote <- pcor.result$sr1[, 1]
r_PID <- pcor.result$sr2[, 1]
ind.clin <- which(nes96$vote == "Clinton")
plot(r_PID[ind.clin], r_vote[ind.clin],
  xlab = expression(r[PID]),
  ylab = expression(r[VOTE]), col = "darkblue", ylim = c(-4, 4),
  xlim = c(-4, 4), pch = 1, cex.lab = 1.2
)
points(r_PID[-ind.clin], r_vote[-ind.clin], col = "red", pch = 4)
abline(v = 0, h = 0, lty = c(2, 2))
abline(v = -1, lty = 3)
abline(v = 1, lty = 3)
legend("topleft",
  legend = c("Clinton", "Dole"),
  col = c("darkblue", "red"), pch = c(1, 4)
)
smoothingSpline.Q <- smooth.spline(r_PID, r_vote, spar = 1.3)
lines(smoothingSpline.Q, col = 1, lwd = 2, lty = 2)

# histogram of projection (Figure 8)
proj.var1 <- cbind(r_PID, r_vote) %*% c(1, sqrt(3)) / 2
hist(proj.var1, breaks = 30, probability = T, xlab = "", main = "")
lines(density(proj.var1), lwd = 2, col = 2, lty = 2)

# Q-Q plot (Figure 9)
qqplot(rnorm(length(r_PID)), r_PID,
  xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5),
  xlab = "Theoretical Quantiles",
  ylab = expression(paste("Sample Quantiles (", r[PID], ")"))
)
abline(0, 1, lty = 2, col = 2)
qqplot(rnorm(length(r_vote)), r_vote,
  xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5),
  xlab = "Theoretical Quantiles",
  ylab = expression(paste("Sample Quantiles (", r[VOTE], ")"))
)
abline(0, 1, lty = 2, col = 2)

### multivariate analysis (5 variables) ###
y3 <- nes96clean$selfLR
y4 <- nes96clean$ClinLR
y5 <- nes96clean$DoleLR
# marginal correlation
dta <- dplyr::select(
  nes96clean, vote.num, PID, selfLR, ClinLR, DoleLR,
  income.num, age, edu.year
) %>%
  dplyr::mutate(
    vote.num = as.factor(vote.num), PID = as.factor(PID),
    selfLR = as.factor(selfLR), ClinLR = as.factor(ClinLR), DoleLR = as.factor(DoleLR)
  )
head(dta)

Mcor_tau <- cor(
  dplyr::select(
    nes96clean, vote.num,
    PID, selfLR, ClinLR, DoleLR
  ),
  method = "kendall"
)

# regression model
fit.selfLR <- polr(selfLR ~ income.num + age + edu.year,
  data = dta, method = "probit"
)
fit.ClinLR <- polr(ClinLR ~ income.num + age + edu.year,
  data = dta, method = "probit"
)
fit.DoleLR <- polr(ClinLR ~ income.num + age + edu.year,
  data = dta, method = "probit"
)

# surrogate residual
Pcor_5v <- Pcor_SR2(
  models = list(fit.vote, fit.PID, fit.selfLR, fit.ClinLR, fit.DoleLR),
  asso_type = c("PartialAsso"),
  cor_method = c("kendall"), rep_num = 30
)

Pcor_5v$Mcor_phi
# Pcor_5v$Mboot_phi
dim(Pcor_5v$boot_SRs)
colnames(Pcor_5v$SRs) <- c("VOTE", "PID", "selfLR", "ClinLR", "DoleLR")
ggpairs.resid(resid_Mat = as.data.frame(Pcor_5v$SRs), colour = "blue")


pcor.result <- Pcor_SR(
  y4, y5, X,
  model1 = fit.ClinLR, model2 = fit.DoleLR,
  link1 = "probit", link2 = "probit", n.avg = 30, method = "kendall"
)
pcor_tau <- pcor.result$phi
pcor_tau

# Old method to draw plot matrix ------------------------------------------

M <- 30
set.seed(2020)
r_vote <- SR(y1, X, model = fit.vote, ndraw = M)
r_PID <- SR(y2, X, model = fit.PID, ndraw = M)
r_selfLR <- SR(y3, X, model = fit.selfLR, ndraw = M)
r_ClinLR <- SR(y4, X, model = fit.ClinLR, ndraw = M)
r_DoleLR <- SR(y5, X, model = fit.DoleLR, ndraw = M)

# association matrix
phi_tau5 <- matrix(apply(sapply(1:M, function(j) cor(cbind(r_vote[, j], r_PID[, j], r_selfLR[, j], r_ClinLR[, j], r_DoleLR[, j]), method = "kendall")), 1, mean), 5, 5)
# bootstrap standard error (parallel computing)
set.seed(2020)
B <- 2000
n <- nrow(nes96clean)
bootind <- sample(n, n * B, replace = T)
bootdat <- lapply(1:B, function(xx) nes96clean[bootind[((xx - 1) * n):(xx * n)], ])
ncore <- detectCores() - 2
cl <- makeCluster(ncore)
registerDoParallel(cl)
Pcor.est <- foreach(i = 1:B, .combine = "rbind", .export = c("polr")) %dopar% {
  y1.b <- bootdat[[i]]$vote.num
  y2.b <- bootdat[[i]]$PID
  y3.b <- bootdat[[i]]$selfLR
  y4.b <- bootdat[[i]]$ClinLR
  y5.b <- bootdat[[i]]$DoleLR
  x.b <- bootdat[[i]][c("income.num", "age", "edu.year")]

  fit1 <- glm(y1.b ~ as.matrix(x.b), family = binomial(link = "probit"))
  fit2 <- polr(as.factor(y2.b) ~ as.matrix(x.b), method = "probit")
  fit3 <- polr(as.factor(y3.b) ~ as.matrix(x.b), method = "probit")
  fit4 <- polr(as.factor(y4.b) ~ as.matrix(x.b), method = "probit")
  fit5 <- polr(as.factor(y5.b) ~ as.matrix(x.b), method = "probit")

  r1 <- SR(y1.b, x.b, model = fit1, ndraw = M)
  r2 <- SR(y2.b, x.b, model = fit2, ndraw = M)
  r3 <- SR(y3.b, x.b, model = fit3, ndraw = M)
  r4 <- SR(y4.b, x.b, model = fit4, ndraw = M)
  r5 <- SR(y5.b, x.b, model = fit5, ndraw = M)

  phi_all <- sapply(1:M, function(j) cor(cbind(r1[, j], r2[, j], r3[, j], r4[, j], r5[, j]), method = "kendall"))
  apply(phi_all, 1, mean)
}
stopCluster(cl)
phi_sd <- matrix(apply(Pcor.est, 2, sd), 5, 5)


# pairwise scatter plot figure
pairs(data.frame(VOTE = r_vote[, 1], PID = r_PID[, 1], selfLR = r_selfLR[, 1], ClinLR = r_ClinLR[, 1], DoleLR = r_DoleLR[, 1]),
  panel = panel.smooth, span = 0.8, lower.panel = NULL, pch = "."
)

######################################

## general association measure and 3-D plot for VOTE and PID
library("copula")
library("plotly")
library("copBasic")

resi <- cbind(r_vote, r_PID)
empC <- C.n(pobs(resi), resi)
empFG1 <- pobs(resi)[, 1] * pobs(resi)[, 2]

WolfCop <- wolfCOP(para = data.frame(r_vote, r_PID), as.sample = TRUE)
WolfCop

## 3-D copula plot
v1 <- v2 <- seq(0, 1, length.out = 100)
aa <- matrix(0, length(v1), length(v1))
for (i in 1:length(v1)) {
  for (j in 1:length(v1)) {
    aa[i, j] <- C.n(t(as.matrix(c(v1[i], v2[j]))), resi) - v1[i] * v2[j]
  }
}

options(Viewer = NULL)
plot_ly(x = v1, y = v2, z = 12 * aa) %>%
  add_surface() %>%
  layout(scene = list(xaxis = list(title = "u"), yaxis = list(title = "v"), zaxis = list(title = "12(C(u,v)-uv)")))

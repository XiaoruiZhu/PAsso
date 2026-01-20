# various internal functions


##### function to obtain probability matrix of adjacent-categorical model ####
p_adj_cate <- function(Z) {
  k <- ncol(Z)
  p1_pj <- p1_pj_inv <- Z
  ZZ <- 0
  for (j in 1:k) {
    ZZ <- ZZ + Z[, j]
    p1_pj[, j] <- exp(ZZ)
    p1_pj_inv[, j] <- 1 / p1_pj[, j]
  }
  p1 <- 1 / (rowSums(p1_pj_inv) + 1)
  pj <- p1_pj
  for (j in 1:k) {
    pj[, j] <- p1 / p1_pj[, j]
  }
  cbind(1 - rowSums(pj), pj)
}

#### surrogate residual ####

SR.acat <- function(y, X, alpha, beta, ndraw = 1) {
  y <- Y1
  X <- X
  alpha <- alpha1
  beta <- beta1
  ndraw <- 1

  y <- as.integer(y)
  n <- length(y)
  z <- sapply(alpha[2:(length(alpha) - 1)], function(a) a + as.matrix(X) %*% beta)
  p.acat <- p_adj_cate(z)
  F.acat <- t(apply(p.acat, 1, cumsum))
  F.acat <- cbind(0, F.acat)
  if (min(y) == 0) {
    S <- sapply(1:n, function(k) runif(ndraw, F.acat[k, y[k] + 1], F.acat[k, y[k] + 2]))
  } else {
    S <- sapply(1:n, function(k) runif(ndraw, F.acat[k, y[k]], F.acat[k, y[k] + 1]))
  }
  S - 1 / 2
}


##### Probability based residual #####

PR.acat <- function(y, X, alpha, beta) {
  z <- sapply(alpha[2:(length(alpha) - 1)], function(a) a + as.matrix(X) %*% beta)
  p.acat <- p_adj_cate(z)
  n <- length(y)
  # calculate probability based residual based on fitted value (prbabilities)
  pyej <- p.acat[cbind(1:n, y)]
  pysj <- sapply(1:n, function(x) sum(p.acat[x, 1:y[x]])) - pyej
  PR <- -1 + 2 * pysj + pyej
  # PR<- 1-2*pysj-pyej

  return(PR)
} # end of the function of probability-based residual


##### Generalized residual (Franses and Paap 2001) added for revision #####
GR.acat <- function(y, X, alpha, beta) {
  y <- as.integer(y)
  n <- length(y)
  z <- sapply(alpha[2:(length(alpha) - 1)], function(a) a + as.matrix(X) %*% beta)
  p.acat <- p_adj_cate(z)
  F.acat <- t(apply(p.acat, 1, cumsum))
  F.acat <- cbind(0, F.acat)

  Pj <- sapply(1:n, function(x) p.acat[x, y[x]])
  Fj <- sapply(1:n, function(x) F.acat[x, y[x] + 1])
  Fj_1 <- sapply(1:n, function(x) F.acat[x, y[x]])

  fj <- dnorm(qnorm(Fj))
  fj_1 <- dnorm(qnorm(Fj_1))
  return((fj_1 - fj) / Pj)
} # end of function

##### deviance residual (-2*loglik) #####
DR.acat <- function(y, X, alpha, beta) {
  n <- length(y)
  z <- sapply(alpha[2:(length(alpha) - 1)], function(a) a + as.matrix(X) %*% beta)
  p.acat <- p_adj_cate(z)
  return(-2 * log(p.acat[cbind(1:n, y)]))
}


# Running models ----------------------------------------------------------
library(VGAM)

load(file = "data/df_ParA.RData")
summary(df_ParA$data)

Y1 <- df_ParA$data$Y1
Y2 <- df_ParA$data$Y2
X <- df_ParA$data$X

# fit.adj1<- vglm(Y1~X, family=acat(link = "logitlink", reverse=TRUE, parallel=TRUE))

fit.adj1 <- vglm(Y1 ~ X, family = acat(reverse = TRUE, parallel = TRUE))


fit1 <- VGAM::vglm(Y1 ~ X,
  family =
    VGAM::cumulative(
      link = "logit", reverse = TRUE,
      parallel = TRUE
    )
)
fit.adj2 <- vglm(Y2 ~ X, family = acat(reverse = TRUE, parallel = TRUE))

fit1@family@vfamily[1]
fit1@misc$reverse
head(fit.adj1@x[, -1])
head(fit.adj1@fitted.values)

head(t(apply(fit.adj1@fitted.values, 1, cumsum)))

head(class(as.matrix(fit.adj1@x)[, -1]))

summary(Y1)

rbind(
  coef(fit.adj1), coef(fit1),
  c(df_ParA$alpha1[c(-1, -7)], df_ParA$beta1)
)
# rbind(coef(fit.adj1), coef(fit1)[-length(coef(fit1))],
#       c(df_ParA$alpha1[c(-1, -7)], df_ParA$beta1))
# + I(X ^ 2)
rbind(coef(fit.adj2), c(df_ParA$alpha2[c(-1, -5)], df_ParA$beta2))


SR1 <- SR.acat(Y1, X, alpha = alpha1, beta = beta1)
SR2 <- SR.acat(Y2, X, alpha = alpha2, beta = beta2)

## obtain SBC residuals (Li and Shepherd 2012 JASA/Biometrika)
PR1 <- PR.acat(Y1, X, alpha = alpha1, beta = beta1)
PR2 <- PR.acat(Y2, X, alpha = alpha2, beta = beta2)

## obtain generalized residuals (Franses and Paap 2001 book)
GR1 <- GR.acat(Y1, X, alpha = alpha1, beta = beta1)
GR2 <- GR.acat(Y2, X, alpha = alpha2, beta = beta2)

## obtain deviance residuals
DR1 <- DR.acat(Y1, X, alpha = alpha1, beta = beta1)
DR2 <- DR.acat(Y2, X, alpha = alpha2, beta = beta2)


# Use generate_residuals  -------------------------------------------------
# Fit cumulative link models

load(file = "data/df_ParA.RData")
summary(df_ParA$data)

fit_clm1 <- VGAM::vglm(Y1 ~ X,
  data = df_ParA$data, family =
    VGAM::cumulative(link = "logit", reverse = TRUE, parallel = TRUE)
)
fit_clm2 <- VGAM::vglm(Y2 ~ X,
  data = df_ParA$data, family =
    VGAM::cumulative(link = "logit", reverse = TRUE, parallel = TRUE)
)
coef(fit_clm1)

# fit_clm1 <- VGAM::vglm(Y1 ~ X, data = df_ParA$data, family =
#                          VGAM::acat(reverse=TRUE, parallel=TRUE))
# fit_clm2 <- VGAM::vglm(Y2 ~ X, data = df_ParA$data, family =
#                          VGAM::acat(reverse=TRUE, parallel=TRUE))
# coef(fit_clm1)
# fit_clm1@family@infos()$link

# SR1 <- generate_residuals(fit_clm1, method = "latent",boot_id = NULL)
# SR2 <- generate_residuals(fit_clm2, method = "latent", boot_id = NULL)

SR1 <- generate_residuals(fit_clm1,
  method = c("jitter"),
  jitter.scale = c("response"),
  boot_id = NULL
)
SR2 <- generate_residuals(fit_clm2,
  method = c("jitter"),
  jitter.scale = c("response"),
  boot_id = NULL
)

## obtain SBC residuals (Li and Shepherd 2012 JASA/Biometrika)
PR1 <- generate_residuals(fit_clm1, method = "Sign", boot_id = NULL)
PR2 <- generate_residuals(fit_clm2, method = "Sign", boot_id = NULL)

## obtain generalized residuals (Franses and Paap 2001 book)
GR1 <- generate_residuals(fit_clm1, method = "General", boot_id = NULL)
GR2 <- generate_residuals(fit_clm2, method = "General", boot_id = NULL)

## obtain deviance residuals
DR1 <- generate_residuals(fit_clm1, method = "Deviance", boot_id = NULL)
DR2 <- generate_residuals(fit_clm2, method = "Deviance", boot_id = NULL)


## visualize residual vs. residual
par(mfrow = c(2, 2))
par(mar = c(4, 4.8, 2.5, 1.5))

plot(PR1, PR2,
  pch = ".", main = "Sign-based Residuals",
  xlab = expression(paste(R[1]^"ALT")),
  ylab = expression(paste(R[2]^"ALT"))
)
plot(GR1, GR2,
  pch = ".", main = "Generalized Residuals",
  xlab = expression(paste(R[1]^"ALT")),
  ylab = expression(paste(R[2]^"ALT")), xlim = c(-4, 4), ylim = c(-4, 4)
)
plot(DR1, DR2,
  pch = ".", main = "Deviance Residuals",
  xlab = expression(paste(R[1]^"ALT")),
  ylab = expression(paste(R[2]^"ALT"))
)
plot(SR1, SR2,
  pch = ".", main = "Surrogate Residuals", xaxt = "n", yaxt = "n",
  xlab = expression(R[1]), ylab = expression(R[2]),
  xlim = c(-1 / 2, 1 / 2), ylim = c(-1 / 2, 1 / 2)
)
axis(1, at = seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))
axis(2, at = seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))

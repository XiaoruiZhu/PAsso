### functions needed for real data analysis ###

######################################################
#### functions for generating surrogate residual  ####
######################################################
library(truncdist)

## redefine gumbel distribution ##
# For log-log link
pgumbel <- function(q, location = 0, scale = 1) {
  q <- (q - location) / scale
  exp(-exp(-q))
}
qgumbel <- function(p, location = 0, scale = 1) {
  -scale * log(-log(p)) + location
}
rgumbel <- function(n, location = 0, scale = 1) {
  qgumbel(runif(n, min = 0, max = 1), location = location, scale = scale)
}

# For complimentary log-log link
pGumbel <- function(q, location = 0, scale = 1) {
  q <- (q - location) / scale
  1 - exp(-exp(q))
}
qGumbel <- function(p, location = 0, scale = 1) {
  scale * log(-log(1 - p)) + location
}
rGumbel <- function(n, location = 0, scale = 1) {
  qGumbel(runif(n, min = 0, max = 1), location = location, scale = scale)
}

## modified rtrunc()
rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
  qtrunc(runif(n, min = 0, max = 1), spec, a = a, b = b, ...)
}
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  G.a <- G(a, ...)
  G.b <- G(b, ...)
  pmin(pmax(a, Gin(G(a, ...) + p * (G(b, ...) - G(a, ...)), ...)), b)
}
#####################################################

SR <- function(y, X, model, alpha = NULL, beta = NULL, ndraw = 1, dist = "norm") {
  y <- as.numeric(y)
  if (min(y) == 0) y <- y + 1

  if (missing(model)) {
    if ((is.null(alpha) | is.null(beta))) {
      stop("Alpha and beta must be provided if model is not provided.")
    } else {
      # alpha<- c(-Inf,alpha,Inf)
      if (missing(X) | sum(X^2) == 0) {
        Lhat <- 0
      } else {
        Lhat <- as.vector(as.matrix(X) %*% beta)
      }
    }
  } else {
    if (length(table(y)) == 2) {
      alpha <- c(-Inf, -model$coeff[1], Inf)
      beta <- model$coeff[-1]
    } else {
      alpha <- c(-Inf, model$zeta, Inf)
      beta <- model$coeff
    }

    if (missing(X) | sum(X^2) == 0) {
      Lhat <- 0
    } else {
      Lhat <- as.vector(as.matrix(X) %*% beta)
    }
  }

  n <- length(y)
  alpha.mat <- matrix(rep(alpha, n), nrow = n, byrow = TRUE)
  lower <- alpha.mat[cbind(1:n, y)]
  upper <- alpha.mat[cbind(1:n, y + 1)]

  # S<- matrix(rtnorm(n*ndraw, Lhat, 1, lower = lower, upper = upper), n, ndraw)
  if (dist == "norm") {
    S <- matrix(rtrunc(n * ndraw, spec = dist, a = lower, b = upper, mean = Lhat), n, ndraw) # n*ndraw numbers with every n matches a, b, mean
  } else {
    if (dist %in% c("logis", "Gumbel", "gumbel")) {
      S <- matrix(rtrunc(n * ndraw, spec = dist, a = lower, b = upper, loc = Lhat), n, ndraw)
    } else {
      stop("Distribution not supported! Supported distributions currently are 'norm', 'logis', 'Gumbel'.")
    }
  }

  r <- S - Lhat
  # plot(z1,r_x)
  return(r)
} # end of the function
################################


## function to obtain \hat{\phi}}
Pcor_SR <- function(y1, y2, x, model1, model2, alpha1 = NULL, beta1 = NULL, alpha2 = NULL, beta2 = NULL, link1 = "probit", link2 = "probit", n.avg = 100, method = "kendall") {
  if (link1 == "probit") distribution1 <- "norm"
  if (link1 == "logit" | link1 == "logistic") {
    link1 <- "logistic"
    distribution1 <- "logis"
  }
  if (link1 == "cloglog") distribution1 <- "Gumbel"
  if (link1 == "loglog") distribution1 <- "gumbel"
  if (link1 == "cauchit") distribution1 <- "cauchy"

  if (link2 == "probit") distribution2 <- "norm"
  if (link2 == "logit" | link2 == "logistic") {
    link2 <- "logistic"
    distribution2 <- "logis"
  }
  if (link2 == "cloglog") distribution2 <- "Gumbel"
  if (link2 == "loglog") distribution2 <- "gumbel"
  if (link2 == "cauchit") distribution2 <- "cauchy"

  n <- length(y1)
  dat <- data.frame(y1, y2, x)

  if (!is.null(alpha1) & !is.null(alpha2) & !is.null(beta1) & !is.null(beta2) & missing(model1) & missing(model2)) {
    sr1 <- SR(y1, x, alpha = alpha1, beta = beta1, ndraw = n.avg, dist = distribution1)
    sr2 <- SR(y2, x, alpha = alpha2, beta = beta2, ndraw = n.avg, dist = distribution2)
    fit1 <- fit2 <- NULL
  } else {
    if (is.null(c(alpha1, alpha2)) & is.null(c(beta1, beta2)) & !missing(model1) & !missing(model2)) {
      sr1 <- SR(y1, x, model = model1, ndraw = n.avg, dist = distribution1)
      sr2 <- SR(y2, x, model = model2, ndraw = n.avg, dist = distribution2)
    } else {
      stop("Either models or alpha's and beta's should be supplied")
    }
  }

  # calculate \hat{\phi}
  if (method == "kendall") {
    cor.est <- sapply(1:n.avg, function(yy) cor(sr1[, yy], sr2[, yy], method = "kendall"))
  }
  if (method == "pearson") {
    cor.est <- sapply(1:n.avg, function(yy) cor(sr1[, yy], sr2[, yy]))
  }
  if (method == "wolfsigma") {
    cor.est <- t(sapply(1:n.avg, function(yy) wolfCOP(para = data.frame(cbind(sr1[, yy], sr2[, yy])), as.sample = TRUE)))
  }
  return(list(phi = mean(cor.est), phi_all = cor.est, dist1 = distribution1, dist2 = distribution2, sr1 = sr1, sr2 = sr2))
} ## end of function
########################

## function to perform bootstrap inference
Pcor_SR_boot <- function(y1, y2, x, B = 2000, link1 = "probit", link2 = "probit", n.avg = 100, method = "kendall", H0 = 0) {
  pcor.est <- rep(NA, B)
  if (length(unique(y1)) > 2 & length(unique(y2)) > 2) {
    for (i in 1:B) {
      tryCatch(
        {
          index <- sample(length(y1), replace = T)
          y1.b <- y1[index]
          y2.b <- y2[index]
          if (is.vector(x)) x.b <- x[index] else x.b <- x[index, ]
          fit1 <- polr(as.factor(y1.b) ~ x.b, method = link1)
          fit2 <- polr(as.factor(y2.b) ~ x.b, method = link2)
          pcor.est[i] <- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1 = link1, link2 = link2, n.avg = n.avg, method = method)$phi
        },
        error = function(e) {}
      )
    }
  } else {
    if (length(unique(y1)) == 2 & length(unique(y2)) > 2) {
      for (i in 1:B) {
        tryCatch(
          {
            index <- sample(length(y1), replace = T)
            y1.b <- y1[index]
            y2.b <- y2[index]
            if (is.vector(x)) x.b <- x[index] else x.b <- x[index, ]
            fit1 <- glm(y1.b ~ x.b, family = binomial(link = link1))
            fit2 <- polr(as.factor(y2.b) ~ x.b, method = link2)
            pcor.est[i] <- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1 = link1, link2 = link2, n.avg = n.avg, method = method)$phi
          },
          error = function(e) {}
        )
      }
    } else {
      for (i in 1:B) {
        tryCatch(
          {
            index <- sample(length(y1), replace = T)
            y1.b <- y1[index]
            y2.b <- y2[index]
            if (is.vector(x)) x.b <- x[index] else x.b <- x[index, ]
            fit1 <- glm(y1.b ~ x.b, family = binomial(link = link1))
            fit2 <- glm(y2.b ~ x.b, family = binomial(link = link2))
            pcor.est[i] <- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1 = link1, link2 = link2, n.avg = n.avg, method = method)$phi
          },
          error = function(e) {}
        )
      }
    }
  }
  Bstar <- sum(!is.na(pcor.est))
  pcor.est1 <- pcor.est[!is.na(pcor.est)]

  if (H0 == 0) {
    pval <- 2 * min(mean(pcor.est1 > 0), mean(pcor.est1 < 0))
  } else {
    pval <- 2 * min(mean(pcor.est1 < H0), mean(pcor.est1 > -1 * H0), 0.5)
  }

  return(list(est = mean(pcor.est1), sd = sd(pcor.est1), pval = pval, CI = quantile(pcor.est1, probs = c(0.025, 0.975)), Bcount = Bstar))
}

#' Surrogate Residual Generating
#'
#' This function simulates surrogate residuals for ordinal variables by adjusting
#' a set of covariates X using cumulative link models. Various links are considered.
#' This version incorporates binary response from glm model.
#' @param y ordinal response variable
#' @param X covariates that need to be adjusted
#' @param model The cumulative link model formula
#' @param alpha The cumulative intercepts
#' @param beta The coefficients in front of covariates
#' @param ndraw The number of simulation for generate residuals
#' @param dist The argument to specify distribution
#'
#' @return A vector contains surrogate residuals for the binary response y
#' @export
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

  # S <- matrix(rtnorm(n*ndraw, Lhat, 1, lower = lower, upper = upper), n, ndraw)
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

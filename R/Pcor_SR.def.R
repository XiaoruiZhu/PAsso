#' Pcor_SR.def
#'
#' @param var a string vector that specifies response variables.
#' @param confounder a string vector that specifies confounders.
#' @param links a string vector contains default link function for fitting models on
#' the input response variables.
#' @param asso_type a character string indicating whether the partial association coefficient
#' \code{"PartialAsso"} (default) or the marginal correlation \code{"Marginal"} is to be
#' computed. They can be abbreviated.
#' @param cor_method correlation coefficient method. Three choices \code{c("kendall", "pearson", "wolfsigma")},
#' and the default is \code{"kendall"}
#' @param surro_method Character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"jitter"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{surrogate approach;}
#'   \item{\code{jitter}}{jittering approach;}
#'   \item{\code{Sign}}{Sign-based residuals;}
#'   \item{\code{General}}{Generalized residuals (Franses and Paap 2001);}
#'   \item{\code{Deviance}}{Deviance residuals (-2*loglik).}
#' }
#' @param rep_num the repeat times for getting the average as estimator of partial
#' association pho. It is the \code{"nsim"} in \code{"resids()"} function
#' @param ... Additional optional arguments. (say \code{"method"} and \code{"jitter.scale"}
#' in \code{"resids"} function. Just keep them as default if not specified. Currently ignored.)
#' @param data the data.frame including responses and confounders.
#'
#' @import Rfast
#'
#' @return A list of class \code{c("PartialCor")} containing
#' the partial correlation matrix and multiple repeats if \code{rep_num} > 1.
#' This object has "arguments" attribute to save c(asso_type, cor_method, surro_method).
#' The list contains:
#' \describe{MatCor=MatCor, rep_MatCor=rep_MatCor, rep_SRs=rep_SRs, models
#'   \item{\code{MatCor}}{The estimated correlation matrix(average of \code{rep_MatCor}) of partial association controling confounders;}
#'   \item{\code{rep_MatCor}}{The replications of estimated correlation matrix;}
#'   \item{\code{rep_SRs}}{The replications of surrogate residuals if partial association is applied;}
#'   \item{\code{models}}{The list stores all input models.}
#' }
#' @export
#'
#' @examples
#' # import data
#' data(nes96)
#' Pcortest_1 <- Pcor_SR.def(var = c("vote.num", "PID"),
#' confounder = c("income.num", "age", "edu.year"),
#' asso_type = "PartialAsso", cor_method = "kendall",
#' rep_num = 30, data = nes96, links=c("probit", "probit"))
#'
Pcor_SR.def <- function(var, confounder, data,
                        links=c("probit", "probit"),
                        asso_type = c("PartialAsso", "Marginal"),
                        cor_method = c("kendall", "pearson", "wolfsigma"),
                        surro_method = c("latent", "jitter", "Sign", "General", "Deviance"),
                        rep_num=100, ...){

  # var <- c("vote.num", "PID")
  # confounder <- c("income.num", "age", "edu.year")
  # asso_type = "PartialAsso"; cor_method = "kendall";
  # rep_num = 30;
  # data = nes96; links=c("probit", "probit")

  n_var <- length(var) # Only work for two response variables
  if (n_var>2) { stop("More than two responses, please calculate them pair by pair!") }
  n_confounder <- length(confounder)

  # Pull out ingredients for cooking --------------------------------------------
  # mods_y <- as.numeric(sapply(models, function(mod) model.frame(mod)[,1]))
  mods_y <- as.matrix(data[,var])

  formulaAll <-
    sapply(var,
           function(X) paste(X, paste(confounder, collapse = "+"), sep = "~"))

  for (i in var) { # Change all response variables to factor type to avoid error!
    data[,i] <- as.factor(data[,i])
  }

  # Use different functions for if response has differen levels!
  if (length(unique(data[,var[1]])) > 2) { # If response has more than 2 levels, use "polr", otherwise "glm".
    fit_1 <- do.call("polr", list(formula = as.formula(formulaAll[1]),
                                 method = links[1], data = quote(data)))
    # FIXED: formula now is shown as what it is (by "do.call" and "quote")!
  } else {
    fit_1 <- do.call("glm",
                     list(formula = as.formula(formulaAll[1]),
                          family = quote(binomial(link = links[1])),
                          data = quote(data)))
  }
  # summary(fit_1)

  if (length(unique(data[,var[2]])) > 2) { # If response has more than 2 levels, use "polr", otherwise "glm".
    fit_2 <- do.call("polr", list(formula = as.formula(formulaAll[2]),
                                  method = links[2], data = quote(data)))
  } else {
    fit_2 <- do.call("glm",
                     list(formula = as.formula(formulaAll[2]),
                          family = quote(binomial(link = links[2])),
                          data = quote(data)))
  }
  # summary(fit_2)


  # Match arguments -------------------------------------------------
  if (missing(asso_type)) { asso_type <- "PartialAsso" }
  asso_type <- match.arg(asso_type)

  if (missing(cor_method)) { cor_method <- "kendall" }
  cor_method <- match.arg(cor_method)

  if (missing(surro_method)) { surro_method <- "latent" }
  surro_method <- match.arg(surro_method)

  # Initialize the cor_func for different methods
  cor_func <- switch(cor_method,
                     kendall = function(X) cor(X, method = "kendall"),
                     pearson = function(X) cor(X),
                     wolfsigma = function(X)
                       t(copBasic::wolfCOP(para = data.frame(X), as.sample = TRUE)))

  # Combine Models together & Display Warnings --------------------------------------------------------
  models <- list(fit_1, fit_2)
  mods_n <- sapply(models, nobs)
  mods_ys_names <- var # sapply(models, function(mod) colnames(model.frame(mod))[1])

  if (!all(mods_n==mods_n[1])) { # if data are not same, stop!
    stop("Stop due to different data in the models!")
  }

  # mods_xs <- sapply(models, function(mod) all.vars(mod$formula))
  # Above line does not work universally, need to be adjusted, but if input confounder, it's ok.
  # if (all(mods_xs==mods_xs[1])){ # [FEATRUE]: This need to check in the future!
  # IF the var and confounder are inputs, it should be fine.
  #   stop("Stop due to different controlling covariates in these models!")
  # }


  # Main Body ---------------------------------------------------------------
  if (asso_type=="Marginal") { # marginal association (Default Kendall's tau)

    MatCor <- rep_MatCor <- cor(mods_y, method = cor_method)

  } else { # Partial Association based on surrogate residuals
    # Only simulate once to get correlation
    # surro_method <- "latent"
    MatCor <- matrix(1, nrow = n_var, ncol = n_var)

    if (rep_num==1) { # Use just one replication of SR to calcualte partial correlation!
      rep_SRs <- sapply(models, function(mod)
        resids(mod, method = surro_method, nsim = 1))
      colnames(rep_SRs) <- mods_ys_names

      if (cor_method == "wolfsigma") { # wolfsigma only return one value!
        MatCor[upper.tri(MatCor)] <- MatCor[lower.tri(MatCor)] <- cor_func(rep_SRs)
      } else {
        MatCor <- cor_func(rep_SRs)
      }
      rep_MatCor <- MatCor # Save Replication again!

    } else { # Repeat rep_num times to get SRs and calcualte average partial correlation!
      rep_SRs <- sapply(models, function(mod)
        attr(resids(mod, method = surro_method, nsim = rep_num), "boot_reps"),
        simplify = "array") # Use array to save surrogate residuals, rep_SRs[,i,] is i-th simulation
      # dim(rep_SRs)
      dimnames(rep_SRs)[[3]] <- mods_ys_names

      rep_MatCor_temp <- apply(X = rep_SRs, MARGIN = 2, cor_func) # For each copy, calculate correlation
      MatCor <- matrix(apply(X = rep_MatCor_temp, MARGIN = 1, mean), nrow = n_var, ncol = n_var)

      # ABOVE: Take mean as the estimate of partial correlation matrix!

      rep_MatCor <- array(rep_MatCor_temp, dim = c(n_var, n_var, rep_num),
                           dimnames = list(mods_ys_names, mods_ys_names))
    }
  }

  # Give names for correlation matrix and "rep_SRs" array
  colnames(MatCor) <- rownames(MatCor) <-
    # colnames(sd_MatCor) <- rownames(sd_MatCor) <-
    mods_ys_names

  PartialCor <- list(MatCor=MatCor, rep_MatCor=rep_MatCor, rep_SRs=rep_SRs,
                     models=models, data=data, mods_n=mods_n, mods_y=mods_y, cor_func=cor_func)
  attr(PartialCor, "arguments") <- c(asso_type, cor_method, surro_method)
  attr(PartialCor, "var") <- var
  attr(PartialCor, "confounder") <- confounder

  class(PartialCor) <- c("PartialCor")

  return(PartialCor)

  } ## end of function

print.PartialCor <- function(object) {
  cat("The partial correlation coefficient matrix: \n")
  print(object$MatCor)

  cat("--------------------------------------------\n")
  cat("\nThe fitted model of first response variable: \n")
  print(object$models[[1]])

  cat("--------------------------------------------\n")
  cat("\nThe fitted model of second response variable: \n")
  print(object$models[[2]])
}


#' The Correlation and partial correlation coefficients between
#' ordinal responses after adjusting for a set of covariates.
#'
#' This function is mainly designed for conducting the partial association analysis
#' bewteen two or more ordinal response variables after adjusting for a set of
#' covariates. It can also conduct marginal correlation analyses without adjusting
#' covariates. It generates an object that has correlation coefficients matrix,
#' and some attributes ("arguments" saves c(association, method, resids.method)).
#' The attribute "responses" contains the names of response variables.
#'
#' @param responses A string vector that specifies response variables. It requires to be equal
#' or greater than two variables in the data frame.
#' @param adjustments A string vector specifies covariates/confounders that need to
#' be adjusted.
#' @param data the data.frame including responses and adjustments.
#' @param association A character string indicating whether the partial association
#' coefficient \code{"partial"} (default) or the marginal correlation \code{"marginal"}
#' is to be computed. They can be abbreviated.
#' @param models A string vector contains default link function of fitting models with
#' respect to each response variable.
#' @param method A string argument to specify correlation coefficient method.
#' Three choices \code{c("kendall", "pearson", "wolfsigma")}. The default is
#' \code{"kendall"}
#' @param resids.method A character string specifying which method to use to generate the
#' residuals. Current options are \code{"latent"} and
#' \code{"jitter"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{surrogate approach;}
#'   \item{\code{jitter}}{jittering approach;}
#'   \item{\code{Sign}}{Sign-based residuals;}
#'   \item{\code{General}}{Generalized residuals (Franses and Paap 2001);}
#'   \item{\code{Deviance}}{Deviance residuals (-2*loglik).}
#' }
#' @param fitted.models A list contains all the models (S3 objects) you want to
#' assess for the partial association between ordinal responses after adjusting
#' for a set of covariates covariates. All of these models should be applied to the
#' same dataset, having same covariates, same sample size etc. The models in this
#'list can be an object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param rep_num A number to specify repeat times of simulation of surrogate residuls
#' such that the partial correlation coefficients are calculated repeatly. The final
#' correlation coefficents the average of all partial correlation coefficients.
#' It is the \code{"nsim"} argument in \code{"resids()"} function.
#' @param ... Additional optional arguments. (say \code{"method"} and \code{"jitter.scale"}
#' in \code{"resids"} function. Just keep them as default if not specified. Currently ignored.)
#'
#' @import matrixStats MASS
#'
#' @return An object of class \code{"PAsso"} is a list containing at least the following
#' components if \code{"partial"} is the association argument. It contains the partial
#' correlation matrix and multiple repeats if \code{rep_num} > 1. This object has "arguments"
#' attributes saved as c(association, method, resids.method).
#' The list contains:
#' \describe{
#'   \item{\code{corr}}{The estimated correlation matrix(average of \code{rep_MatCorr}) of partial association after
#'   adjusting confounders;}
#'   \item{\code{rep_corr}}{The replications of estimated correlation matrix;}
#'   \item{\code{rep_SRs}}{The replications of surrogate residuals if partial association is applied;}
#'   \item{\code{fitted.models}}{The list stores all fitted.models;}
#'   \item{\code{data}}{The data.frame of original dataset;}
#'   \item{\code{mods_n}}{The sample size of each fitted model;}
#'   \item{\code{cor_func}}{The correlation function after assign different method.}
#' }
#' @references
#' Dungang Liu, Shaobo Li and Yan Yu. Assessing Partial Association Between Ordinal
#' Variables: A General Framework
#'
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted). URL
#' http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20
#'
#' @export
#'
#' @examples
#' #
#' # Adjacent Categories Regression Model Example to compare different residuals
#' #
#'
#' data("df_ParA")
#' summary(df_ParA$data)
#' fit_clm1 <- VGAM::vglm(Y1 ~ X, family =
#'                       VGAM::cumulative(link = "logit",reverse=TRUE,parallel = TRUE),
#'                       data = df_ParA$data)
#' fit_clm2 <- VGAM::vglm(Y2 ~ X, family =
#'                        VGAM::cumulative(link = "logit",reverse=TRUE,parallel = TRUE),
#'                        data = df_ParA$data)
#' SR1 <- resids(fit_clm1, method = "latent",boot_id = NULL)
#' SR2 <- resids(fit_clm2, method = "latent", boot_id = NULL)
#'
#' ## obtain SBC residuals (Li and Shepherd 2012 JASA/Biometrika)
#' PR1 <- resids(fit_clm1, method = "Sign", boot_id = NULL)
#' PR2 <- resids(fit_clm2, method = "Sign", boot_id = NULL)
#'
#' ## obtain generalized residuals (Franses and Paap 2001 book)
#' GR1 <- resids(fit_clm1, method = "General", boot_id = NULL)
#' GR2 <- resids(fit_clm2, method = "General", boot_id = NULL)
#'
#' ## obtain deviance residuals
#' DR1 <- resids(fit_clm1, method = "Deviance", boot_id = NULL)
#' DR2 <- resids(fit_clm2, method = "Deviance", boot_id = NULL)
#'
#' ## visualize residual vs. residual
#' par(mfrow=c(2,2))
#' par(mar=c(4, 4.8, 2.5, 1.5))
#'
#' plot(PR1, PR2, pch=".", main = "Sign-based Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")))
#' plot(GR1, GR2, pch=".", main = "Generalized Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")), xlim = c(-4,4), ylim=c(-4,4))
#' plot(DR1, DR2, pch='.', main = "Deviance Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")))
#' plot(SR1, SR2, pch=".", main = "Surrogate Residuals", xaxt="n", yaxt="n",
#'      xlab = expression(R[1]), ylab = expression(R[2]),
#'      xlim = c(-1/2,1/2), ylim=c(-1/2,1/2))
#' axis(1, at=seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))
#' axis(2, at=seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))
#'
#' # Import nes96 data in "parasol"
#' data(nes96)
#' # Parial association:
#' PAsso_1 <- corr(responses = c("vote.num", "PID"),
#'                 adjustments = c("income.num", "age", "edu.year"),
#'                 data = nes96,
#'                 association = c("partial"),
#'                 models = c("probit", "probit"),
#'                 method = c("kendall"),
#'                 resids.method = "latent",
#'                 fitted.models = NULL, rep_num = 100)
#'
#' # Marginal association:
#' MAsso_1 <- corr(responses = c("vote.num", "PID"),
#'                 adjustments = c("income.num", "age", "edu.year"),
#'                 data = nes96,
#'                 association = c("marginal"))
#'
#' # Compare marginal correlation with partial correlation.
#' PAsso_1
#' MAsso_1
#'
#'
corr <- function(responses, adjustments, data,
                 association = c("partial", "marginal"),
                 models = c("probit", "probit"),
                 method = c("kendall", "pearson", "wolfsigma"),
                 resids.method = c("latent", "jitter", "Sign", "General", "Deviance"),
                 fitted.models = NULL,
                 rep_num = 100, ...){

  # TEST HEADER:
  # responses <- c("vote.num", "PID")
  # adjustments <- c("income.num", "age", "edu.year")
  # association = "partial"; method = "kendall"; resids.method = "latent"
  # rep_num = 30; data = nes96;
  # models = c("probit", "probit")

  # Match arguments -------------------------------------------------
  call <- match.call()
  if (missing(association)) { association <- "partial" }
  association <- match.arg(association)

  if (missing(method)) { method <- "kendall" }
  method <- match.arg(method)

  if (missing(resids.method)) { resids.method <- "latent" }
  resids.method <- match.arg(resids.method)

  # Initialize the cor_func for different methods
  cor_func <- switch(method,
                     kendall = function(X) cor(X, method = "kendall"),
                     pearson = function(X) cor(X),
                     wolfsigma = function(X)
                       t(copBasic::wolfCOP(para = data.frame(X), as.sample = TRUE)))

  # Pull out ingredients for cooking --------------------------------------------
  if (!missing(fitted.models) & is.list(fitted.models)) { # If fitted.models is imported as a list, then done!
    responses <- sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])
    n_responses <- length(responses)

    mods_ys <- sapply(fitted.models, function(mod) model.frame(mod)[,1])
    colnames(mods_ys) <- c(responses)
    conf_temp <- lapply(fitted.models, function(mod) colnames(model.frame(mod)[1,-1]))

    if (length(unique.default(conf_temp)) != 1L) { # If data are not same, stop!
      stop("The imported fitted.models have different confounders!")
    } else { # If confounders are same, use them!
      adjustments <- conf_temp[[1]]; n_adjustments <- length(adjustments)
    }

    data <- data.frame(mods_ys, model.frame(fitted.models[[1]])[,-1]) # Make sure data is data.frame to avoid issue!
    colnames(data) <- c(responses, adjustments)

  } else { # Start to deal with "responses" and "adjustments"
    if ((length(responses) != length(models)) | ((!missing(responses)) & missing(models))) {
      # models_Q <- utils::menu(choices = c("Yes", "No"),
      #                         title="Not all models of the responses are specified. Will you sse same model for all responses?")
      # if (models_Q) {
      # message("Use same probit model for all models!")
      models <- rep("probit", length(responses))

    }

    n_responses <- length(responses)
    n_adjustments <- length(adjustments)

    # mods_ys <- as.numeric(sapply(models, function(mod) model.frame(mod)[,1]))
    mods_ys <- as.matrix(data[,responses])

    formulaAll <-
      sapply(responses,
             function(X) paste(X, paste(adjustments, collapse = "+"), sep = "~"))

    for (i in responses) { # Change all response variables to factor type to avoid error!
      data[,i] <- as.factor(data[,i])
    }

    # Combine Models together -------------------------------
    # Use two different functions (polr, glm) if response has differen levels!
    # NEEDED FEATURE: Test for other packages.
    fitted.models <- list()
    for (i in 1:n_responses) {
      if (length(unique(data[,responses[i]])) > 2) { # If response has more than 2 levels, use "polr", otherwise "glm".
        fitted_temp <- do.call("polr", list(formula = as.formula(formulaAll[i]),
                                            method = models[i], data = quote(data)))
        # assign(x = paste("fitted", responses[i], sep = "_"),
        # fitted_temp)
        fitted.models[[i]] <- fitted_temp
        # FIXED: formula now is shown as what it is (by "do.call" and "quote")!
      } else {
        fitted_temp <- do.call("glm",
                               list(formula = as.formula(formulaAll[i]),
                                    family = quote(binomial(link = models[i])),
                                    data = quote(data)))
        # assign(x = paste("fitted", responses[i], sep = "_"),
        # fitted_temp)
        fitted.models[[i]] <- fitted_temp
        # summary(fitted_vote.num)
        # identical(fitted_temp, fitted_vote.num)
      }
    }
  }

  mods_n <- sapply(fitted.models, nobs)
  # mods_ys_names <- responses # sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])

  if (length(unique.default(mods_n)) != 1L) { # if sample sizes of all models object are not same, stop!
    stop("Stop due to different data in the models!")
  }

  # Main Body ---------------------------------------------------------------
  if (association=="marginal") { # marginal association (Default Kendall's tau)

    MatCorr <- cor(mods_ys, method = method)
    attr(MatCorr, "arguments") <- c(association, method)
    attr(MatCorr, "responses") <- responses
    class(MatCorr) <- c("MAsso", class(MatCorr))

    return(MatCorr)

  } else { # Partial Association based on surrogate residuals
    # Only simulate once to get correlation
    MatCorr <- matrix(1, nrow = n_responses, ncol = n_responses, dimnames = list(responses, responses))

    if (rep_num==1) { # Use just one replication of SR to calcualte partial correlation!
      rep_SRs <- sapply(fitted.models, function(mod)
        resids(mod, method = resids.method, nsim = 1))
      colnames(rep_SRs) <- responses

      if (method == "wolfsigma") { # wolfsigma only return one value!
        MatCorr[upper.tri(MatCorr)] <- MatCorr[lower.tri(MatCorr)] <- cor_func(rep_SRs)
      } else {
        MatCorr <- cor_func(rep_SRs)
      }
      rep_MatCorr <- MatCorr # Save Replication again!

    } else { # Repeat rep_num times to get SRs and calcualte average partial correlation!
      rep_SRs <- sapply(fitted.models, function(mod)
        attr(resids(mod, method = resids.method, nsim = rep_num), "boot_reps"),
        simplify = "array") # Use array to save surrogate residuals, rep_SRs[,i,] is i-th simulation
      # dim(rep_SRs)
      dimnames(rep_SRs)[[3]] <- responses

      rep_MatCorr_temp <- apply(X = rep_SRs, MARGIN = 2, cor_func) # For each copy, calculate correlation
      MatCorr <- matrix(apply(X = rep_MatCorr_temp, MARGIN = 1, mean),
                       nrow = n_responses, ncol = n_responses, dimnames = list(responses, responses))
      # ABOVE: Take mean as the estimate of partial correlation matrix!

      rep_MatCorr <- array(rep_MatCorr_temp, dim = c(n_responses, n_responses, rep_num),
                           dimnames = list(responses, responses))
    }

    PartialAsso <- list(corr=MatCorr, rep_corr=rep_MatCorr, rep_SRs=rep_SRs,
                       fitted.models=fitted.models, data=data, mods_n=mods_n,
                       cor_func=cor_func)
    attr(PartialAsso, "arguments") <- c(association, method, resids.method)
    attr(PartialAsso, "responses") <- responses
    attr(PartialAsso, "adjustments") <- adjustments

    class(PartialAsso) <- c("PAsso", class(PartialAsso))
    return(PartialAsso)
  }
} ## end of function


#' @rdname print
#' @method print PAsso
#'
#' @export
print.PAsso <- function(x, ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n")
  print(signif(x$corr, ...))
}


#' @rdname summary
#' @method summary PAsso
#'
#' @export
summary.PAsso <- function(object, ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n")
  print(signif(object$corr, ...))

  cat("--------------------------------------------\n")
  cat("\nThe fitted models of the response variables are: \n", sep = "")
  print(object$fitted.models)
}


#' @rdname summary
#' @method summary MAsso
#'
#' @export
summary.MAsso <- function(object, ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n")
  print(signif(object, ...))

  cat("--------------------------------------------\n")
  cat("\nThe fitted models of the response variables are: \n", sep = "")
  print(object$fitted.models)
}



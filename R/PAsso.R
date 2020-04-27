#' Partial association analysis between ordinal responses after adjusting for a set of covariates
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
#' @param uni.model A character string specifying the universal model setting for all
#' responses. Default \code{"probit"} refers to cumulative probit model. \code{"logit"}
#' refers to cumulative logit model.
#' @param models A string vector contains default link functions of fitting models with
#' respect to each response variable. If \code{"uni.model"} is provided, this one will
#' be generated automaticlly.
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
#' list can be an object of class \code{\link[ordinal]{clm}},
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
#' @import MASS
#'
#' @return An object of class \code{"PAsso"} is a list containing at least the following
#' components if \code{"partial"} is the association argument. It contains the partial
#' correlation matrix and multiple repeats if \code{rep_num} > 1. This object has "arguments"
#' attribute saved as c(association, method, resids.method), "responses" attribute, and
#' "adjustments" attribute.
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
#' Liu, Dungang, Li, Shaobo, Yu, Yu, and Moustaki, Irini. Assessing partial association between
#' ordinal variables: quantification, visualization, and hypothesis testing, \emph{Journal of the
#' American Statistical Association}, Revision under review.
#'
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association}.
#' \url{http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20}
#'
#' @importFrom ggplot2  aes_string geom_abline geom_boxplot geom_point
#'
#' @importFrom ggplot2 geom_smooth ggplot ggtitle guides labs xlab ylab
#'
#' @importFrom stats .checkMFClasses lowess median model.frame model.matrix
#'
#' @importFrom stats model.response nobs pbinom pcauchy plogis pnorm ppoints
#'
#' @importFrom stats predict qcauchy qlogis qnorm qqline qqplot qqnorm quantile
#'
#' @importFrom stats qunif runif
#'
#' @export
#'
#' @examples
#' # Import nes2016_pre data in "PAsso"
#' data(nes2016_pre)
#' # Parial association:
#' PAsso_1 <- PAsso(responses = c("Prevote.num", "PID"),
#'                 adjustments = c("income.num", "age", "edu.year"),
#'                 data = nes2016_pre,
#'                 association = c("partial"),
#'                 models = c("probit", "probit"),
#'                 method = c("kendall"))
#'
#' PAsso_1
#'
#' @useDynLib PAsso
PAsso <- function(responses, adjustments, data,
                 association = c("partial", "marginal"),
                 uni.model = c("probit", "logit"),
                 models = NULL,
                 method = c("kendall", "pearson", "wolfsigma"),
                 resids.method = c("latent", "jitter", "Sign", "General", "Deviance"),
                 fitted.models = NULL,
                 rep_num = 30, ...){

  # TEST HEADER:
  # data=boot_data
  # responses = attr(object, "responses")
  # adjustments = attr(object, "adjustments")
  # models = attr(object, "models")
  # association = arguments[1]
  # method = arguments[2]
  # resids.method = arguments[3]

  # responses <- c("vote.num", "PID")
  # adjustments <- c("income.num", "age", "edu.year")
  # association = "partial"; method = "kendall"; resids.method = "latent"
  # rep_num = 30; data = nes2016;
  # models = c("probit", "probit")

  n_resp <- ifelse(missing(responses), length(fitted.models), length(responses))

  # Match arguments -------------------------------------------------
  call <- match.call()

  if (missing(uni.model) & missing(models)) {
    uni.model <- "probit"
    models <- rep(uni.model, n_resp)
    uni.model <- match.arg(uni.model)
  }

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



  # Main Body ---------------------------------------------------------------
  if (association=="marginal") { # marginal association (Default Kendall's tau)
    # Pull out ingredients for cooking marginal association --------------------------------------------
    if (missing(responses) & association=="marginal") {
      stop("For marginal association analysis, responses are required!")
    }

    mods_ys <- apply(X = data[,responses], MARGIN = 2, FUN = as.numeric)
    # Save numeric response for calculating marginal corr

    MatCorr <- cor(mods_ys, method = method)
    attr(MatCorr, "arguments") <- c(association, method)
    attr(MatCorr, "responses") <- responses
    class(MatCorr) <- c("MAsso", class(MatCorr))

    return(MatCorr)

  } else { # Partial Association based on surrogate residuals
    # Pull out ingredients for cooking partial association --------------------------------------------
    if (!missing(fitted.models) & is.list(fitted.models)) { # If fitted.models is imported as a list, then done!

      conf_temp0 <- lapply(fitted.models, function(mod) colnames(model.frame(mod)[1,]))

      responses <- lapply(conf_temp0, `[[`, 1)
      n_responses <- length(responses)

      mods_ys <- sapply(fitted.models, function(mod) as.numeric(model.frame(mod)[,1]))
      colnames(mods_ys) <- c(responses)
      # Save numeric response for calculating marginal corr

      conf_temp <- lapply(conf_temp0, function(mod) mod[-1])


      # responses <- sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])
      # n_responses <- length(responses)
      #
      # mods_ys <- sapply(fitted.models, function(mod) model.frame(mod)[,1])
      # colnames(mods_ys) <- c(responses)
      # conf_temp <- lapply(fitted.models, function(mod) colnames(model.frame(mod)[1,-1]))

      if (length(unique.default(conf_temp)) != 1L) { # If data are not same, stop!
        stop("The imported fitted.models have different confounders!")
      } else { # If confounders are same, use them!
        adjustments <- conf_temp[[1]]; n_adjustments <- length(adjustments)
      }

      data <- data.frame(mods_ys, model.frame(fitted.models[[1]])[,-1]) # Make sure data is data.frame to avoid issue!
      colnames(data) <- c(responses, adjustments)

      # Obtain models(links)
      distName <- c()
      for (i in 1:n_responses) {
        distName[i] <- getDistributionName(fitted.models[[i]])
      }
      models <- sapply(X = distName,
                       FUN = function(x) switch(x,
                                                "logis" = "logit",
                                                "norm" = "probit",
                                                "gumbel" = "loglog",
                                                "Gumbel" = "cloglog",
                                                "cauchy" = "cauchit"))

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
      mods_ys <- apply(X = data[,responses], MARGIN = 2, FUN = as.numeric)
      # Save numeric response for calculating marginal corr

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
                                              Hess = TRUE, # Need this to draw coefficients and t-values as output
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

      # Transfer the data again to obtain numeric responses!!!
      data <- data.frame(mods_ys, model.frame(fitted.models[[1]])[,-1]) # Make sure data is data.frame to avoid issue!
      colnames(data) <- c(responses, adjustments)
    }

    mods_n <- sapply(fitted.models, nobs)
    # mods_ys_names <- responses # sapply(fitted.models, function(mod) colnames(model.frame(mod))[1])

    if (length(unique.default(mods_n)) != 1L) { # if sample sizes of all models object are not same, stop!
      stop("Stop due to different data in the models!")
    }


    # Simulate surrogate residuals for calculating the correlation -----------------------------------
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

    # Add marginal association!
    Marg_corr <- cor(mods_ys, method = method)

    PartialAsso <- list(corr=MatCorr, rep_corr=rep_MatCorr, rep_SRs=rep_SRs,
                       fitted.models=fitted.models, data=data, mods_n=mods_n,
                       cor_func=cor_func, Marg_corr=Marg_corr)
    attr(PartialAsso, "arguments") <- c(association, method, resids.method)
    attr(PartialAsso, "responses") <- responses
    attr(PartialAsso, "adjustments") <- adjustments
    attr(PartialAsso, "models") <- models

    class(PartialAsso) <- c("PAsso", class(PartialAsso))
    return(PartialAsso)
  }
} ## end of function


#' @title Print partial association matrix
#' @rdname print
#' @method print PAsso
#'
#' @export
print.PAsso <- function(x, digits = max(2, getOption("digits")-2), ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n")

  # x$corr[lower.tri(x$corr)] <- NA

  # print.default(format(x$corr, digits = max(2, (digits))),
                # print.gap = 2, na.print = "",
                # quote = FALSE, ...)

  temp <- format(x$corr, digits = max(2, (digits)), ...)
  temp[lower.tri(temp)] <- NA

  print.default(temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)
}

#' @title Summary of partial association analysis
#' @description This function summarizes the partial association analysis by
#' providing partial association matrix, marginal association matrix, and a
#' matrix of models' estimation.
#'
#' @rdname summary
#' @method summary PAsso
#'
#' @export
summary.PAsso <- function(object, digits = max(3L, getOption("digits")-2L), ...) {
  cat("-------------------------------------------- \n")
  cat("The partial correlation coefficient matrix: \n\n")
  # print(signif(object$corr, ...))

  temp <- format(object$corr, digits = max(2, (digits)), ...)
  temp[lower.tri(temp)] <- NA

  print.default(temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)

  cat("-------------------------------------------- \n")
  cat("The marginal correlation coefficient matrix: \n\n")
  Marg_temp <- format(object$Marg_corr, digits = max(2, (digits)), ...)
  Marg_temp[lower.tri(Marg_temp)] <- NA

  print.default(Marg_temp,
                print.gap = 2, na.print = "",
                quote = FALSE, ...)

  cat("\n--------------------------------------------\n")
  cat("--------------------------------------------\n")
  cat("\nThe coefficients of fitted models are: \n\n", sep = "")
  # print(object$fitted.models)

  # object=PAsso_1; digits=max(3L, getOption("digits")-2L)

  responses <- attr(object, "responses")
  adjustments <- attr(object, "adjustments")
  n_resp <- length(responses)
  n_adju <- length(adjustments)
  n_samp <- nobs(object$fitted.models[[1]])

  coefs_table <- matrix(NA, nrow = n_adju*3, ncol = n_resp)
  colnames(coefs_table) <- responses
  rownames_coefs_table <- rep(NA, n_adju*2)
  rownames_coefs_table[seq(1, n_adju*3, by = 3)] <- adjustments
  rownames_coefs_table[seq(2, n_adju*3, by = 3)] <- rep("Std. Error", n_adju)
  rownames_coefs_table[seq(3, n_adju*3, by = 3)] <- rep("---", n_adju)

  rownames(coefs_table) <- rownames_coefs_table

  for (i in 1:n_resp) {
    # i <- 1
    # Obtain results from the summary output
    sumry <- summary(object$fitted.models[[i]])$coefficients

    if (inherits(object$fitted.models[[i]], "polr")) {
      # Obtain coefficients and standard error
      coefs_se <- sumry[1:(n_adju), ]
      # Obtain p-value
      temp_p <- 2*pt(-abs(coefs_se[,3]), df = n_samp-1)
      `Pr` <- format.pval(temp_p, digits = max(1L, min(5L, digits - 1L)),
                  eps = .Machine$double.eps)
    } else if (inherits(object$fitted.models[[i]], "glm")) {
      # Obtain coefficients and standard error
      coefs_se <- sumry[-1,1:3]
      # Obtain p-value
      temp_p <- sumry[-1,4]
      `Pr` <- format.pval(temp_p, digits = max(1L, min(5L, digits - 1L)),
                             eps = .Machine$double.eps)
    }
    Signif <- symnum(temp_p, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))

    coefs_se <- format(round(coefs_se, digits = digits),
                       digits = digits)
    coefs_se <- cbind(coefs_se, `Pr`, Signif)

    coefs_table[seq(1, n_adju*3, by = 3),i] <- paste(coefs_se[,1], coefs_se[,5], sep = "")
    coefs_table[seq(2, n_adju*3, by = 3),i] <- coefs_se[,2]
  }
  print.default(coefs_table,
                print.gap = 2, na.print = "",
                quote = FALSE)

  # Add stars of significance level --------------------------------------------------
  if ((w <- getOption("width")) <
      nchar(sleg <- attr(Signif, "legend"))) {
    sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
  }
  cat("Signif. codes:  ", sleg, sep = "",
      fill = w + 4 + max(nchar(sleg, "bytes") - nchar(sleg)))
}




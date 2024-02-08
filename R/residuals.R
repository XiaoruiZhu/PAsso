################################################################################
# S3 generic functions are developed for the class of "clm", "lrm", "orm",
# "polr"
#################################################################################

#' @title Extract Model Residuals
#'
#' @description A generic function to simulate surrogate residuals for cumulative link
#' regression models using the latent method described in Liu and Zhang (2017).
#'
#' It also support the sign-based residuals (Li and Shepherd, 2010), generalized
#' residuals (Franses and Paap, 2001), and deviance residuals for cumulative link
#' regression models.
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param type The type of residuals which should be returned. The alternatives
#' are: "surrogate" (default), "sign", "general", and "deviance". Can be abbreviated.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals;}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#'
#' @param jitter When the \code{type = "surrogate"}, this argument is a character string
#' specifying which method to use to generate the surrogate response values. Current
#' options are \code{"latent"} and \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{latent approach;}
#'   \item{\code{uniform}}{jittering uniform approach.}
#' }
#'
#' @param jitter.uniform.scale When the \code{jitter = "uniform"}, this is a character
#' string specifying the scale on which to perform the jittering whenever
#' \code{jitter = "uniform"}. Current options are \code{"response"} and \code{"probability"}.
#' Default is \code{"response"}.
#'
#' @param nsim An integer specifying the number of replicates to use.
#' Default is \code{1L} meaning one simulation only of residuals.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @return A numeric vector of class \code{c("numeric", "resids")} containing
#' the simulated surrogate residuals. Additionally, if \code{nsim} > 1,
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{draws}}{A matrix with \code{nsim} columns, one for each
#'   is a replicate of the surrogate residuals. Note, they correspond
#'   to the original ordering of the data;}
#'   \item{\code{draws_id}}{A matrix  with \code{nsim} columns. Each column
#'   contains the observation number each surrogate residuals corresponds to in
#'   \code{draws}. (This is used for plotting purposes.)}
#' }
#'
#' @note
#' Surrogate response values require sampling from a continuous distribution;
#' consequently, the result will be different with every call to
#' \code{surrogate}. The internal functions used for sampling from truncated
#' distributions are based on modified versions of
#' \code{\link[truncdist]{rtrunc}} and \code{\link[truncdist]{qtrunc}}.
#'
#' For \code{"glm"} objects, only the \code{binomial()} family is supported.
#'
#' @references
#' Liu, D., Li, S., Yu, Y., & Moustaki, I. (2020). Assessing partial association between
#' ordinal variables: quantification, visualization, and hypothesis testing. \emph{Journal
#' of the American Statistical Association}, 1-14. \doi{10.1080/01621459.2020.1796394}
#'
#' Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal regression models:
#' A surrogate approach. \emph{Journal of the American Statistical Association}, 113(522), 845-854.
#' \doi{10.1080/01621459.2017.1292915}
#'
#' Li, C., & Shepherd, B. E. (2010). Test of association between two ordinal variables
#' while adjusting for covariates. \emph{Journal of the American Statistical Association},
#' 105(490), 612-620. \doi{10.1198/jasa.2010.tm09386}
#'
#' Franses, P. H., & Paap, R. (2001). \emph{Quantitative models in marketing research}.
#' Cambridge University Press. \doi{10.1017/CBO9780511753794}
#'
#' @importFrom stats residuals
#' @name residuals
#' @method residuals clm
#' @export
#'
#' @examples
#' # Example 1
#' # Load data with binary response
#' data(ANES2016)
#' # Fit glm model with binomial logit model
#' fit.prevote <- glm(PreVote.num ~ age + edu.year + income.num,
#'                    data = ANES2016, family = "binomial")
#'
#' # Simulate surrogate residuals
#' res1 <- residuals(fit.prevote,
#'                   type = "surrogate",
#'                   jitter="latent",
#'                   jitter.uniform.scale="response")
#' attr(res1,"arguments")
#'
residuals.clm <- function(object,
                          type = c("surrogate", "sign", "general", "deviance"),
                          jitter = c("latent", "uniform"),
                          jitter.uniform.scale = c("probability", "response"),
                          nsim = 1L, ...) {

  # Sanity check
  if (!inherits(object, c("clm", "glm", "lrm", "orm", "polr"))) {
    stop(deparse(substitute(object)), " should be of class \"clm\", \"glm\", ",
         "\"lrm\", \"orm\", or \"polr\".")
  }

  # Match arguments
  type <- match.arg(type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Switch different type of residuals, need to extract "latent" and "jitter" from "surrogate"
  if (type == "surrogate") {
    if (jitter == "uniform") {
      # Issue warning for jittering surrogate method
      message("Jittering uniform is an experimental feature, use at your own risk!")
    }
    gene_method <- jitter

  } else { # When type is "sign", "general", or "deviance", extract "method" for generate_residuals()
    gene_method <- type
  }

  # Generate surrogate response values
  # Original way, has issue to deal with S4 vglm.
  r <- generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale)


  # Multiple samples
  if (nsim > 1L) {  # multiple draws
    draws <- draws_id <- matrix(nrow = nobs(object), ncol = nsim)
    for(i in seq_len(nsim)) {
      # draws_id[, i] <- sample(nobs(object), replace = TRUE)
      # BUG FIXED: Above original code is not correct! Replicate to get many draws of residuals!
      draws_id[, i] <- seq_along(getResponseValues(object))
      draws[, i] <-
        generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale,
                           draws_id = draws_id[, i, drop = TRUE])
    }
    attr(r, "draws") <- draws
    attr(r, "draws_id") <- draws_id
  }
  attr(r, "names") <- NULL
  attr(r, "arguments") <- c(type, jitter, jitter.uniform.scale)

  # Return residuals
  class(r) <- c("numeric", "resid")
  r

}


################################################################################
# GLM has its own Generic function for cumulative link models;
# I combine my new features with the recent version of residuals.glm (stats 3.6.3).
# 1. It can generate "surrogate", "sign", "general", and "deviance" residuals.
# 3. It has two approaches to generate surrogate residuals: latent and uniform(jittering).
################################################################################


#' @name residuals
#' @method residuals ord
#'
#' @return A "resid" object with attributes. It contains a vector or a matrix (nsim>1) of
#' residuals for the adjacent categories model.
#' @export
residuals.ord <- function (
    object,
    type = c("surrogate", "sign", "general", "deviance",
             "pearson", "working", "response", "partial"),
    jitter = c("latent", "uniform"),
    jitter.uniform.scale = c("probability", "response"),
    nsim = 1L, ...)
{
  type <- match.arg(type)

  # Sanity check: add our new types of residuals to residuals.glm()
  if (type %in% c("surrogate", "sign", "general")) {
    res <- residuals.clm(object = object,
                         type = type,
                         jitter = jitter,
                         jitter.uniform.scale = jitter.uniform.scale,
                         nsim = nsim,...)
  } else if (type %in% c("deviance", "pearson", "working", "response", "partial")) { # Keep original below: residuals.glm()
    res <- stats::residuals.glm(object = object,
                                type = type,
                                ...)
  }
  return(res)
}

#' @rdname residuals
#' @method residuals lrm
#' @export
residuals.lrm <- residuals.clm



#' @rdname residuals
#' @method residuals orm
#' @export
residuals.orm <- residuals.clm



#' @rdname residuals
#' @method residuals polr
#' @export
residuals.polr <- residuals.clm


#' @rdname residuals
#' @method residuals glm
#' @export
residuals.glm <- residuals.ord

#' p_adj_cate
#'
#' @param Z A numerical vector that inputs the latent variable for generating probabilities of adjacent
#' categories regression model.
#'
#' @return A matrix (n by level of respones plus 1) of probabilities of the adjacent categories model.
#'
#' @keywords internal
p_adj_cate <- function(Z){
  k <- ncol(Z)
  p1_pj <- p1_pj_inv<- Z
  ZZ <- 0
  for(j in 1:k){
    ZZ <- ZZ+Z[,j]
    p1_pj[,j] <- exp(ZZ)
    p1_pj_inv[,j] <- 1/p1_pj[,j]
  }
  p1 <- 1/(rowSums(p1_pj_inv)+1)
  pj <- p1_pj
  for(j in 1:k){
    pj[,j] <- p1/p1_pj[,j]
  }
  cbind(1-rowSums(pj), pj)
}


#' generate_residuals_acat
#'
#' @param y A vector inputs the response variable.
#'
#' @param X A data.frame inputs the covariates.
#' @param alpha A vector provides the estimated intercepts of adjacent categories model. If the response
#' has k levels, there should be k+1 numbers in this alpha argument with the k-1 estimated intercepts.
#' The lower bound and upper bound are "-Inf" and "Inf".
#' @param beta A vector provides the estimated coefficients.
#' @param nsim A number to specify the replication of residuals.
#'
#' @return A vector or a matrix (nsim>1) of residuals for the adjacent categories model.
#'
#' @keywords internal
generate_residuals_acat <- function(y, X, alpha, beta, nsim=1){
  # y = y; X = X; alpha = alphas; beta = betas; nsim=1
  # alpha <- matrix(alpha, nrow=1)
  # beta <- matrix(beta, nrow=1)
  n <- length(y)
  z <- sapply(alpha[2:(length(alpha)-1)], function(a) a + tcrossprod(as.matrix(X), beta))
  p_acat <- p_adj_cate(z)
  F_acat <- t(apply(p_acat, 1, cumsum))
  F_acat <- cbind(0, F_acat)
  if(min(y)==0){
    R <- sapply(1:n, function(k) runif(nsim, F_acat[k,y[k]+1], F_acat[k,y[k]+2]))
  }else{
    R <- sapply(1:n, function(k) runif(nsim, F_acat[k,y[k]], F_acat[k,y[k]+1]))
  }
  R - 1/2
}

#' This is a function to deal with the vglm object in S4.
#'
#' @param object An object of class \code{\link[VGAM]{vglm}}.
#'
#' @param type The type of residuals which should be returned. The alternatives
#' are: "surrogate" (default), "sign", "general", and "deviance". Can be abbreviated.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals;}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#' @param jitter A character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{latent approach;}
#'   \item{\code{uniform}}{jittering uniform approach.}
#' }
#' @param jitter.uniform.scale A character string specifying the scale on which to perform
#' the jittering whenever \code{jitter = "uniform"}. Current options are
#' \code{"response"} and \code{"probability"}. Default is \code{"response"}.
#' @param nsim An integer specifying the number of replicates to use.
#' Default is \code{1L} meaning one simulation only of residuals.
#' @param ... Additional optional arguments.
#'
#' @return A "resid" object with attributes. It contains a vector or a matrix (nsim>1) of
#' residuals for the adjacent categories model.
#'
#' @export
#' @keywords internal
residualsAcat <- function(object,
                          type = c("surrogate", "sign", "general", "deviance"),
                          jitter = c("latent", "uniform"),
                          jitter.uniform.scale = c("probability", "response"),
                          nsim = 1L, ...)
{
  # object <- fit
  # type = "surrogate"
  # jitter = "latent"
  # jitter.uniform.scale = "probability"
  # nsim = 1L

  # Match arguments -------------
  type <- match.arg(type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Switch different type of residuals, need to extract "latent" and "jitter" from "surrogate" -------------
  if (type == "surrogate") {
    if (jitter == "uniform") {
      # Issue warning for jittering surrogate method
      message("Jittering uniform is an experimental feature, use at your own risk!")
    }
    gene_method <- jitter
  } else { # When type is "sign", "general", or "deviance", extract "method" for generate_residuals()
    gene_method <- type
  }

  # Generate surrogate response values -------------
  coefs <- coef(object)
  alphas <- matrix(c(-Inf, coefs[1:(ncat(object)-1)], Inf), nrow = 1)
  betas <- matrix(coefs[-c(1:(ncat(object)-1))], nrow = 1)
  # message(alphas, "  ", betas)
  y <- getResponseValues(object)

  # X <- as.matrix(model.frame(object)[,-1])
  X <- as.matrix(object@x[,-1]) # This just work for the "vglm" adjacent categories model!

  r <- generate_residuals_acat(y = y, X = X, alpha = alphas, beta = betas, nsim=1)

  # Multiple samples -------------
  if (nsim > 1L) {  # multiple draws
    draws_id <- matrix(seq(nobs(object)), nrow = nobs(object), ncol = nsim, byrow = F)

    draws <-
      t(generate_residuals_acat(y = y, X = X, alpha = alphas, beta = betas, nsim = nsim))

    attr(r, "draws") <- draws
    attr(r, "draws_id") <- draws_id
  }
  attr(r, "names") <- NULL
  attr(r, "arguments") <- c(type, jitter, jitter.uniform.scale)

  # Return residuals
  class(r) <- c("numeric", "resid")

  return(r)
}


################################################################################
# S4 generic functions are developed for the class of "vgam", or "vglm".
#################################################################################

#' @title Extract Residuals from models of VGAM
#'
#' @description A internal function to simulate surrogate residuals for models fitted by
#' \code{\link[VGAM]{vglm}} and \code{\link[VGAM]{vglm}}. Now, this one support \code{"vglm"} and
#' \code{"vgam"}, and adjacent categories regression model by \code{"vglm"}. This one may
#' need to update to support more models from VGAM.
#'
#' @keywords internal
#' @export
residualsVGAM <- function(object, type = c("surrogate", "sign", "general", "deviance"),
                          jitter = c("latent", "uniform"),
                          jitter.uniform.scale = c("probability", "response"),
                          nsim = 1L, ...) {
  # Sanity check
  if (!inherits(object, c("vgam", "vglm"))) {
    stop(deparse(substitute(object)), " should be of class \"vgam\", or \"vglm\".")
  }
  # Match arguments
  type <- match.arg(type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Switch different type of residuals, need to extract "latent" and "jitter" from "surrogate"
  if (type == "surrogate") {
    if (jitter == "uniform") {
      # Issue warning for jittering surrogate method
      message("Jittering uniform is an experimental feature, use at your own risk!")
    }
    gene_method <- jitter

  } else { # When type is "sign", "general", or "deviance", extract "method" for generate_residuals()
    gene_method <- type
  }

  if (isS4(object) & inherits(object, "vglm") & (attr(slot(object, "family"), "vfamily")[1]=="acat")) {
    # If object is S4 and adjacent categories regression model, need to use "residualsAcat" instead!
    r <- residualsAcat(object, ...)

    # Below is another way to solve the S4 issue, but redundant now.
    # use_func <- "residualsAcat"
    # r <- do.call(what = use_func,
    #              args = list(object = object,
    #                          type = type, jitter = jitter,
    #                          jitter.uniform.scale = jitter.uniform.scale,
    #                          nsim = nsim, ...)
    # )
  } else {
    r <- generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale, ...)
  }

  # Multiple samples
  if (nsim > 1L) {  # multiple draws
    draws <- draws_id <- matrix(nrow = nobs(object), ncol = nsim)

    # If object is S4, need to use "residualsAcat" instead!
    if (isS4(object) & inherits(object, "vglm") & (attr(slot(object, "family"), "vfamily")[1]=="acat")) {

      for(i in seq_len(nsim)) {
        # draws_id[, i] <- sample(nobs(object), replace = TRUE)
        # BUG FIXED: Above original code is not correct! Use same IDs to replicate and get many draws of residuals!
        draws_id[, i] <- seq_along(getResponseValues(object))

        suppressMessages(
          draws[, i] <- residualsAcat(object = object,
                                      type = type,
                                      jitter = jitter,
                                      jitter.uniform.scale = jitter.uniform.scale,
                                      nsim = nsim,...)
        )

        # Below is another way to solve the S4 issue, but redundant now.
        # use_func <- "residualsAcat"
        #
        # draws_id[, i] <- do.call(what = use_func,
        #              args = list(object = object,
        #                          type = type, jitter = jitter,
        #                          jitter.uniform.scale = jitter.uniform.scale,
        #                          nsim = nsim, ...)
        # )

      }

    } else {

      for(i in seq_len(nsim)) {
        # draws_id[, i] <- sample(nobs(object), replace = TRUE)
        # BUG FIXED: Above original code is not correct! Replicate to get many draws of residuals!
        draws_id[, i] <- seq_along(getResponseValues(object))

        draws[, i] <-
          generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale,
                             draws_id = draws_id[, i, drop = TRUE])
      }

    }

    attr(r, "draws") <- draws
    attr(r, "draws_id") <- draws_id
  }

  attr(r, "names") <- NULL
  attr(r, "arguments") <- c(type, jitter, jitter.uniform.scale)

  # Return residuals
  class(r) <- c("numeric", "resid")
  r

}

#' ref the vglm' adjacent categories regression model to using \code{"residualsAcat"} function.
#' @param object An object of class \code{\link[VGAM]{vglm}}.
#' @param type The type of residuals which should be returned. The alternatives
#' are: "surrogate" (default), "sign", "general", and "deviance". Can be abbreviated.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals;}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#'
#' @param jitter When the \code{type = "surrogate"}, this argument is a character string
#' specifying which method to use to generate the surrogate response values. Current
#' options are \code{"latent"} and \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{latent approach;}
#'   \item{\code{uniform}}{jittering uniform approach.}
#' }
#'
#' @param jitter.uniform.scale When the \code{jitter = "uniform"}, this is a character
#' string specifying the scale on which to perform the jittering whenever
#' \code{jitter = "uniform"}. Current options are \code{"response"} and \code{"probability"}.
#' Default is \code{"response"}.
#'
#' @param ... Additional optional arguments (Include all arguments in \code{"residuals"} function).
#'
#' @rdname residuals
#' @method residuals vglm
#' @export
setMethod("residuals",  "vglm",
          definition = residualsVGAM)

#' ref the vgam' adjacent categories regression model to using \code{"residualsAcat"} function.
#' @param object An object of class \code{\link[VGAM]{vgam}}.
#' @param ... Additional optional arguments (Include all arguments in \code{"residuals"} function).
#'
#' @rdname residuals
#' @method residuals vgam
#' @export
setMethod("residuals",  "vgam",
          definition = residualsVGAM)



#' Extract Residuals from a Partial Association Analysis
#'
#' @param object An object of class \code{PAsso}.
#'
#' @param draw_id A number refers to the i-th draw of residuals.
#'
#' @param ... Additional optional arguments.
#'
#' @return A matrix of class \code{c("matrix", "resids")} containing
#' the simulated surrogate residuals used for the partial association
#' analysis in \code{PAsso}. Additionally, if \code{rep_num} > 1 in \code{PAsso},
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{draws}}{An array contains all draws of residuals.}
#' }
#'
#'
#' @importFrom stats residuals
#' @name residuals
#' @method residuals PAsso
#'
#' @export
#'
#' @examples
#' # Example 2
#' # residuals() function can also work for PAsso object
#' # Load data
#' data("ANES2016")
#' PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
#'                  adjustments = c("income.num", "age", "edu.year"),
#'                  data = ANES2016)
#'
#' # Extract surrogate residuals from the PAsso object
#' res1 <- residuals(PAsso_1)
#'
residuals.PAsso <- function(object, draw_id=1, ...) {
  if ((draw_id>=1) & (draw_id<=dim(object$rep_SRs)[2])) {
    resids_PAsso <- object$rep_SRs[,draw_id,]
    # resids_PAsso <- PAsso_1$rep_SRs[,1,]
  } else {
    stop("The draw_id is out of bound.")
  }
  attr(resids_PAsso, "draws") <- object$rep_SRs
  attr(resids_PAsso, "arguments") <- attr(object, "arguments")
  # Return residuals
  class(resids_PAsso) <- c("resid", class(resids_PAsso))
  return(resids_PAsso)
}

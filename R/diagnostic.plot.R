#' @title Residual-based diagnostic plots
#'
#' @description A set of visualization tools for the diagnostic of the fitted model in
#' the partial association analysis. It can provides a plot matrix including Q-Q plots,
#' residual-vs-fitted plots, residual-vs-covariate plots of all the fitted models.
#' This function also support the direct diagnostic of the cumulative link regression model
#' in the class of \code{\link[ordinal]{clm}}, \code{\link[stats]{glm}}, \code{\link[rms]{lrm}},
#' \code{\link[rms]{orm}}, \code{\link[MASS]{polr}}. Currently, \code{\link[VGAM]{vglm}}
#' is not supported.
#'
#' @param object The object in the support classes (This function is mainly designed
#' for \code{\link[PAsso]{PAsso}}).
#'
#' @param output A character string specifying what type of output to plot. Default is
#' \code{"qq"} which produces a plot matrix with quantile-quantile plots of the residuals.
#' \code{"fitted"} produces a plot matrix between residuals and all corresponding fitted responses.
#' \code{"covariates"} produces a plot matrix between residuals and corresponding covariate.
#'
#' @param ... This function is based on a modified version of \code{"autoplot"} function in
#' \code{"sure"}. Additional optional arguments can be passed onto for drawing various plots.
#'
#' @return A list of \code{"ggplot"} objects. For class "PAsso", the very last element
#' of the list is the combined diagnostic plots.
#'
#' @name diagnostic.plot
#'
#' @export
#'
#' @examples
#' # Import data for partial association analysis
#' data("nes2016")
#'
#' PAsso_3v <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
#'                   adjustments = c("income.num", "age", "edu.year"),
#'                   data = nes2016, uni.model = "probit",
#'                   method = c("kendall"),
#'                   resids.type = "surrogate", jitter = "latent")
#'
#' diag_p1 <- diagnostic.plot(object = PAsso_3v, output = "qq")
#' diag_p2 <- diagnostic.plot(object = PAsso_3v, output = "fitted")
#' diag_p3 <- diagnostic.plot(object = PAsso_3v, output = "covariate")
#'
#' # Simply diagnose a model
#' # Fit cumulative link models
#'
#' fit1 <- ordinal::clm(Prevote.num ~ income.num + age + edu.year, data = nes2016, link = "logit")
#'
#' # diagnostic.plot
#' plot_qq_1 <- diagnostic.plot(object = fit1, output = "qq")
#' plot_fit_1 <- diagnostic.plot(object = fit1, output = "fitted")
#' plot_cov_1 <- diagnostic.plot(object = fit1, output = "covariate")
#'
diagnostic.plot <- function(object, ...) {
  UseMethod("diagnostic.plot")
}

#' @rdname diagnostic.plot
#' @method diagnostic.plot default
#' @export
diagnostic.plot.default <- function(
  object, ...
){
  warn_str <- paste("diagnostic.plot does not know how to handle object of class ",
                    class(object),
                    "and can only be used on classes PAsso, PAsso.test, resid, clm, glm, lrm, orm, polr.")
  warning(paste(strwrap(warn_str), collapse = "\n"))
}

#' @rdname diagnostic.plot
#' @method diagnostic.plot resid
#' @export
diagnostic.plot.resid <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.resid(object=object, ...)
}


#' @rdname diagnostic.plot
#' @method diagnostic.plot PAsso
#' @export
diagnostic.plot.PAsso <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  # object = PAsso_2; output = "covariate"

  # What type of output plot to produce
  output <- match.arg(output, several.ok = FALSE)
  rep_SRs <- object$rep_SRs
  resp_name <- attr(object, "responses")

  n_resp <- length(resp_name)
  nCol <- floor(sqrt(n_resp))
  plot_list <- list()

  # return a matrix-plot including diagnostics of models.
  if (output == "qq") {
    for (i in 1:n_resp) {
      plot_list[[i]] <-
        # autoplot(rep_SRs[,1,i], output = output,
        autoplot(object$fitted.models[[i]], output = output,
                 resp_name = resp_name[i], ...)
    }
    # Save the combined plot
    print(do.call("grid.arrange", c(plot_list, ncol=nCol)))

  } else if (output == "fitted") {

    for (i in 1:n_resp) {
      plot_list[[i]] <-
        autoplot(object$fitted.models[[i]], output = output, resp_name = resp_name[i],
                 alpha = 0.5, ...)
    }
    # Save the combined plot
    print(do.call("grid.arrange", c(plot_list, ncol=nCol)))

  } else {
    adjust_name <- attr(object, "adjustments")
    n_adjust <- length(adjust_name)
    t_lenght <- n_resp*n_adjust
    adjust_id <- rep(1:n_adjust, times=n_resp) # make index for covariate name in the for loop
    resp_id <- rep(1:n_resp, each=n_adjust) # make index for response name in the for loop

    for (i in 1:(t_lenght)) {
      plot_list[[i]] <-
        autoplot(object$fitted.models[[resp_id[i]]], output = "covariate",
                 x = object$data[,adjust_name[adjust_id[i]]],
                 xlab = adjust_name[adjust_id[i]],
                 # resp_name = resp_name[resp_id[i]], ...)
                 resp_name = resp_name[resp_id[i]])

      if (i %% n_adjust != 1) { # First plot of each response, draw ylab, otherwise, no ylab
        plot_list[[i]] <- plot_list[[i]] + ylab("")
      }
      if ((i-1) %/% n_adjust != (n_resp-1)) { # last row of plot(last response), draw xlab, otherwise, no xlab
        plot_list[[i]] <- plot_list[[i]] + xlab("")
      }
    }
    # Save the combined plot
    print(do.call("grid.arrange", c(plot_list, ncol=n_adjust)))
  }

  return(plot_list)
}

#' @rdname diagnostic.plot
#' @method diagnostic.plot clm
#' @export
diagnostic.plot.clm <- autoplot.clm

#' @rdname diagnostic.plot
#' @method diagnostic.plot glm
#' @export
diagnostic.plot.glm <- autoplot.glm


#' @rdname diagnostic.plot
#' @method diagnostic.plot lrm
#' @export
diagnostic.plot.lrm <- autoplot.lrm


#' @rdname diagnostic.plot
#' @method diagnostic.plot orm
#' @export
diagnostic.plot.orm <- autoplot.orm


#' @rdname diagnostic.plot
#' @method diagnostic.plot polr
#' @export
diagnostic.plot.polr <- autoplot.polr

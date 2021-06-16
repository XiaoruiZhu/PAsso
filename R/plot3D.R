#' @title 3-D P-P plot and false color level plot for the inspection of the partial association analysis
#'
#' @description A list of 3-D P-P plots (or false color level plots when \code{type = "contour"})
#' for the inspection of the partial association analysis. Each plot is either 3-D P-P plot or level plot
#' from an empirical copula trained from the surrogate residuals of a pair of responses.
#'
#' @param object The input object should be a "PAsso" class that is generated
#' by "PAsso" or "test".
#' @param y1 A string to specify the first response for the 3D plot.
#' @param y2 A string to specify the second response for the 3D plot. If either one of the
#' y1 or y2 is missing. The \code{plot3D} will draw 3D plots or level plots for all pairs of responeses.
#' @param ... Additional optional arguments.
#'
#' @details All the plots are based on surrogate residuals generated from \code{"residuals"}
#' function in \code{sure}. Graphics are designed based on
#' PAsso and \code{"plotly"}.
#'
#' @return If response y1 or y2 is not specified, a list of \code{"plotly"} objects includes
#' all pairs of responses will be returned (with name "response 1 v.s. response 2" etc.). If
#' responses y1 and y2 are specified, returns a 3D plot as \code{"plotly"} object.
#'
#' @rdname plot3D
#' @importFrom copula C.n pobs
#' @importFrom plotly plotly_build plot_ly add_surface layout
#' @importFrom dplyr %>%
#' @export plot3D
#'
#' @examples
#' # Did not run this to save time
#' # data("ANES2016")
#' # PAsso_3v <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
#' #                   adjustments = c("income.num", "age", "edu.year"),
#' #                   data = ANES2016)
#'
#' # plot3D(PAsso_3v, y1="PID", y2="selfLR")
#'
plot3D <- function(object, y1, y2, ...) {
  UseMethod("plot3D")
}

#' @rdname plot3D
#' @method plot3D default
#' @export
plot3D.default <- function(object, y1, y2, ...){
  warning(paste("plot3D does not know how to handle object of class ",
                class(object),
                "and can only be used on classes PAsso"))
}

#' @param plot_list The list to save plots. Each 3D plot is between one pair of covariates.
#'
#' @param rep_SRs The surrogate responses array saved in the PAsso object.
#' @param m The index of first covariate. Since the names of covariates are empty, need to input index number.
#' @param n The index of second covariate.
#' @param plot_titles The title of the plot.
#' @param type A character string specifying the trace type (e.g. "surface3D", "contour"). "contour" creates a
#' 2D contour plot between u and v.
#'
#' @keywords internal
#' This is an internal function to draw one plotly (3D surface or contour plot).
plot3D_one <- function(plot_list, rep_SRs, m, n, plot_titles, type = c("surface3D", "contour")) {

  type <- match.arg(type)

  resi <- cbind(rep_SRs[,1,m], rep_SRs[,1,n])
  empC <- C.n(pobs(resi), resi)
  empFG1 <- pobs(resi)[,1]*pobs(resi)[,2]

  ## 3-D copula plot
  u <- v <- seq(0, 1, length.out = 100)
  aa <- matrix(0, length(u), length(u))
  for(i in 1:length(u)){
    for(j in 1:length(u)){
      aa[i,j] <- C.n(t(as.matrix(c(u[i], v[j]))), resi)-u[i]*v[j]
    }
  }

  # options(Viewer=NULL)
  # name <- paste("", i_plot, sep = "_")
  plot_list[[plot_titles]] <-
    if (type == "surface3D") {
      plotly_build(plot_ly(x = u, y = v, z = 12*aa) %>%
                     add_surface() %>%
                     layout(scene = list(xaxis= list(title= "u"),
                                         yaxis= list(title= "v"),
                                         zaxis= list(title= "12(C(u,v)-uv)")),
                            title = paste("3-D P-P Plot: ", plot_titles, sep = "")))
    } else if (type == "contour") {
      plotly_build(plot_ly(x = u,
                           y = v,
                           z = 12*aa,
                           type = "contour") %>%
                     layout(xaxis = list(title= "u"),
                            yaxis = list(title= "v")) %>%
                     add_annotations(
                       text = paste("Contour Plot: ", plot_titles, sep = ""),
                       x = 0.5,
                       y = 1,
                       yref = "paper",
                       xref = "paper",
                       xanchor = "middle",
                       yanchor = "top",
                       showarrow = FALSE,
                       font = list(size = 15))
        )


    }

  return(plot_list)
}

#' @param object A PAsso class of object.
#'
#' @param y1 A string to specify the first response for the 3D plot.
#' @param y2 A string to specify the second response for the 3D plot. If either one of the
#' y1 or y2 is missing. The \code{plot3D} will draw 3D plots for all pairs of responses.
#' @param type A character string specifying the trace type (e.g. "surface3D", "contour"). "contour" creates a
#' 2D contour plot between u and v.
#' @param ... Additional optional arguments.
#'
#' @rdname plot3D
#' @method plot3D PAsso
#' @export
plot3D.PAsso <- function(
  object,
  y1, y2,
  type = c("surface3D", "contour"),
  ...
) {
  # object = PAsso_2; y1 = "selfLR"; y2 = "PID"

  if (!inherits(object, "PAsso")) stop("Input object must be 'PAsso' class.")

  type <- match.arg(type)

  rep_SRs <- object$rep_SRs
  resp_name <- attr(object, "responses")
  n_resp <- length(resp_name)

  if (missing(y1) | missing(y2)) { # If no input for the responses pair, draw all.
    all_pairs <- apply(expand.grid(resp_name, resp_name), 1, paste, collapse=" v.s. ")
    mat_pairs <- matrix(all_pairs, nrow = n_resp, byrow = F)
    plot_titles <- mat_pairs[upper.tri(mat_pairs)]

    plot_list <- list()
    m <- 1; n <- 2

    for (i_plot in 1:(n_resp*(n_resp-1)/2)) {
      plot_list <- plot3D_one(plot_list = plot_list, rep_SRs = rep_SRs,
                              m = m, n = n,
                              plot_titles = plot_titles[i_plot],
                              type = type)
      # Update response ---------------------------------------------------------
      if (n == n_resp) { m <- m + 1; n <- m + 1 } else n <- n + 1
    }
  } else {
    m <- which(y1 == resp_name)
    n <- which(y2 == resp_name)

    i_plot <- plot_titles <- apply(expand.grid(y1, y2), 1, paste, collapse=" v.s. ")
    plot_list <- list()

    plot_list <- plot3D_one(plot_list = plot_list, rep_SRs = rep_SRs,
                            m = m, n = n,
                            plot_titles = plot_titles,
                            type = type)
    plot_list <- plot_list[[1]]
  }

  return(plot_list)
}

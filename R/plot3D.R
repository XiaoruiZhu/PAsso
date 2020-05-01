#' @title 3-D P-P plot for the inspection of the partial association analysis
#'
#' @description A list of 3-D P-P plots for the inspection of the partial
#' association analysis. Each plot is the 3-D P-P plot from an empirical
#' copula trained from the surrogate residuals of a pair of responses.
#'
#' @param object The input object should be a "PAsso" class that is generated
#' by "Passo" or "test".
#'
#' @details All the plots are based on surrogate residuals generated from \code{"resides"}
#' function in \code{\link{sure}}. Graphics are designed based on
#' \code{\link{PAsso}} and \code{"plotly"}.
#'
#' @return A list of \code{"plotly"} objects with name "plot_1", "plot_2" etc.
#'
#' @rdname plot3D
#' @export plot3D
#'
#' @examples
#' data("nes2016")
#' PAsso_3v <- PAsso(responses = c("vote.num", "PID", "selfLR"),
#'                  adjustments = c("income.num", "age", "edu.year"),
#'                  data = nes2016)
#'
#' plot_all <- plot3D(PAsso_3v)
#' plot_all$plot_1
#'
plot3D <- function(object) {
  UseMethod("plot3D")
}

#' @rdname plot3D
#' @method plot3D default
#' @export
plot3D.default <- function(object, ...){
  warning(paste("plot3D does not know how to handle object of class ",
                class(object),
                "and can only be used on classes PAsso"))
}


#' @rdname plot3D
#' @method plot3D PAsso
#' @export
plot3D.PAsso <- function(
  object, ...
) {
  # object <- PAsso_2

  if (!inherits(object, "PAsso")) stop("Input object must be 'PAsso' class.")

  rep_SRs <- object$rep_SRs
  resp_name <- attr(object, "responses")
  n_resp <- length(resp_name)

  all_pairs <- apply(expand.grid(resp_name, resp_name), 1, paste, collapse=" v.s. ")
  mat_pairs <- matrix(all_pairs, nrow = n_resp, byrow = F)
  plot_titles <- mat_pairs[upper.tri(mat_pairs)]

  plot_list <- list()
  m <- 1; n <- 2

  for (i_plot in 1:(n_resp*(n_resp-1)/2)) {
    resi <- cbind(rep_SRs[,1,m], rep_SRs[,1,n])
    empC <- copula::C.n(pobs(resi), resi)
    empFG1 <- copula::pobs(resi)[,1]*copula::pobs(resi)[,2]

    ## 3-D copula plot
    v1 <- v2 <- seq(0, 1, length.out = 100)
    aa <- matrix(0, length(v1), length(v1))
    for(i in 1:length(v1)){
      for(j in 1:length(v1)){
        aa[i,j] <- copula::C.n(t(as.matrix(c(v1[i], v2[j]))), resi)-v1[i]*v2[j]
      }
    }

    # options(Viewer=NULL)
    name <- paste("plot", i_plot, sep = "_")
    plot_list[[name]] <-
      plotly::plotly_build(plotly::plot_ly(x = v1, y = v2, z = 12*aa) %>%
      plotly::add_surface() %>%
      plotly::layout(scene = list(xaxis= list(title= "u"),
                                 yaxis= list(title= "v"),
                                 zaxis= list(title= "12(C(u,v)-uv)")),
                     title = paste("3-D P-P Plot: ", plot_titles[i_plot], sep = "")))


    # Update response ---------------------------------------------------------
    if (n == n_resp) { m <- m + 1; n <- m + 1 } else n <- n + 1

  }
  return(plot_list)
}

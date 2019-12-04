#' Partial regression plot matrix
#'
#' A pairwise plot matrix reveals the partial association between ordinal variables.
#' All the plots are based on surrogate residuals generated from \code{"resides"} function.
#' Graphics are designed based on \code{\link[ggplot2]{ggplot2}} and
#' \code{\link[GGally]{GGally}}.
#'
#' @param resid_Mat An matrix contains residuals for all ordinal responses.
#'
#'
#' @param color spedify the points color
#' @param ... Additional optional arguments to be passed onto
#' \code{\link[sure]{resids}}.
#' @return A \code{"GGally"} object.
#'
#' @rdname ggpairs.resid
#'
#' @export
#'
#' @examples
#' # See ?resids for an example
#' ?resids
ggpairs.resid <- function(
  resid_Mat,
  color="blue", ...
) {
 GGally::ggpairs(resid_Mat,
         upper = list(
           continuous = GGally::wrap("smooth_loess", colour=color)),
   # diag = list(continuous = wrap("barDiag", colour = color)),
   lower = list(continuous = GGally::wrap("cor"))
 )
}


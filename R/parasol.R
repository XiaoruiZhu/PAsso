#' parasol: An R package for assessing the Partial ASSociation between
#' Ordinal variables. It includes quantification, visualization, and hypothesis testing.
#'
#' The \code{parasol}: a package for assessing the PARtial ASsociation between
#' OrdinaL variables for fitted ordinal and general (e.g., binary)
#' regression models of class
#' \code{\link[ordinal]{clm}}, \code{\link[stats]{glm}},
#' \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' The development version can be found on GitHub:
#' \url{https://github.com/XiaoruiZhu/parasol}. As of right now, \code{parasol} exports the
#' following functions:
#' \itemize{
#'   \item{\code{Pcor_SR2}} - calculate partial association based on
#'   surrogate residuals (or others in the future);
#'   \item{\code{resids}} - construct (surrogate-based) residuals;
#'   \item{\code{surrogate}} - construct sorrogate latent variables;
#'   \item{\code{autoplot}} - plot diagnostics using
#'   \code{\link[ggplot2]{ggplot2}}-based graphics;
#'   \item{\code{ggpairs.resid}} - draw partial regression plot matrix;
#'   \item{\code{gof}} - simulate p-values from a goodness-of-fit test.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
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
#' @docType package
#'
#' @name parasol
#' @useDynLib parasol
#' @importFrom Rcpp sourceCpp
#'
NULL

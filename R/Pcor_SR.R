#' Pcor_SR
#' function to obtain \eqn{hat{phi}}
#'
#' @param y1 ordinal response 1
#' @param y2 ordinal response 2
#' @param x covariates
#' @param link1 link function of ordinal response 1
#' @param link2 link function of ordinal response 2
#' @param model1 The model 1 between ordinal response 1 and covariates
#' @param model2 The model 2 between ordinal response 1 and covariates
#' @param alpha1 Alpha of model 1
#' @param beta1 Betas of model 1
#' @param alpha2 Alpha of model 2
#' @param beta2 Betas of model 2
#' @param n.avg The number of simulation for generate residuals
#' @param method The method to calculate correlation
#'
#' @return A list contains ...
#' \itemize{
##'  \item{"parameter 1"}{Stuff}
##'  \item{"parameter 2"}{Stuff}
##' }
#' @export
Pcor_SR<- function(y1, y2, x, model1, model2,
                   alpha1=NULL, beta1=NULL, alpha2=NULL, beta2=NULL,
                   link1="probit", link2="probit", n.avg=100, method="kendall"){
  if(link1=="probit") distribution1="norm"
  if(link1=="logit" | link1=="logistic"){link1="logistic"; distribution1="logis"}
  if(link1=="cloglog") distribution1="Gumbel"
  if(link1=="loglog") distribution1="gumbel"
  if(link1=="cauchit") distribution1="cauchy"

  if(link2=="probit") distribution2="norm"
  if(link2=="logit" | link2=="logistic"){link2="logistic"; distribution2="logis"}
  if(link2=="cloglog") distribution2="Gumbel"
  if(link2=="loglog") distribution2="gumbel"
  if(link2=="cauchit") distribution2="cauchy"

  n<- length(y1)
  dat<- data.frame(y1, y2, x)

  if(!is.null(alpha1) & !is.null(alpha2) & !is.null(beta1) & !is.null(beta2) & missing(model1) & missing(model2)){
    sr1<- SR(y1, x, alpha=alpha1, beta=beta1, ndraw = n.avg, dist = distribution1)
    sr2<- SR(y2, x, alpha=alpha2, beta=beta2, ndraw = n.avg, dist = distribution2)
    fit1<- fit2<- NULL
  }else{
    if(is.null(c(alpha1, alpha2)) & is.null(c(beta1, beta2)) & !missing(model1) & !missing(model2)){
      sr1<- SR(y1, x, model=model1, ndraw = n.avg, dist = distribution1)
      sr2<- SR(y2, x, model=model2, ndraw = n.avg, dist = distribution2)
    }else{
      stop("Either models or alpha's and beta's should be supplied")
    }
  }

  # calculate \hat{\phi}
  if(method=="kendall"){
    cor.est<- sapply(1:n.avg, function(yy) cor(sr1[,yy], sr2[,yy], method = "kendall"))
  }
  if(method=="pearson"){
    cor.est<- sapply(1:n.avg, function(yy) cor(sr1[,yy], sr2[,yy]))
  }
  if(method=="wolfsigma"){
    cor.est<- t(sapply(1:n.avg, function(yy) wolfCOP(para = data.frame(cbind(sr1[,yy], sr2[,yy])), as.sample = TRUE)))
  }
  return(list(phi=mean(cor.est), phi_all=cor.est, dist1=distribution1, dist2=distribution2, sr1=sr1, sr2=sr2))

}## end of function

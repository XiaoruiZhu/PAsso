#' Pcor_SR_boot function to perform bootstrap inference
#'
#' @param y1 ordinal response 1
#' @param y2 ordinal response 2
#' @param x covariates
#' @param B bootstrap times
#' @param link1 link function of ordinal response 1
#' @param link2 link function of ordinal response 2
#' @param n.avg The number of simulation for generate residuals
#' @param method The method to calculate correlation
#' @param H0 The mu when conduct hypothesis testing
#'
#' @import MASS
#' @return A list contains ...
#' @export
#'
Pcor_SR_boot<- function(y1, y2, x, B=2000, link1="probit", link2="probit", n.avg=100, method="kendall", H0=0){

  pcor.est<- rep(NA, B)
  if(length(unique(y1))>2 & length(unique(y2))>2){
    for(i in 1:B){
      tryCatch({
        index<- sample(length(y1), replace=T)
        y1.b<- y1[index]
        y2.b<- y2[index]
        if(is.vector(x)) x.b<- x[index] else x.b<- x[index,]
        fit1<- polr(as.factor(y1.b)~ x.b, method=link1)
        fit2<- polr(as.factor(y2.b)~ x.b, method=link2)
        pcor.est[i]<- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1=link1, link2=link2, n.avg=n.avg, method=method)$phi
      }, error=function(e){})
    }
  }else{
    if(length(unique(y1))==2 & length(unique(y2))>2){
      for(i in 1:B){
        tryCatch({
          index<- sample(length(y1), replace=T)
          y1.b<- y1[index]
          y2.b<- y2[index]
          if(is.vector(x)) x.b<- x[index] else x.b<- x[index,]
          fit1<- glm(y1.b~ x.b, family = binomial(link =link1))
          fit2<- polr(as.factor(y2.b)~ x.b, method=link2)
          pcor.est[i]<- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1=link1, link2=link2, n.avg=n.avg, method=method)$phi
        }, error=function(e){})
      }
    }else{
      for(i in 1:B){
        tryCatch({
          index<- sample(length(y1), replace=T)
          y1.b<- y1[index]
          y2.b<- y2[index]
          if(is.vector(x)) x.b<- x[index] else x.b<- x[index,]
          fit1<- glm(y1.b~ x.b, family = binomial(link =link1))
          fit2<- glm(y2.b~ x.b, family = binomial(link =link2))
          pcor.est[i]<- Pcor_SR(y1.b, y2.b, x.b, model1 = fit1, model2 = fit2, link1=link1, link2=link2, n.avg=n.avg, method=method)$phi
        }, error=function(e){})
      }
    }
  }
  Bstar<- sum(!is.na(pcor.est))
  pcor.est1<- pcor.est[!is.na(pcor.est)]

  if(H0==0){
    pval<- 2*min(mean(pcor.est1>0), mean(pcor.est1<0))
  }else{
    pval<- 2*min(mean(pcor.est1<H0), mean(pcor.est1> -1*H0), 0.5)
  }

  return(list(est=mean(pcor.est1), sd=sd(pcor.est1), pval=pval, CI=quantile(pcor.est1, probs = c(0.025, 0.975)), Bcount=Bstar))
}

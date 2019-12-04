#' Pcor_SR
#'
#' A function to obtain partial associations between several ordinal variables
#' after adjusting for a set of covariates.
#'
#' @param models a list contains all the models (S3 objects) you want to
#' assess for the partial association between ordinal responses after adjusting
#' for a set of covariates covariates.
#' They all should be applied to same dataset, having same covariates, same sample size etc.
#' The models in this list can be an object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#' @param asso_type a character string indicating whether the partial association coefficient
#' \code{"PartialAsso"} (default) or the marginal correlation \code{"Marginal"} is to be
#' computed. They can be abbreviated.
#' @param S.E a logical argument indicates whether conduct bootstrap for standard error.
#' It can be very slow so a double-check question should be answered before move on.
#' @param boot_S.E the bootstrap times for calculating the standard error of the partial
#' association pho.
#' @param rep_num the repeat times for getting the average as estimator of partial
#' association pho. It is the \code{"nsim"} in \code{"resids()"} function
#' @param ... Additional optional arguments. (say \code{"method"} and \code{"jitter.scale"}
#' in \code{"resids"} function. Just keep them as default if not specified. Currently ignored.)
#' @param cor_method correlation coefficient method
#'
#' @return A list contains ...
#' \itemize{
##'  \item{"Mcor_phi"}{a correlation matrix of all responses in these models}
##'  \item{"SRs"}{a matrix contains all surrogate residuals of each model}
##'  \item{"Mboot_phi"}{Stuff}
##'  \item{"boot_SRs"}{an array contains all bootstrap surrogate residuals of each model.
##'  First dimension is sample size, second is bootstrap times, third is number of models.}
##' }
#' @export
#'
Pcor_SR <- function(models, asso_type = c("PartialAsso", "Marginal"),
                   cor_method = c("kendall", "pearson", "wolfsigma"),
                   rep_num=10, S.E=FALSE, boot_S.E=100, ...){
  # model1 = fit.vote; model2 = fit.PID
  # models <- list(model1, model2)
  # num_models <- length(models)
  # models
  # class(models[[1]])
  # class(model1)
  #
  # asso_type = "PartialAsso"; cor_method = "kendall";
  # rep_num = 30; boot_S.E = 100;


  # distributions <- sapply(models, function(mod) getDistributionName(mod))

  # Match arguments -------------------------------------------------
  if (missing(asso_type)) { asso_type <- "PartialAsso" }
  asso_type <- match.arg(asso_type)

  if (missing(cor_method)) { cor_method <- "kendall" }
  cor_method <- match.arg(cor_method)

  # Initialize the cor_func for different methods
  cor_func <- switch(cor_method,
                     kendall = function(X) cor(X, method = "kendall"),
                     pearson = function(X) cor(X),
                     wolfsigma = function(X) t(copBasic::wolfCOP(para = data.frame(X),
                                                       as.sample = TRUE)))
  if (S.E) {
    StillBoot <- utils::menu(c("Yes", "No"), title="Using bootstrap to calculate standard errors of phi could be slow. \n
         Do you really want to continue?")
    if (StillBoot) { # If still run bootstrap for S.E, throw messages, otherwise stop.
      S.E <- TRUE
      message("Runinng...")
    } else {
      stop("Please change the S.E to FALSE and do it again!", call. = FALSE)
    }
  }

  # Warnings --------------------------------------------------------
  mods_n <- sapply(models, nobs)
  mods_ys <- sapply(models, function(mod) colnames(model.frame(mod))[1])
  # mods_xs <- sapply(models, function(mod) all.vars(mod$formula))

  if (!all(mods_n==mods_n[1])) { # if data are not
    stop("Stop due to different data in the models!")
  }
  # if (all(mods_xs==mods_xs[1])){ # [FEATRUE]: This need to check in the future!
  #   stop("Stop due to different controlling covariates in these models!")
  # }

  # Pull out ingredients for cooking --------------------------------------------
  mods_y <- sapply(models, function(mod) model.frame(mod)[,1])

  # Main Body ---------------------------------------------------------------
  if (asso_type=="Marginal") { # marginal association (Default Kendall's tau)

    Mcor_phi <- cor(mods_y, method = cor_method)

    if (S.E) { # [FEATURE]: standard error from bootstrap, leave it in the future!
      # phi_sd_boot<- sd(sapply(1:boot_S.E, function(b){
      #   index <- sample(n, replace=T)
      #   cor(y1[index], y2[index], method = cor_method)
      # }))
    }
  } else { # Partial Association based on surrogate residuals

    # Only simulate once to get correlation
    SRs <- sapply(models, function(mod)
      resids(mod, method = c("latent"), nsim = 1))

    Mcor_phi <- cor_func(SRs)

    # [BUG]: bootstrap "rep_num" times and use mean as estimator of phi,
    #        the results are problematic, leave it in the future!
    boot_SRs <- sapply(models, function(mod)
      attr(resids(mod, method = c("latent"), nsim = rep_num), "boot_reps"),
      simplify = "array")
    # dim(boot_SRs)

    Mboot_phi_temp <- apply(X = boot_SRs, MARGIN = 2, cor_func)
    Mboot_phi <- matrix(rowMeans(Mboot_phi_temp),
                        nrow = length(mods_ys), ncol = length(mods_ys),
                        dimnames = list(mods_ys, mods_ys))
    # dim(Mboot_phi)

    if (S.E & !missing(boot_S.E)) { # boot_S.E=2000
      message("Under Developing!")
    }
  }

  # Give names for correlation matrix and "boot_SRs" array
  colnames(Mcor_phi) <- rownames(Mcor_phi) <- colnames(SRs) <- mods_ys

  dimnames(boot_SRs)[[3]] <- mods_ys

  return(
    list(Mcor_phi=Mcor_phi, SRs=SRs, Mboot_phi=Mboot_phi, boot_SRs=boot_SRs)
  )

} ## end of function

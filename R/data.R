#' A list from an adjacent categories model
#'
#' This list contains dataset and model coefficients. This example is used to
#' illustrate the association between the residual variables when the
#' two ordinal variables Y1 and Y2 are partially independent. Data simulated from an
#' adjacent categories regression model with an ordered (preferably) factor response.
#' beta1 = 1, beta2 = -1,
#' alpha1 = (-Inf, -3, -2, 0, 2, 3, Inf), alpha2 = (-Inf, -2, 0, 2, Inf)
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 10000 rows and 11 variables.
#' \itemize{
#'   \item \code{X} A continuous predictor variable.
#'   \item \code{Y1} The response variable; an ordered factor.
#'   \item \code{Y2} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang, Li, Shaobo, Yu, Yan, and Moustaki, Irini. Assessing partial association between
#' ordinal variables: quantification, visualization, and hypothesis testing, \emph{Journal of the
#' American Statistical Association}, Revision under review.
#'
#' @name df_AdjCat
#'
#' @usage
#' data(df_AdjCat)
#'
#' @examples
#' #
#' # Adjacent Categories Regression Model Example to compare different residuals.
#' # After adjusting covariates, the association between residuals variables should be independdent.
#' # Surrogate residuals has this property, whereas other types of residuals do not.
#' #
#'
#' data("df_AdjCat")
#' summary(df_AdjCat$data)
#' fit_clm1 <- VGAM::vglm(Y1 ~ X, family =
#'                       VGAM::cumulative(link = "logit",reverse=TRUE,parallel = TRUE),
#'                       data = df_AdjCat$data)
#' fit_clm2 <- VGAM::vglm(Y2 ~ X, family =
#'                        VGAM::cumulative(link = "logit",reverse=TRUE,parallel = TRUE),
#'                        data = df_AdjCat$data)
#' SR1 <- residuals(object = fit_clm1, type = "surrogate", surr.method = "latent", boot_id = NULL)
#' SR2 <- residuals(fit_clm2, type = "surrogate", surr.method = "latent", boot_id = NULL)
#'
#' ## obtain SBC residuals (Li and Shepherd 2012 JASA/Biometrika)
#' PR1 <- residuals(fit_clm1, type = "sign", boot_id = NULL)
#' PR2 <- residuals(fit_clm2, type = "sign", boot_id = NULL)
#'
#' ## obtain generalized residuals (Franses and Paap 2001 book)
#' GR1 <- residuals(fit_clm1, type = "general", boot_id = NULL)
#' GR2 <- residuals(fit_clm2, type = "general", boot_id = NULL)
#'
#' ## obtain deviance residuals
#' DR1 <- residuals(fit_clm1, type = "deviance", boot_id = NULL)
#' DR2 <- residuals(fit_clm2, type = "deviance", boot_id = NULL)
#'
#' ## visualize residual vs. residual
#' par(mfrow=c(2,2))
#' par(mar=c(4, 4.8, 2.5, 1.5))
#'
#' plot(PR1, PR2, pch=".", main = "sign-based Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")))
#' plot(GR1, GR2, pch=".", main = "generalized Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")), xlim = c(-4,4), ylim=c(-4,4))
#' plot(DR1, DR2, pch='.', main = "deviance Residuals",
#'      xlab = expression(paste(R[1]^"ALT")),
#'      ylab = expression(paste(R[2]^"ALT")))
#' plot(SR1, SR2, pch=".", main = "Surrogate Residuals", xaxt="n", yaxt="n",
#'      xlab = expression(R[1]), ylab = expression(R[2]),
#'      xlim = c(-1/2,1/2), ylim=c(-1/2,1/2))
#' axis(1, at=seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))
#' axis(2, at=seq(-0.5, 0.5, 0.25), labels = seq(-0.5, 0.5, 0.25))
#'
NULL

#' US 2016 national election study with pre-election interview only (Clean)
#'
#' A subset of 2,188 participants of the 2016 American National Election Time Series Study,
#' which was to track the enduring social trend and record the political moment of 2016 (DeBell, 2018).
#' This study consisted of two surveys with same population. The pre-election interview
#' was during the weeks before the 2016 general election, including 4,271 respondents in
#' total. The post-election interview is the re-interview during the weeks after the
#' election, including 3,649 respondents (662 respondents did not complete post-interviews).
#'
#' The Pre-election preference is recorded as "IntendVote" and the "Prevote.num"
#' is the numeric of it. Observations with missing values, or "No thought"
#' responses have been removed. Responsdents expressing a voting preference
#' other than Clinton or Trump have been removed.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2188 rows and 10 variables.
#' \itemize{
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{edu.year} Respondent's education year, which is mapped from \code{education}
#'   level: \code{MS}=8, \code{HSdrop}=11, \code{HS}=12, \code{Coll}=14,
#'   \code{CCdeg}=15, \code{BAdeg}=17, \code{MAdeg}=19.
#'
#'   \item \code{education} Respondent's education level.
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'   It is median value of the range of each \code{income} level
#'
#'   \item \code{income.num} Respondent's family income level:
#'
#'   '(01) 01. Under $5,000' = 5,
#'
#'   '(02) 02. $5,000-$9,999',
#'
#'   '(03) 03. $10,000-$12,499',
#'
#'   '(04) 04. $12,500-$14,999',
#'
#'   '(05) 05. $15,000-$17,499',
#'
#'   '(06) 06. $17,500-$19,999',
#'
#'   '(07) 07. $20,000-$22,499',
#'
#'   '(08) 08. $22,500-$24,999',
#'
#'   '(09) 09. $25,000-$27,499',
#'
#'   '(10) 10. $27,500-$29,999',
#'
#'   '(11) 11. $30,000-$34,999',
#'
#'   '(12) 12. $35,000-$39,999',
#'
#'   '(13) 13. $40,000-$44,999',
#'
#'   '(14) 14. $45,000-$49,999',
#'
#'   '(15) 15. $50,000-$54,999',
#'
#'   '(16) 16. $55,000-$59,999',
#'
#'   '(17) 17. $60,000-$64,999',
#'
#'   '(18) 18. $65,000-$69,999',
#'
#'   '(19) 19. $70,000-$74,999',
#'
#'   '(20) 20. $75,000-$79,999',
#'
#'   '(21) 21. $80,000-$89,999',
#'
#'   '(22) 22. $90,000-$99,999',
#'
#'   '(23) 23. $100,000-$109,999',
#'
#'   '(24) 24. $110,000-$124,999',
#'
#'   '(25) 25. $125,000-$149,999',
#'
#'   '(26) 26. $150,000-$174,999',
#'
#'   '(27) 27. $175,000-$249,999',
#'
#'   '(28) 28. $250,000 or more' = 250.
#'
#'   \item \code{PID} Party identification: a numeric variable with value from 1 to 7
#'   representing strong Democrat, \code{strDem} < weak Democrat, \code{weakDem} <
#'   independent Democrat, \code{indDem} < independent independent \code{indind} <
#'   indepedent Republican, \code{indRep} < weak Republican, \code{weakRep} <
#'   strong Republican, \code{strRep}.

#'   \item \code{selfLR} The respondent' self-placement about own left-right in 7 ordinal levels
#'   (from extremely liberal to extremely conservative). \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{TrumpLR} The respondent's opinion about Donald Trump's left-right
#'   placement (same scale as selfLR).
#'
#'   \item \code{ClinLR} The respondent's opinion about Hilary Clinton's left-right
#'   placement (same scale as selfLR).
#'
#'   \item \code{IntendVote} The respondent's voting preference between
#'   Donald Trump and Hilary Clinton two months preceeding the
#'   November election (Pre-election interview). It is a factor with levels
#'   \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Prevote.num} Recode the intend vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{WeightforPreElection} Pre-election weight of a respondent.
#'
#' }
#'
#' @references
#' DeBell, Matthew, Jon A. Krosnick, Katie Gera, David S. Yeager, and Michael P. McDonald.
#' The turnout gap in surveys: Explanations and solutions. \emph{Sociological Methods & Research}, 2018:
#' 0049124118769085.
#'
#' @name nes2016
#'
#' @usage
#' data(nes2016)
#'
#' @examples
#' head(nes2016)
NULL

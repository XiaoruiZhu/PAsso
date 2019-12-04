#' Simulated quadratic data
#'
#' Data simulated from a probit model with a quadratic trend. The data are
#' described in Example 2 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df1
#'
#' @usage
#' data(df1)
#'
#' @examples
#' head(df1)
NULL


#' Simulated heteroscedastic data
#'
#' Data simulated from a probit model with heteroscedasticity. The data are
#' described in Example 4 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df2
#'
#' @usage
#' data(df2)
#'
#' @examples
#' head(df2)
NULL


#' Simulated Gumbel data
#'
#' Data simulated from a log-log model with a quadratic trend. The data are
#' described in Example 3 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df3
#'
#' @usage
#' data(df3)
#'
#' @examples
#' head(df3)
NULL


#' Simulated proportionality data
#'
#' Data simulated from from two separate probit models. The data are described
#' in Example 5 of Liu and Zhang (2017).
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 4000 rows and 2 variables.
#' \itemize{
#'   \item \code{x} The predictor variable.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df4
#'
#' @usage
#' data(df4)
#'
#' @examples
#' head(df4)
NULL


#' Simulated interaction data
#'
#' Data simulated from from an ordered probit model with an interaction effect.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 2000 rows and 3 variables.
#' \itemize{
#'   \item \code{x1} A continuous predictor variable.
#'   \item \code{x2} A factor with two levels: \code{"Control"} and
#'   \code{"Treatment"}.
#'   \item \code{y} The response variable; an ordered factor.
#' }
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association} (accepted).
#'
#' @name df5
#'
#' @usage
#' data(df5)
#'
#' @examples
#' head(df5)
NULL


#' Ordinal Response with Adjacent Categories Probabilities
#'
#' Data simulated from an adjacent categories regression model
#' with an ordered (preferably) factor response.
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
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association}.
#'
#' @name df_ParA
#'
#' @usage
#' data(df_ParA)
#'
#' @examples
#' head(df_ParA)
NULL



#' US 1996 national election study
#'
#' A data with 10 variable subset of the 1996 American National
#' Election Study.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 944 rows and 13 variables.
#' \itemize{
#'   \item \code{popul} population of respondent's location in 1000s of people.
#'
#'   \item \code{TVnews} days in the past week spent watching news on TV.
#'
#'   \item \code{selfLR} Left-Right self-placement of respondent: an ordered
#'   factor with levels \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{ClinLR} Left-Right placement of Bill Clinton (same
#'   scale as selfLR): an ordered factor with levels \code{extLib} <
#'   \code{Lib} < \code{sliLib} < \code{Mod} < \code{sliCon} <
#'   \code{Con} < \code{extCon}.
#'
#'   \item \code{PID} Party identification.
#'
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{educ} Respondent's education,
#'
#'   \item \code{edu.year} Respondent's education (Recoded), edu.year=
#'   recode(educ, 'MS'=8, 'HSdrop'=11, 'HS'=12, 'Coll'=14, 'CCdeg'=15,
#'   'BAdeg'=17, 'MAdeg'=19).
#'
#'   \item \code{income.num} Respondent's family income: an ordered
#'   factor with levels: income.num=recode(income, '$3Kminus'=3, '$3K-$5K'=4,
#'   '$5K-$7K'=6, '$7K-$9K'=8, '$9K-$10K'=9.5, '$10K-$11K'=10.5,
#'   '$11K-$12K'=11.5, '$12K-$13K'=12.5, '$13K-$14K'=13.5, '$14K-$15K'=14.5,
#'   '$15K-$17K'=16, '$17K-$20K'=18.5, '$20K-$22K'=21, '$22K-$25K'=23.5,
#'   '$25K-$30K'=27.5, '$30K-$35K'=32.5, '$35K-$40K'=37.5, '$40K-$45K'=42.5,
#'   '$45K-$50K'=47.5, '$50K-$60K'=55, '$60K-$75K'=67, '$75K-$90K'=80,
#'   '$90K-$105K'=95, '$105Kplus'=105).
#'
#'   \item \code{vote} Expected vote in 1996 presidential election:
#'   a factor with levels \code{Clinton} and \code{Dole}.
#'   vote.num=recode(vote, 'Clinton'=0, 'Dole'=1),
#' }
#'
#' @references
#'
#' @name nes96clean
#'
#' @usage
#' data(nes96clean)
#'
#' @examples
#' head(nes96clean)
NULL

#' Ordinal Response with Adjacent Categories Probabilities
#'
#' Data simulated from an adjacent categories regression model
#' with an ordered (preferably) factor response. beta1 = 1, beta2 = -1,
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
#' Liu, Dungang, Li, Shaobo, Yu, Yu, and Moustaki, Irini. Assessing partial association between
#' ordinal variables: quantification, visualization, and hypothesis testing, \emph{Journal of the
#' American Statistical Association}, Revision under review.
#'
#' @name df_AdjCat
#'
#' @usage
#' data(df_AdjCat)
#'
#' @examples
#' head(df_AdjCat)
NULL

#' Raw dataset of US 2016 national election study
#'
#' A subset of the 2016 American National Election Study.
#' Pre-election preference is recorded as "IntendVote", while the actual vote is "voteResult".
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 4271 rows and 12 variables.
#' \itemize{
#'   \item \code{WeightforPreElection} Pre-election weight of a respondent.
#'
#'   \item \code{WeightforPostElection} Post-election weight of a respondent.
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
#'   \item \code{PID} Party identification: an ordered factor with levels strong
#'   Democrat, \code{strDem} < weak Democrat, \code{weakDem} < independent Democrat,
#'   \code{indDem} < independent independent \code{indind} < indepedent Republican,
#'   \code{indRep} < weak Republican, \code{weakRep} < strong Republican, \code{strRep}.
#'
#'   \item \code{agegroup} Respondent's age category.
#'
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{educ} Respondent's education category: an ordered factor with levels
#'   8 years or less, \code{MS} < high school dropout, \code{HSdrop} < high school
#'   diploma or GED, \code{HS} < some College, \code{Coll} < Community or junior
#'   College degree, \code{CCdeg} < BA degree, \code{BAdeg} < postgraduate degree,
#'   \code{MAdeg}.
#'
#'   \item \code{income} Categorical variable representing the Respondent's family income
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'   income.num = recode(income, '(01) 01. Under $5,000' = 5, '(02) 02. $5,000-$9,999' = 7.5,
#'   '(03) 03. $10,000-$12,499' = 11.25, '(04) 04. $12,500-$14,999' = 13.75,
#'   '(05) 05. $15,000-$17,499' = 16.25, '(06) 06. $17,500-$19,999' = 18.75,
#'   '(07) 07. $20,000-$22,499' = 21.25, '(08) 08. $22,500-$24,999' = 23.75,
#'   '(09) 09. $25,000-$27,499' = 26.25, '(10) 10. $27,500-$29,999' = 28.75,
#'   '(11) 11. $30,000-$34,999' = 32.5, '(12) 12. $35,000-$39,999' = 37.5,
#'   '(13) 13. $40,000-$44,999' = 42.5, '(14) 14. $45,000-$49,999' = 47.5,
#'   '(15) 15. $50,000-$54,999' = 52.5, '(16) 16. $55,000-$59,999' = 57.5,
#'   '(17) 17. $60,000-$64,999' = 62.5, '(18) 18. $65,000-$69,999' = 67.5,
#'   '(19) 19. $70,000-$74,999' = 72.5, '(20) 20. $75,000-$79,999' = 77.5,
#'   '(21) 21. $80,000-$89,999' = 85, '(22) 22. $90,000-$99,999' = 95,
#'   '(23) 23. $100,000-$109,999' = 105, '(24) 24. $110,000-$124,999' = 117.5,
#'   '(25) 25. $125,000-$149,999' = 137.5, '(26) 26. $150,000-$174,999' = 162.5,
#'   '(27) 27. $175,000-$249,999' = 212.5, '(28) 28. $250,000 or more' = 250)
#'
#'   \item \code{IntendVote} The intend vote two months preceeding the
#'   November election (Pre-election interview). It is a factor with levels
#'   \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Prevote.num} Recode the intend vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{voteResult} The actual vote answered in the interview during the
#'   two months following the election (Post-election interview). It is a factor
#'   with levels \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Postvote.num} Recode the actual vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#' }
#'
#'
#' @name nes2016_raw
#'
#' @usage
#' data(nes2016_raw)
#'
#' @examples
#' head(nes2016_raw)
NULL

#' US 2016 national election study with pre-election interview only
#'
#' A subset of the 2016 American National Election Study.
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
#'   \item \code{edu.year} Respondent's education year, which is mapped from raw data
#'   \code{nes2016_raw}. \code{MS}=8, \code{HSdrop}=11, \code{HS}=12, \code{Coll}=14,
#'   \code{CCdeg}=15, \code{BAdeg}=17, \code{MAdeg}=19.
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'   income.num = recode(income, '(01) 01. Under $5,000' = 5, '(02) 02. $5,000-$9,999' = 7.5,
#'   '(03) 03. $10,000-$12,499' = 11.25, '(04) 04. $12,500-$14,999' = 13.75,
#'   '(05) 05. $15,000-$17,499' = 16.25, '(06) 06. $17,500-$19,999' = 18.75,
#'   '(07) 07. $20,000-$22,499' = 21.25, '(08) 08. $22,500-$24,999' = 23.75,
#'   '(09) 09. $25,000-$27,499' = 26.25, '(10) 10. $27,500-$29,999' = 28.75,
#'   '(11) 11. $30,000-$34,999' = 32.5, '(12) 12. $35,000-$39,999' = 37.5,
#'   '(13) 13. $40,000-$44,999' = 42.5, '(14) 14. $45,000-$49,999' = 47.5,
#'   '(15) 15. $50,000-$54,999' = 52.5, '(16) 16. $55,000-$59,999' = 57.5,
#'   '(17) 17. $60,000-$64,999' = 62.5, '(18) 18. $65,000-$69,999' = 67.5,
#'   '(19) 19. $70,000-$74,999' = 72.5, '(20) 20. $75,000-$79,999' = 77.5,
#'   '(21) 21. $80,000-$89,999' = 85, '(22) 22. $90,000-$99,999' = 95,
#'   '(23) 23. $100,000-$109,999' = 105, '(24) 24. $110,000-$124,999' = 117.5,
#'   '(25) 25. $125,000-$149,999' = 137.5, '(26) 26. $150,000-$174,999' = 162.5,
#'   '(27) 27. $175,000-$249,999' = 212.5, '(28) 28. $250,000 or more' = 250)
#'
#'   \item \code{PID} Party identification: a numeric variable with value from 1 to 7
#'   representing strong Democrat, \code{strDem} < weak Democrat, \code{weakDem} <
#'   independent Democrat, \code{indDem} < independent independent \code{indind} <
#'   indepedent Republican, \code{indRep} < weak Republican, \code{weakRep} <
#'   strong Republican, \code{strRep}.

#'   \item \code{selfLR} Left-Right self-placement of respondent: a numeric variable
#'   with value from 1 to 7 representing \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{TrumpLR} Left-Right placement of Donald Trump (same
#'   scale as selfLR), a numeric variable with value from 1 to 7.
#'
#'   \item \code{ClinLR} Left-Right placement of Bill Clinton (same
#'   scale as selfLR): a numeric variable with value from 1 to 7.
#'
#'   \item \code{IntendVote} The intend vote two months preceeding the
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
#'
#' @name nes2016_pre
#'
#' @usage
#' data(nes2016_pre)
#'
#' @examples
#' head(nes2016_pre)
NULL

#' US 2016 national election study with post-election interview only
#'
#' A subset of the 2016 American National Election Study. The actual vote is "voteResult".
#' "Postvote.num" is recoded variable of "voteResult". Observations with missing values,
#' or "No thought" responses have been removed. Responsdents expressing
#' a voting preference other than Clinton or Trump have been removed.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 1892 rows and 10 variables.
#' \itemize{
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{edu.year} Respondent's education year, which is mapped from raw data
#'   \code{nes2016_raw}. \code{MS}=8, \code{HSdrop}=11, \code{HS}=12, \code{Coll}=14,
#'   \code{CCdeg}=15, \code{BAdeg}=17, \code{MAdeg}=19.
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'
#'   \item \code{PID} Party identification: a numeric variable with value from 1 to 7
#'   representing strong Democrat, \code{strDem} < weak Democrat, \code{weakDem} <
#'   independent Democrat, \code{indDem} < independent independent \code{indind} <
#'   indepedent Republican, \code{indRep} < weak Republican, \code{weakRep} <
#'   strong Republican, \code{strRep}.

#'   \item \code{selfLR} Left-Right self-placement of respondent: a numeric variable
#'   with value from 1 to 7 representing \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{TrumpLR} Left-Right placement of Donald Trump (same
#'   scale as selfLR), a numeric variable with value from 1 to 7.
#'
#'   \item \code{ClinLR} Left-Right placement of Bill Clinton (same
#'   scale as selfLR): a numeric variable with value from 1 to 7.
#'
#'   \item \code{voteResult} The actual vote answered in the interview during the
#'   two months following the election (Post-election interview). It is a factor
#'   with levels \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Postvote.num} Recode the actual vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{WeightforPostElection} Post-election weight of a respondent.
#'
#' }
#'
#'
#' @name nes2016_post
#'
#' @usage
#' data(nes2016_post)
#'
#' @examples
#' head(nes2016_post)
NULL

#' US 2016 national election study of respondents without missing in both pre-election and post-election.
#'
#' A subset of the 2016 American National Election Study. Pre-election preference is recorded
#' as "IntendVote", while the actual vote is "voteResult".
#' "Prevote.num" and "Postvote.num" are the numeric for them, respectively. Observations
#' with missing values, or "No thought" responses have been removed. Responsdents expressing
#' a voting preference other than Clinton or Trump have been removed.
#'
#' @docType data
#'
#' @keywords datasets
#'
#' @format A data frame with 1690 rows and 14 variables.
#' \itemize{
#'   \item \code{age} Respondent's age in years.
#'
#'   \item \code{edu.year} Respondent's education year, which is mapped from raw data
#'   \code{nes2016_raw}. \code{MS}=8, \code{HSdrop}=11, \code{HS}=12, \code{Coll}=14,
#'   \code{CCdeg}=15, \code{BAdeg}=17, \code{MAdeg}=19.
#'
#'   \item \code{income.num} Respondent's family income in thousands: an numerical variable.
#'
#'   \item \code{PID} Party identification: a numeric variable with value from 1 to 7
#'   representing strong Democrat, \code{strDem} < weak Democrat, \code{weakDem} <
#'   independent Democrat, \code{indDem} < independent independent \code{indind} <
#'   indepedent Republican, \code{indRep} < weak Republican, \code{weakRep} <
#'   strong Republican, \code{strRep}.

#'   \item \code{selfLR} Left-Right self-placement of respondent: a numeric variable
#'   with value from 1 to 7 representing \code{extLib}: extremely liberal, \code{Lib}: liberal,
#'   \code{sliLib}: slightly liberal, \code{Mod}: moderate, \code{sliCon}:
#'   slightly conservative, \code{Con}: conservative, \code{extCon}: extremely
#'   conservative. \code{extLib} < \code{Lib} < \code{sliLib} <
#'   \code{Mod} < \code{sliCon} < \code{Con} < \code{extCon}.
#'
#'   \item \code{TrumpLR} Left-Right placement of Donald Trump (same
#'   scale as selfLR), a numeric variable with value from 1 to 7.
#'
#'   \item \code{ClinLR} Left-Right placement of Bill Clinton (same
#'   scale as selfLR): a numeric variable with value from 1 to 7.
#'
#'   \item \code{IntendVote} The intend vote two months preceeding the
#'   November election (Pre-election interview). It is a factor with levels
#'   \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Prevote.num} Recode the intend vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{voteResult} The actual vote answered in the interview during the
#'   two months following the election (Post-election interview). It is a factor
#'   with levels \code{HillaryClinton} and \code{DonaldTrump}.
#'
#'   \item \code{Postvote.num} Recode the actual vote to numeric values,
#'   'HillaryClinton'=0, 'DonaldTrump'=1.
#'
#'   \item \code{WeightforPreElection} Pre-election weight of a respondent.
#'
#'   \item \code{WeightforPostElection} Post-election weight of a respondent.
#'
#' }
#'
#'
#' @name nes2016_prepost
#'
#' @usage
#' data(nes2016_prepost)
#'
#' @examples
#' head(nes2016_prepost)
NULL

#' US 1996 national election study
#'
#' A data with 13 variables subset of the 1996 American National Election Study.
#' The data has been cleaned, and a few numeric variables are added.
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
#'   \item \code{income} Respondent's family income
#'
#'   \item \code{income.num} Respondent's family income: an ordered
#'   factor with levels: income.num=recode(income, '3Kminus'=3, '3K-$5K'=4,
#'   '5K−7K'=6, '7K−9K'=8, '9K−10K'=9.5, '10K−11K'=10.5,
#'   '11K−12K'=11.5, '12K−13K'=12.5, '13K−14K'=13.5, '14K−15K'=14.5,
#'   '15K−17K'=16, '17K−20K'=18.5, '20K−22K'=21, '22K−25K'=23.5,
#'   '25K−30K'=27.5, '30K−35K'=32.5, '35K−40K'=37.5, '40K−45K'=42.5,
#'   '45K−50K'=47.5, '50K−60K'=55, '60K−75K'=67, '75K−90K'=80,
#'   '90K−105K'=95, '$105Kplus'=105).
#'
#'   \item \code{vote} Expected vote in 1996 presidential election:
#'   a factor with levels \code{Clinton} and \code{Dole}.
#'
#'   \item \code{vote.num} recode(vote, 'Clinton'=0, 'Dole'=1),
#' }
#'
#'
#' @name nes96
#'
#' @usage
#' data(nes96)
#'
#' @examples
#' head(nes96)
NULL

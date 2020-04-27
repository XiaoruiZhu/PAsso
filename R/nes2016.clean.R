#' nes2016pre.clean
#' is a function to clean the nes2016raw and return the cleaned dataset
#' with the pre-election interview only.
#'
#' @param interview A string to specify the obsertations that want to be include. Default is "pre",
#' which only keep the observations without missing in the pre-election interview. "post" will return
#' observations without missing in the post-election interview. If interview equals to "both",
#' returned dataset only contains the observations without missing in both the pre and post-election
#' interviews.
#'
#' @return a dataset with pre-election intend vote and demographic information.
#' @export
#'
#' @examples
#' nes2016_pre <- nes2016pre.clean()
#'
nes2016.clean <- function(interview=c("pre", "post", "both")){
  if (missing(interview)) { interview <- "pre" }
  interview <- match.arg(interview)

  nes2016_temp <- PAsso::nes2016_raw %>%
    dplyr::mutate(Postvote.num=recode(voteResult, 'HillaryClinton'=0, 'DonaldTrump'=1),
                  Prevote.num=recode(IntendVote, 'HillaryClinton'=0, 'DonaldTrump'=1),
                  PID=recode(PID, 'strDem'=1, 'weakDem'=2, 'indDem'=3, 'indind'=4, 'indRep'=5, 'weakRep'=6, 'strRep'=7),
                  selfLR=recode(selfLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
                  ClinLR=recode(ClintonLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
                  TrumpLR=recode(TrumpLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
                  edu.year=recode(education, 'MS'=8, 'HSdrop'=11, 'HS'=12, 'Coll'=14, 'CCdeg'=15, 'BAdeg'=17, 'MAdeg'=19),
                  income.num=recode(income,
                                    '(01) 01. Under $5,000' = 5,
                                    '(02) 02. $5,000-$9,999' = 7.5,
                                    '(03) 03. $10,000-$12,499' = 11.25,
                                    '(04) 04. $12,500-$14,999' = 13.75,
                                    '(05) 05. $15,000-$17,499' = 16.25,
                                    '(06) 06. $17,500-$19,999' = 18.75,
                                    '(07) 07. $20,000-$22,499' = 21.25,
                                    '(08) 08. $22,500-$24,999' = 23.75,
                                    '(09) 09. $25,000-$27,499' = 26.25,
                                    '(10) 10. $27,500-$29,999' = 28.75,
                                    '(11) 11. $30,000-$34,999' = 32.5,
                                    '(12) 12. $35,000-$39,999' = 37.5,
                                    '(13) 13. $40,000-$44,999' = 42.5,
                                    '(14) 14. $45,000-$49,999' = 47.5,
                                    '(15) 15. $50,000-$54,999' = 52.5,
                                    '(16) 16. $55,000-$59,999' = 57.5,
                                    '(17) 17. $60,000-$64,999' = 62.5,
                                    '(18) 18. $65,000-$69,999' = 67.5,
                                    '(19) 19. $70,000-$74,999' = 72.5,
                                    '(20) 20. $75,000-$79,999' = 77.5,
                                    '(21) 21. $80,000-$89,999' = 85,
                                    '(22) 22. $90,000-$99,999' = 95,
                                    '(23) 23. $100,000-$109,999' = 105,
                                    '(24) 24. $110,000-$124,999' = 117.5,
                                    '(25) 25. $125,000-$149,999' = 137.5,
                                    '(26) 26. $150,000-$174,999' = 162.5,
                                    '(27) 27. $175,000-$249,999' = 212.5,
                                    '(28) 28. $250,000 or more' = 250))

  if (interview == "pre") {
    nes2016 <- nes2016_temp %>%
      dplyr::filter(!is.na(Prevote.num) & (Prevote.num != 3)) %>% # Only keep obs of two votes in prevote
      dplyr::select(-c(education, agegroup, income, voteResult,
                ClintonLR, Postvote.num, WeightforPostElection)) %>% # Remove useless variables
      dplyr::select(c(age, edu.year, income.num,
               PID, selfLR, TrumpLR, ClinLR, IntendVote, Prevote.num,
               WeightforPreElection))
  } else if (interview == "post") {
    nes2016 <- nes2016_temp %>%
      dplyr::filter(!is.na(Postvote.num) & (Postvote.num != 3)) %>% # Only keep obs of two votes in postvote
      dplyr::select(-c(education, agegroup, income, IntendVote,
                ClintonLR,
                Prevote.num, WeightforPreElection)) %>% # Remove useless variables
      dplyr::select(c(age, edu.year, income.num,
               PID, selfLR, TrumpLR, ClinLR, voteResult, Postvote.num,
               WeightforPostElection))
  } else { # Keep both prevote and postvote
    nes2016 <- nes2016_temp %>%
      dplyr::filter(!is.na(Postvote.num) & !is.na(Prevote.num) & (Postvote.num != 3) & (Prevote.num != 3)) %>%
      # Only keep obs of two votes in both prevote and postvote
      dplyr::select(-c(education, agegroup, income, ClintonLR)) %>%
      dplyr::select(c(age, edu.year, income.num,
                      PID, selfLR, TrumpLR, ClinLR,
                      IntendVote, Prevote.num, voteResult, Postvote.num,
                      WeightforPreElection, WeightforPostElection))
  }
  nes2016 %>%
    dplyr::filter(complete.cases(.)) %>% # Only keep obs without missing
    droplevels(.) # Drop useless factors of vote variable IntendVote
}

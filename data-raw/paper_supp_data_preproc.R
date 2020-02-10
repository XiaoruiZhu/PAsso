library(tidyverse)
library(faraway)

## nes96 data pre-processing
data("nes96")

# Frequency table (table 3)
table(nes96$PID, nes96$vote)

nes.clean<- function(){
  nes96 %>%
    mutate(vote.num=recode(vote, 'Clinton'=0, 'Dole'=1),
           PID=recode(PID, 'strDem'=1, 'weakDem'=2, 'indDem'=3, 'indind'=4, 'indRep'=5, 'weakRep'=6, 'strRep'=7),
           selfLR=recode(selfLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
           ClinLR=recode(ClinLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
           DoleLR=recode(DoleLR, 'extLib'=1, 'Lib'=2, 'sliLib'=3, 'Mod'=4, 'sliCon'=5, 'Con'=6, 'extCon'=7),
           edu.year=recode(educ, 'MS'=8, 'HSdrop'=11, 'HS'=12, 'Coll'=14, 'CCdeg'=15, 'BAdeg'=17, 'MAdeg'=19),
           income.num=recode(income, '$3Kminus'=3, '$3K-$5K'=4, '$5K-$7K'=6, '$7K-$9K'=8, '$9K-$10K'=9.5, '$10K-$11K'=10.5,
                             '$11K-$12K'=11.5, '$12K-$13K'=12.5, '$13K-$14K'=13.5, '$14K-$15K'=14.5, '$15K-$17K'=16,
                             '$17K-$20K'=18.5, '$20K-$22K'=21, '$22K-$25K'=23.5, '$25K-$30K'=27.5, '$30K-$35K'=32.5,
                             '$35K-$40K'=37.5, '$40K-$45K'=42.5, '$45K-$50K'=47.5, '$50K-$60K'=55, '$60K-$75K'=67,
                             '$75K-$90K'=80, '$90K-$105K'=95, '$105Kplus'=105))

}


# data pre-processing
nes96 <- nes.clean()

save(nes96, file = "data/nes96.rda", compress='xz')


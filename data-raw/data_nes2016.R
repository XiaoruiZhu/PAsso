library(tidyverse)
library(faraway)
source(file = "data-raw/nes2016.clean.R")

nes2016_raw <- read.csv(file = "D:/Dropbox/1.Research/2.Surrogate/5.R_Packages/PartialAssociation_DL_SL_XZ/3.Data/ANES2016.csv")
summary(nes2016_raw)
nes2016_raw <- nes2016_raw %>% select(-c(X))

save(nes2016_raw, file = "data/nes2016_raw.rda", compress = "xz")

# Frequency table (table 3)
table(nes2016_raw$vote, nes2016_raw$PID)
table(nes2016_raw$IntendVote, nes2016_raw$voteResult)


# pre: Save it as "nes2016", observations without any missing in prevote ------------------------------------------------------------------
nes2016_pre <- nes2016.clean(interview = "pre")
dim(nes2016_pre)
summary(nes2016_pre)
ANES2016 <- nes2016_pre %>% rename(
  PreVote = IntendVote,
  PreVote.num = Prevote.num,
  WeightforPreVote = WeightforPreElection
)

table(ANES2016$PreVote.num)
table(ANES2016$PreVote)
summary(ANES2016)
save(ANES2016, file = "data/ANES2016.rda", compress = "xz")

# Post: Save observations without any missing in postvote -----------------------------------------------------
nes2016_post <- nes2016.clean(interview = "post")

summary(nes2016_post)
table(nes2016_post$Postvote.num)
save(nes2016_post, file = "data/nes2016_post.rda", compress = "xz")

# Both: Save observations without any missing in prevote and postvote -----------------------------------------------------
nes2016_prepost <- nes2016.clean(interview = "both")
dim(nes2016_prepost)
summary(nes2016_prepost)
table(nes2016_prepost$Postvote.num)
save(nes2016_prepost, file = "data/nes2016_prepost.rda", compress = "xz")

# Test
nes2016_test <- nes2016.clean(interview = "post")
dim(nes2016_test)

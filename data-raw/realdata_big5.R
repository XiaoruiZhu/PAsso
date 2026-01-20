setwd("C:/Users/Shaobo/Dropbox/Research/With Dungang/Framework/Rcode")

library(MASS)
library(tidyverse)
library(corrplot)
big5 <- read_csv("Rdata/big5.csv")
big5

big5.clean <- big5 %>%
  filter(age <= 90, country == "US", race == 3) %>%
  mutate(engnat = as.factor(engnat), gender = as.factor(gender), hand = as.factor(hand), source = as.factor(source))
big5.clean

### marginal correlation
t1 <- system.time(corr.tau50 <- cor(big5.clean[, -c(1:7)], method = "kendall"))
corrplot(corr.tau50, method = "color", tl.pos = "n", number.digits = 2, number.cex = 0.45, title = "Marginal association", mar = c(0, 0, 1.5, 0))

# corrplot.mixed(corr.tau50, upper = "color", number.digits = 2, number.cex = 0.45, title="Kendall's tau", mar=c(0,0,1.5,0))

### partial correlation
source("paper_supp_functions.R")

X <- as.data.frame(big5.clean[, 2:6])
sr.mat <- matrix(0, nrow(X), 50)
for (i in 1:50) {
  y <- pull(big5.clean[, 7 + i])
  dat <- data.frame(y, X)
  fit <- polr(as.factor(y) ~ ., data = dat, method = "probit")
  sr.mat[, i] <- SR(y = y, X = model.matrix(~., X)[, -1], model = fit)
}

t2 <- system.time(Pcorr.tau50 <- cor(sr.mat, method = "kendall"))
corrplot(Pcorr.tau50, method = "color", tl.pos = "n", number.digits = 2, number.cex = 0.45, title = "Partial association", mar = c(0, 0, 1.5, 0))

# save(corr.tau50, Pcorr.tau50, file="Rdata/big5_cor_mat.Rdata")
########################################

### compare marginal and partial association
par(mfrow = c(1, 2))
corrplot(corr.tau50, method = "color", tl.pos = "n", number.digits = 2, number.cex = 0.45, title = "Marginal association", mar = c(0, 0, 1.5, 0))
corrplot(Pcorr.tau50, method = "color", tl.pos = "n", number.digits = 2, number.cex = 0.45, title = "Partial association", mar = c(0, 0, 1.5, 0))

# dev.copy2pdf(file="fig-big5-corrplot.pdf")
#################################

### quantify the change of association strength
load("Rdata/big5_cor_mat.Rdata")
change_norm <- change_mapd <- change_mabs <- change_max <- change_01 <- matrix(0, 5, 5)
for (i in 1:5) {
  for (j in i:5) {
    index_row <- (i - 1) * 10 + c(1:10)
    index_col <- (j - 1) * 10 + c(1:10)

    ## mean absolute percent change
    diff_abs <- abs(corr.tau50[index_row, index_col]) - abs(Pcorr.tau50[index_row, index_col])
    index_pos <- diff_abs > 0
    change_mapd[i, j] <- mean(diff_abs[index_pos] / abs(Pcorr.tau50[index_row, index_col][index_pos]))

    ## mean absolute change
    change_mabs[i, j] <- mean(diff_abs[index_pos])

    ## max absolute change
    change_max[i, j] <- max(diff_abs[index_pos])

    ## % change over 0.1
    change_01[i, j] <- mean(diff_abs[index_pos] > 0.1)

    ## Frobenius norm
    change_norm[i, j] <- sum((corr.tau50[index_row, index_col] - Pcorr.tau50[index_row, index_col])^2)
  }
}

aa <- abs(corr.tau50 - Pcorr.tau50)
max(aa)
ind <- which(aa == max(aa), arr.ind = TRUE)
ind
corr.tau50[43, 46]
Pcorr.tau50[43, 46]
###########################

## testing
X <- as.data.frame(big5.clean[, 2:6])
X <- model.matrix(~., X)[, -1]


# using parallel
library(doParallel)
library(doSNOW)
ncore <- detectCores() - 4
cl <- makeCluster(ncore)
registerDoParallel(cl)

pct <- proc.time()
index <- expand.grid(1:10, 1:10)
boot_sd <- foreach(k = 1:100, .export = c("pull", "polr", "pGumbel", "qGumbel")) %dopar% {
  i <- index[k, 2]
  j <- index[k, 1]
  y1 <- pull(big5.clean[, 7 + i]) # E - Extraversion
  # y2<- pull(big5.clean[,7+20+j])  # A - Agreeableness
  y2 <- pull(big5.clean[, 7 + 30 + j]) # C - Concientiousness
  Pcor_SR_boot(y1, y2, X, B = 1000, link1 = "probit", link2 = "probit", n.avg = 5, method = "kendall", H0 = 0)
}
proc.time() - pct
stopCluster(cl)
################################################

cl <- makeCluster(ncore)
registerDoParallel(cl)

pct <- proc.time()
index <- expand.grid(1:10, 1:10)
boot_sd1 <- foreach(k = 1:100, .export = c("pull", "polr", "pGumbel", "qGumbel")) %dopar% {
  i <- index[k, 2]
  j <- index[k, 1]
  y1 <- pull(big5.clean[, 7 + i]) # E - Extraversion
  # y2<- pull(big5.clean[,7+20+j])  # A - Agreeableness
  y2 <- pull(big5.clean[, 7 + 30 + j]) # C - Concientiousness
  Pcor_SR_boot(y1, y2, X, B = 1000, link1 = "probit", link2 = "probit", n.avg = 5, method = "kendall", H0 = 0.1)
}
proc.time() - pct
stopCluster(cl)

# save(boot_sd, boot_sd1, file="Rdata/big5_test_01_E_C.Rdata")
# ## for loop
# pct<- proc.time()
# for(i in 1:10){
#   y1<- pull(big5.clean[,7+i])  # E
#   y2<- pull(big5.clean[,7+20+1])  # A1
#   boot_sd[[i]]<- Pcor_SR_boot(y1, y2, X, B=1000, link1="probit", link2="probit", n.avg=5, method="kendall", H0=0.1)
# }
# proc.time()-pct
# save(boot_sd, file = "Rdata/big5_inf_EA1.Rdata")

# load("Rdata/big5_test_01_E_C.Rdata")
pval0 <- est0 <- std0 <- Bcount0 <- matrix(0, 10, 10)
pval1 <- est1 <- std1 <- Bcount1 <- matrix(0, 10, 10)
for (k in 1:100) {
  i <- index[k, 2]
  j <- index[k, 1]

  est0[i, j] <- round(boot_sd[[k]]$est, 3)
  std0[i, j] <- round(boot_sd[[k]]$sd, 3)
  pval0[i, j] <- round(boot_sd[[k]]$pval, 3)
  Bcount0[i, j] <- round(boot_sd[[k]]$Bcount, 3)

  est1[i, j] <- round(boot_sd1[[k]]$est, 3)
  std1[i, j] <- round(boot_sd1[[k]]$sd, 3)
  pval1[i, j] <- round(boot_sd1[[k]]$pval, 3)
  Bcount1[i, j] <- round(boot_sd1[[k]]$Bcount, 3)
}
rownames(pval0) <- rownames(est0) <- rownames(std0) <- paste("E", 1:10, sep = "")
colnames(pval0) <- colnames(est0) <- colnames(std0) <- paste("C", 1:10, sep = "")
rownames(pval1) <- rownames(est1) <- rownames(std1) <- paste("E", 1:10, sep = "")
colnames(pval1) <- colnames(est1) <- colnames(std1) <- paste("C", 1:10, sep = "")

par(mfrow = c(1, 2))
corrplot((pval0 <= 0.05) * -0.5, method = "circle", number.digits = 2, number.cex = 1, title = "", mar = c(0, 0, 0, 2), cl.pos = "n")
title(main = bquote(atop(paste("Testing of ", phi == 0, " between "), "Extraversion and Concientiousness")), line = 1.5)
corrplot((pval1 <= 0.05) * -0.5, method = "circle", number.digits = 2, number.cex = 1, title = "", mar = c(0, 0, 0, 2), cl.pos = "n")
title(main = bquote(atop(paste("Testing of ", abs(phi) < 0.1, " between "), "Extraversion and Concientiousness")), line = 1.5)

# dev.copy2pdf(file="fig_big5_testing_E_C.pdf")

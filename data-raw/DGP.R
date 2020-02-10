# Data Generating Process for Surrogate Residual --------------------------
ordinalize <- function(z, threshold) {
  sapply(z, FUN = function(x) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && x > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
}


# Function to simulate quadratic data
simQuadraticData <- function(n = 2000) {
  threshold <- c(0, 4, 8)
  x <- runif(n, min = 1, max = 7)
  z <- 16 - 8 * x + 1 * x ^ 2 + rnorm(n)  # rlnorm(n)
  y <- sapply(z, FUN = function(zz) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && zz > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
  data.frame("y" = as.ordered(y), "x" = x)
}

# Simulate data
set.seed(977)
df1 <- simQuadraticData(n = 2000)
table(df1$y)


# Function to simulate heteroscedastic data
simHeteroscedasticData <- function(n = 2000) {
  threshold <- c(0, 30, 70, 100)
  x <- runif(n, min = 2, max = 7)
  y <-   sapply(36 + 4 * x + rnorm(n, sd = x ^ 2), FUN = function(zz) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && zz > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
  data.frame("y" = as.ordered(y), "x" = x)
}

# Simulate heteroscedastic data
set.seed(108)
df2 <- simHeteroscedasticData(n = 2000)
table(df2$y)


# Function to simulate data with Gumbel errors on the linear predictor scale
simGumbelData <- function(n = 2000) {
  x <- runif(n, min = 1, max = 7)
  z <- 16 - 8 * x + 1 * x ^ 2 + sure:::rgumbel(n)
  y <- ordinalize(z, threshold = c(0, 4, 8))
  data.frame("y" = as.ordered(y), "x" = x)
}

# Simulate data
set.seed(977)
df3 <- simGumbelData(n = 2000)
table(df3$y)

# Function to simulate the data from Example 5 in Dungang and Zhang (2017).
simProportionalityData <- function(n = 2000) {
  x <- runif(n, min = -3, max = 3)
  z1 <- 0 - 1 * x + rnorm(n)
  z2 <- 0 - 1.5 * x + rnorm(n)
  y1 <- ordinalize(z1, threshold = c(-1.5, 0))
  y2 <- ordinalize(z2, threshold = c(1, 3))
  data.frame("y" = as.ordered(c(y1, y2)), "x" = c(x, x))
}

# Simulate data
set.seed(977)
df4 <- simProportionalityData(n = 2000)
table(df4$y)


# Function to simulate data from an ordered probit model with an interaction
# term
simInteractionData <- function(n = 2000) {
  threshold <- c(0, 20, 40)
  x1 <- runif(n, min = 1, max = 7)
  x2 <- gl(2, n / 2, labels = c("Control", "Treatment"))
  z <- 16 - 5 * x1 + 3 * (x2 == "Treatment") + 10 * x1 * (x2 == "Treatment") +
    rnorm(n)
  y <- sapply(z, FUN = function(zz) {
    ordinal.value <- 1
    index <- 1
    while(index <= length(threshold) && zz > threshold[index]) {
      ordinal.value <- ordinal.value + 1
      index <- index + 1
    }
    ordinal.value
  })
  data.frame("y" = as.ordered(y), "x1" = x1, "x2" = x2)
}

# Simulate data
set.seed(922)
df5 <- simInteractionData(n = 2000)
pairs(df5)

# Save data sets
save(df1, file = "data/df1.rda", compress='xz')
save(df2, file = "data/df2.rda", compress='xz')
save(df3, file = "data/df3.rda", compress='xz')
save(df4, file = "data/df4.rda", compress='xz')
save(df5, file = "data/df5.rda", compress='xz')

# Data Generating Process for Partial Association --------------------------
library(multinomRob)

###########################
n <- 10000
beta1 <- 1
beta2 <- -1
alpha1 <- c(-Inf, -3, -2, 0, 2, 3, Inf)
alpha2 <- c(-Inf, -2, 0, 2, Inf)
X <- rnorm(n)

Z1 <- sapply(alpha1[2:6], function(alpha) alpha+beta1*X)
Z2 <- sapply(alpha2[2:4], function(alpha) alpha+beta2*X)

## function to obtain probability matrix of adjacent-categorical model
p_adj_cate <- function(Z){
  k<- ncol(Z)
  p1_pj<- p1_pj_inv<- Z
  ZZ<- 0
  for(j in 1:k){
    ZZ<- ZZ+Z[,j]
    p1_pj[,j]<- exp(ZZ)
    p1_pj_inv[,j]<- 1/p1_pj[,j]
  }
  p1<- 1/(rowSums(p1_pj_inv)+1)
  pj<- p1_pj
  for(j in 1:k){
    pj[,j]<- p1/p1_pj[,j]
  }
  cbind(1-rowSums(pj), pj)
}
############################

p1.mat <- p_adj_cate(Z1)
p2.mat <- p_adj_cate(Z2)
Y1 <- as.ordered(sapply(1:n, function(x) which(rmultinom(1, 1, p1.mat[x,])==1)))
Y2 <- as.ordered(sapply(1:n, function(x) which(rmultinom(1, 1, p2.mat[x,])==1)))
df_ParA <- list(data=data.frame(Y1=Y1, Y2=Y2, X=X, Z1=Z1, Z2=Z2),
                alpha1=alpha1, alpha2=alpha2, beta1=beta1, beta2=beta2)

save(df_ParA, file = "data/df_ParA.rda", compress='xz')

library(VGAM)
fit.adj1<- vglm(Y1~X, family=acat(reverse=TRUE, parallel=TRUE))
fit.adj2<- vglm(Y2~X, family=acat(reverse=TRUE, parallel=TRUE))
summary(fit.adj1)

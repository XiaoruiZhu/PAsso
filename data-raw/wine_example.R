library(ordinal)
?wine
data(wine, package = "ordinal") # load wine data set
summary(wine)
summary(wine$bottle)
wine.clm <- clm(rating ~ temp + contact, data = wine, link = "probit")

library(parasol)
test1 <- corr(responses = c("rating", "temp"),
              adjustments = c("contact", "judge"),
              data = wine, association = "partial")
test1
summary(test1)

test1_test <- corr.test(object = test1, boot_SE = 200, H0 = 0)
test1_test

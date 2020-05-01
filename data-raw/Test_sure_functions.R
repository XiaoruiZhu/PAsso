
# Test "residuals" generic function -----------------------------------------------
# Generate data from a quadratic probit model

set.seed(101)
n <- 2000
x <- runif(n, min = -3, max = 6)
z <- 10 + 3 * x - 1 * x^2 + rnorm(n)
y <- ifelse(z <= 0, yes = 0, no = 1)

# Scatterplot matrix
pairs(~ x + y + z)

# Setup for side-by-side plots
par(mfrow = c(1, 2))

# Misspecified mean structure
fm1 <- glm(y ~ x, family = binomial(link = "probit"))
scatter.smooth(x, y = residuals_glm(fm1),
               main = "Misspecified model",
               ylab = "Surrogate residual",
               lpars = list(lwd = 3, col = "red2"))
abline(h = 0, lty = 2, col = "blue2")

# Correctly specified mean structure
fm2 <- glm(y ~ x + I(x ^ 2), family = binomial(link = "probit"))
scatter.smooth(x, y = residuals_glm(fm2),
               main = "Correctly specified model",
               ylab = "Surrogate residual",
               lpars = list(lwd = 3, col = "red2"))
abline(h = 0, lty = 2, col = "blue2")


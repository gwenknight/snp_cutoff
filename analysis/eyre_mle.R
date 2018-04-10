set.seed(1001)

N <- 100

 x <- rnorm(N, mean = 3, sd = 2)

 mean(x)
 sd(x)



LL <- function(mu, sigma) {
     R = dnorm(x, mu, sigma)
     #
     -sum(log(R))
 }

library(stats4)

mle(LL, start = list(mu = 1, sigma=1))
mle(LL, start = list(mu = 1, sigma=1), method = "L-BFGS-B", lower = c(-Inf, 0),
      upper = c(Inf, Inf))
#### Code to replicate  (de Silva 2016 analysis)
library(stats4)

set.seed(1001)

### Likelihood from paper

f <- function(u){ ((mu*(t+2*u))^s) * exp(-u*(2*mu + 1 / N))  }

LL <- function(mu, N) {
  (exp(-mu*t)/N*factorial(s)) * integrate(f, 0, Inf, rel.tol = 1e-15)
}

# data
t = c(1,1,2,3,1,2,1) # distribution of times between samples
s = c(10,10,20,30,10,20,10) # SNP differences

t = 1
s = 10

# try to fit
mle(LL, start = list(mu = 1, N=10))






### example

N <- 100

 x <- rnorm(N, mean = 3, sd = 2)

 mean(x)
 sd(x)



LL <- function(mu, sigma) {
     R = dnorm(x, mu, sigma)
     #
     -sum(log(R))
 }



mle(LL, start = list(mu = 1, sigma=1))
mle(LL, start = list(mu = 1, sigma=1), method = "L-BFGS-B", lower = c(-Inf, 0),
      upper = c(Inf, Inf))

rm(list=ls())

library(dplyr)
library(psych)
library(AER)
library(gmm)

set.seed(1045)

obs <- 100001


# generate tfp as a stochastic AR(1) process


rho <- 0.5

shkvar <- 0.05

shk <- rnorm(obs, mean = 0, sd = sqrt(shkvar))

stat <- rnorm(obs, mean = 0, sd = sqrt(shkvar/(1-rho^2)))

a <- stat

for (i in 2:obs) {
  a[i] <- rho*a[i-1] + shk[i]
}


# checks
#
# dat <- data.frame(id = c(1:obs), a = a, laga = lag(a))
# 
# reg <- lm(a ~ laga, data = dat)
# 
# summary(reg)
# 
# var(a)
# var(stat)
# var(shk)


# generate log rental price of capital as a stochastic iid process


ufix <- log(0.1)

urndvar <- 0.1

u <- rnorm(obs, mean=ufix, sd=sqrt(urndvar))


# generate log wage as a stochastic MA(1) process


wfix <- log(1)

wrndvar <- 0.05

wrnd <- rnorm(obs, mean = 0, sd = sqrt(wrndvar))

phi <- 1

w <- wfix + wrnd + phi*lag(wrnd)


# checks
#
# wfix
# mean(w, na.rm = TRUE)
# var(wrnd)
# var(w, na.rm = TRUE)
# 
# wr <- wrnd
# wr1 <- lag(wrnd)
# 
# dat <- data.frame(id = c(1:obs), w = w, wr = wr, wr1 = wr1)
# 
# reg <- lm(w ~ wr + wr1, data = dat)
# 
# summary(reg)


# specify output elasticity parameters (DRS, sum to less than 1)

bx <- 0.6

bk <- 0.3


# generate log of capital (predetermined input)

rat <- bk/(1-bx)

fac <- 1/(1-rat)

pow <- 1/(1-bx)

tau <- bx^(bx*pow) - bx^pow

logetau <- log(tau) - bx*pow*wfix - phi*bx*pow*lag(wrnd) + ((bx*pow)^2)*wrndvar/2

k <- fac*(log(rat) + logetau - u + pow*rho*lag(a) + 0.5*shkvar*pow^2)


# generate log of revenue

r <- bx*pow*(log(bx) - w) + rat*k + pow*a


# generate revenue share for capital

sharek <- exp(u+k-r)


# arithmetic, geometric, harmonic means of revenue share for capital

mean(sharek, na.rm = TRUE)

geometric.mean(sharek, na.rm = TRUE)

harmonic.mean(sharek, na.rm = TRUE)


# generate log of labour (flexible input)

x <- pow*(log(bx) - w + bk*k + a)


# generate revenue share for labour

sharex <- exp(w+x-r)


# arithmetic, geometric, harmonic means of revenue share for labour

mean(sharex, na.rm = TRUE)

geometric.mean(sharex, na.rm = TRUE)

harmonic.mean(sharex, na.rm = TRUE)


# gross output and intermediate inputs - Leontief

thetav <- 1
thetam <- 1

theta <- thetam/thetav

i <- r - log(theta)


# introduce measurement error in log revenue

merrvar <- 0.1

merr <- rnorm(obs, mean = 0, sd = sqrt(merrvar))

r <- r + merr


# re-calculate revenue shares

sharek <- exp(u+k-r)

sharex <- exp(w+x-r)

# means of revenue share for capital

mean(sharek, na.rm = TRUE)

geometric.mean(sharek, na.rm = TRUE)

harmonic.mean(sharek, na.rm = TRUE)

# means of revenue share for labour

mean(sharex, na.rm = TRUE)

geometric.mean(sharex, na.rm = TRUE)

harmonic.mean(sharex, na.rm = TRUE)


# checks

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, a = a, w = w)

# infeasible ols using true tfp

olsinf <- lm(r ~ x + k + a, data = dat)

summary(olsinf)

# feasible ols omitting tfp

ols <- lm(r ~ x + k, data = dat)

summary(ols)




# ACF


# FIRST STAGE

# generate variables

xx <- x^2

kk <- k^2

ii <- i^2

xi <- x*i

ki <- k*i

xk <- x*k

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, i = i, xx = xx, kk = kk, ii = ii, xi = xi, ki = ki, xk = xk)

# ACF first stage regression (ols)

acf1 <- lm(r ~ x + k + i + xx + kk + ii + xk + xi + ki, data = dat, na.action = na.exclude)

summary(acf1)

# tests

linearHypothesis(acf1, c("i = 1"), test = "Chisq")

linearHypothesis(acf1, c("x = 0", "k = 0", "xx = 0", "kk = 0", "ii = 0", "xk = 0", "xi = 0", "ki = 0"), test = "Chisq")


# more parsimonious first stage regression

acf1a <- lm(r ~ i, data = dat, na.action = na.exclude)

summary(acf1a)

# test

linearHypothesis(acf1a, c("i = 1"), test = "Chisq")



# SECOND STAGE


# obtain fitted values from the first stage regression

rhat <- predict(acf1a)

# generate lagged variables

r1 <- lag(r)
r2 <- lag(r, n = 2L)

rhat1 <- lag(rhat)
rhat2 <- lag(rhat, n=2L)

x1 <- lag(x)
x2 <- lag(x, n = 2L)

k1 <- lag(k)
k2 <- lag(k, n = 2L)

dat <- data.frame(id = c(1:obs), r = r, rhat = rhat, x = x, k = k, r1 = r1, rhat1 = rhat1, x1 = x1, k1 = k1, r2 = r2, rhat2 = rhat2, x2 = x2, k2 = k2)


# unrestricted dynamic specification


ivacf2 <- ivreg(rhat ~ x + x1 + k + k1 + rhat1 | k + k1 + k2 + x1 + x2 + rhat1 + rhat2, data = dat)

summary(ivacf2)


# compare

ivlagxdyn <- ivreg(r ~ x + x1 + k + k1 + r1 | k + k1 + k2 + x1 + x2 + r2, data = dat)

summary(ivlagxdyn)



# restricted dynamic specification


varsacf2 <- cbind(rhat, rhat1, x, x1, k, k1, rhat2, x2, k2)

varsacf2 <- na.omit(varsacf2)

momsacf2 <- function(theta,x) {
  mk <- x[,5]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk1 <- x[,6]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk2 <- x[,9]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mx1 <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mx2 <- x[,8]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mrh1 <- x[,2]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mrh2 <- x[,7]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4]
  f <- cbind(mk, mk1, mk2, mx1, mx2, mrh1, mrh2, mc)
  return(f)
}

acf2 <- gmm(momsacf2, varsacf2, c(rho = 0, bx = 0.5, bk = 0.5, b0 = 0))

summary(acf2)


# compare

varsbb <- cbind(r, r1, x, x1, k, k1, r2, x2, k2)

varsbb <- na.omit(varsbb)

momsbb <- function(theta,x) {
  mk <- x[,5]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk1 <- x[,6]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk2 <- x[,9]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mx1 <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mx2 <- x[,8]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mr2 <- x[,7]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4]
  f <- cbind(mk, mk1, mk2, mx1, mx2, mr2, mc)
  return(f)
}

bb <- gmm(momsbb, varsbb, c(rho = 0, bx = 0.5, bk = 0.5, b0 = 0))

summary(bb)



# IMPOSE LEONTIEF SPECIFICATION

# generate lagged variables

i1 <- lag(i)
i2 <- lag(i, n = 2L)

dat <- data.frame(id = c(1:obs), i = i, x = x, k = k, i1 = i1, x1 = x1, k1 = k1, i2 = i2, x2 = x2, k2 = k2)


# unrestricted dynamic specification


ivacflf <- ivreg(i ~ x + x1 + k + k1 + i1 | k + k1 + k2 + x1 + x2 + i1 + i2, data = dat)

summary(ivacflf)



# MIS-SPECIFIED


# introduce measurement error in log intermediate inputs

merrintvar <- 0.2

merrint <- rnorm(obs, mean = 0, sd = sqrt(merrintvar))

i <- i + merrint


# FIRST STAGE

# generate variables

ii <- i^2

xi <- x*i

ki <- k*i

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, i = i, xx = xx, kk = kk, ii = ii, xi = xi, ki = ki, xk = xk, x1 = x1, k1 = k1, i1 = i1)

acf1mis <- lm(r ~ x + k + i + xx + kk + ii + xk + xi + ki, data = dat, na.action = na.exclude)

summary(acf1mis)

# tests

linearHypothesis(acf1mis, c("i = 1"), test = "Chisq")

linearHypothesis(acf1mis, c("x = 0", "k = 0", "xx = 0", "kk = 0", "ii = 0", "xk = 0", "xi = 0", "ki = 0"), test = "Chisq")

# alternative test

acf1misa <- lm(r ~ x + k + i + xx + kk + ii + xk + xi + ki + x1 + k1 + i1, data = dat, na.action = na.exclude)

summary(acf1misa)

linearHypothesis(acf1misa, c("x1 = 0", "k1 = 0", "i1 = 0"), test = "Chisq")


# SECOND STAGE


# obtain fitted values from the first stage regression

rhat <- predict(acf1mis)

rhat1 <- lag(rhat)
rhat2 <- lag(rhat, n=2L)

dat <- data.frame(id = c(1:obs), r = r, rhat = rhat, x = x, k = k, r1 = r1, rhat1 = rhat1, x1 = x1, k1 = k1, r2 = r2, rhat2 = rhat2, x2 = x2, k2 = k2)


# unrestricted dynamic specification


ivacf2mis <- ivreg(rhat ~ x + x1 + k + k1 + rhat1 | k + k1 + k2 + x1 + x2 + rhat1 + rhat2, data = dat)

summary(ivacf2mis)


# restricted dynamic specification


varsacf2mis <- cbind(rhat, rhat1, x, x1, k, k1, rhat2, x2, k2)

varsacf2mis <- na.omit(varsacf2mis)

acf2mis <- gmm(momsacf2, varsacf2mis, c(rho = 0, bx = 0.5, bk = 0.5, b0 = 0))

summary(acf2mis)

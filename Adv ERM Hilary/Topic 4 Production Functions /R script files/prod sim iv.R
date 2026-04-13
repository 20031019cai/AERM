rm(list=ls())

library(dplyr)
library(psych)
library(AER)
library(gmm)

set.seed(1045)

obs <- 100001


# generate tfp as a stochastic AR(1) process


rho <- 0

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




# 2SLS using log wage as the instrument for log labour
# consistent only if rho=0, wrndvar>0, phi =/=0 and/or urndvar>0

ivwage <- ivreg(r ~ x + k | k + w, data = dat)

summary(ivwage)


# replicate using gmm function

vars <- cbind(r, x, k, w)

vars <- na.omit(vars)

moments <- function(theta, x) {
  mk <- x[,3]*(x[,1] - theta[1]*x[,2] - theta[2]*x[,3] - theta[3])
  mw <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*x[,3] - theta[3])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*x[,3] - theta[3]
  f <- cbind(mk, mw, mc)
  return(f)
}

gmmwage <- gmm(moments, vars, c(bx = 0, bk = 0, b0 = 0))

summary(gmmwage)




# 2SLS using lagged labour as the instrument for log labour
# consistent only if rho=0, wrndvar>0, phi =/=0 and urndvar>0

x1 <- lag(x)

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, x1 = x1)

ivlagx <- ivreg(r ~ x + k | k + x1, data = dat)

summary(ivlagx)





# gmm allowing for AR(1) autocorrelation in tfp

# (1) using wage as instrument
# consistent only if wrndvar>0, phi =/=0 and/or urndvar>0


# generate lagged variables

r1 <- lag(r)
r2 <- lag(r, n = 2L)

# x1 <- lag(x)
x2 <- lag(x, n = 2L)

k1 <- lag(k)
k2 <- lag(k, n = 2L)

w1 <- lag(w)
w2 <- lag(w, n = 2L)


# unrestricted dynamic specification


dat <- data.frame(id = c(1:obs), r = r, r1 = r1, r2 = r2, x = x, x1 = x1, x2 = x2, k = k, k1 = k1, k2 = k2, w = w, w1 = w1, w2 = w2)

ivwagedyn <- ivreg(r ~ x + x1 + k + k1 + r1 | k + k1 + k2 + x1 + w + w1 + w2 + r2, data = dat)

summary(ivwagedyn)


# restricted dynamic specification


varsbbw <- cbind(r, r1, x, x1, k, k1, r2, x2, k2, w, w1, w2)

varsbbw <- na.omit(varsbbw)

momsbbw <- function(theta,x) {
  mk <- x[,5]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk1 <- x[,6]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mk2 <- x[,9]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mx1 <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mw <- x[,10]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mw1 <- x[,11]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mw2 <- x[,12]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mr2 <- x[,7]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3]*(x[,5] - theta[1]*x[,6]) - theta[4]
  f <- cbind(mk, mk1, mk2, mx1, mw, mw1, mw2, mr2, mc)
  return(f)
}

bbwage <- gmm(momsbbw, varsbbw, c(rho = 0, bx = 0.5, bk = 0.5, b0 = 0))

summary(bbwage)





# (2) using lagged labour as instrument
# consistent only if wrndvar>0, phi =/=0 and urndvar>0


# unrestricted dynamic specification


ivlagxdyn <- ivreg(r ~ x + x1 + k + k1 + r1 | k + k1 + k2 + x1 + x2 + r2, data = dat)

summary(ivlagxdyn)


# restricted dynamic specification


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




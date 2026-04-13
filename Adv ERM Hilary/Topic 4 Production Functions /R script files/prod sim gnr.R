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



# GNR


# log of revenue share for the flexible input

logsharex <- w + x - r


# geometric mean using lm

data <- data.frame(id = c(1:obs), logsharex = logsharex)

gnr1 <- lm(logsharex ~ 1, data = dat, na.action = na.exclude)

summary(gnr1)

logbxhat <- coef(summary(gnr1))["(Intercept)","Estimate"]

bxhat <- exp(logbxhat)


# geometric mean using nls

gnr1a <- nls(logsharex ~ log(bx), data = dat, start = list(bx = 0.5))

summary(gnr1a)


# alternative

# bxhat <- geometric.mean(sharex, na.rm = TRUE)



# SECOND STAGE using bxhat
# requires varurnd > 0


y <- r - bxhat*x

y1 <- lag(y)

y2 <- lag(y, n = 2L)

k1 <- lag(k)

k2 <- lag(k, n = 2L)

dat <- data.frame(id = c(1:obs), y = y, k = k, y1 = y1, k1 = k1, y2 = y2, k2 = k2)


# unrestricted dynamic specification


ivgnr2 <- ivreg(y ~ k + k1 + y1 | k + k1 + k2 + y2, data = dat)

summary(ivgnr2)


# compare

x1 <- lag(x)
x2 <- lag(x, n = 2L)

r1 <- lag(r)
r2 <- lag(r, n = 2L)

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, r1 = r1, x1 = x1, k1 = k1, r2 = r2, x2 = x2, k2 = k2)


ivlagxdyn <- ivreg(r ~ x + x1 + k + k1 + r1 | k + k1 + k2 + x1 + x2 + r2, data = dat)

summary(ivlagxdyn)



# restricted dynamic specification


varsgnr2 <- cbind(y, y1, k, k1, y2, k2)

varsgnr2 <- na.omit(varsgnr2)

momsgnr2 <- function(theta,x) {
  mk <- x[,3]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mk1 <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mk2 <- x[,6]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  my2 <- x[,5]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4])  - theta[3]
  f <- cbind(mk, mk1, mk2, my2, mc)
  return(f)
}

gnr2 <- gmm(momsgnr2, varsgnr2, c(rho = 0, bk = 0.5, b0 = 0))

summary(gnr2)


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



# SECOND STAGE using both bxhat and the first stage residuals
# requires varurnd > 0


# obtain residuals

merrhat <- residuals(gnr1)

# corrected dependent variable

ycorr <- r - bxhat*x + merrhat

ycorr1 <- lag(ycorr)

ycorr2 <- lag(ycorr, n = 2L)


dat <- data.frame(id = c(1:obs), ycorr = ycorr, k = k, ycorr1 = ycorr1, k1 = k1, ycorr2 = ycorr2, k2 = k2)


# unrestricted dynamic specification


ivgnr2a <- ivreg(ycorr ~ k + k1 + ycorr1 | k + k1 + k2 + ycorr1 + ycorr2, data = dat)

summary(ivgnr2a)


# compare

summary(ivgnr2)

summary(ivlagxdyn)


# restricted dynamic specification


varsgnr2a <- cbind(ycorr, ycorr1, k, k1, ycorr2, k2)

varsgnr2a <- na.omit(varsgnr2a)

momsgnr2a <- function(theta,x) {
  mk <- x[,3]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mk1 <- x[,4]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mk2 <- x[,6]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  myc1 <- x[,2]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  myc2 <- x[,5]*(x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4]) - theta[3])
  mc <- x[,1] - theta[1]*x[,2] - theta[2]*(x[,3] - theta[1]*x[,4])  - theta[3]
  f <- cbind(mk, mk1, mk2, myc1, myc2, mc)
  return(f)
}

gnr2a <- gmm(momsgnr2a, varsgnr2a, c(rho = 0, bk = 0.5, b0 = 0))

summary(gnr2a)


# compare

summary(gnr2)

summary(bb)
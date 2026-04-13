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

urndvar <- 0

u <- rnorm(obs, mean = ufix, sd = sqrt(urndvar))


# generate log wage as a stochastic MA(1) process


wfix <- log(1)

wrndvar <- 0.1

wrnd <- rnorm(obs, mean = 0, sd = sqrt(wrndvar))

phi <- 0

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




# OP FIRST STAGE
# requires rho>0, phi=0, wrndvar>0, urndvar=0


# generate additional explanatory variables

fk <- lead(k)

kk <- k^2
fkfk <- fk^2
kfk <- k*fk

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, fk = fk, kk = kk, fkfk = fkfk, kfk = kfk, a = a)

# op first stage regression (ols)

op1 <- lm(r ~ x + k + fk + kk + fkfk + kfk, data = dat, na.action = na.exclude)

summary(op1)

# inspect decision rule for capital

kdr <- lm(fk ~ k + a, data = dat)

summary(kdr)

# more parsimonious first stage regression (special case)

linearHypothesis(op1, c("kk = 0", "fkfk = 0", "kfk = 0"), test = "Chisq")

op1a <- lm(r ~ x + k + fk, data = dat, na.action = na.exclude)

summary(op1a)


# OP SECOND STAGE


# obtain fitted values from the first stage regression

rhat <- predict(op1a)

# rhat <- predict(op1)


# obtain estimated bx from the first stage regression

bxhat <- coef(summary(op1a))["x","Estimate"]

# bxhat <- coef(summary(op1))["x","Estimate"]

# bxhat <- summary(op1a)$coefficients[2,1]


# generate phihat

phihat <- rhat - bxhat*x

# generate dependent variable

y <- r - bxhat*x

# generate explanatory variables

phihat1 <- lag(phihat)

k1 <- lag(k)

dat <- data.frame(id = c(1:obs), y = y, k = k, phihat1 = phihat1, k1 = k1)

# op second stage regression (nlls)

op2 <- nls(y ~ bk*k + rho*(phihat1 - bk*k1), data = dat, start = list(bk = 0.5, rho = 1))

summary(op2)

# op2l <- lm(y ~ k + phihat1 + k1 - 1, data = dat)
# summary(op2l)


# VARIANTS

# translog first stage

xx <- x^2

xk <- x*k

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, xx = xx, kk = kk, xk = xk, fk = fk, fkfk = fkfk, kfk = kfk, a = a)

# op first stage regression (ols)

op1tr <- lm(r ~ x + k + xx + kk + xk + fk + fkfk + kfk, data = dat)

summary(op1tr)

linearHypothesis(op1tr, c("xx = 0", "kk = 0", "xk = 0"), test = "Chisq")

# quadratic tfp process

dat <- data.frame(id = c(1:obs), y = y, k = k, phihat1 = phihat1, k1 = k1)

# op second stage regression (nlls)

op2qu <- nls(y ~ bk*k + rho1*(phihat1 - bk*k1) + rho2*(phihat1 - bk*k1)^2, data = dat, start = list(bk = 0.5, rho1 = 1, rho2 = 0))

summary(op2qu)


# TESTING (add lagged inputs to first stage regression specification)

x1 <- lag(x)

dat <- data.frame(id = c(1:obs), r = r, x = x, k = k, fk = fk, kk = kk, fkfk = fkfk, kfk = kfk, x1 = x1, k1 = k1)

# op first stage regression (ols)

op1atest <- lm(r ~ x + k + fk + x1 + k1, data = dat)

summary(op1atest)

linearHypothesis(op1atest, c("x1 = 0", "k1 = 0"), test = "Chisq")

# op first stage regression (ols)

op1test <- lm(r ~ x + k + fk + kk + fkfk + kfk + x1 +k1, data = dat)

summary(op1test)

linearHypothesis(op1test, c("x1 = 0", "k1 = 0"), test = "Chisq")





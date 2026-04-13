# Classwork 3
# Generate endogenous Nash equilibrium price and market share data for nested logit model

#Clear and add library

install.packages("ivreg")
install.packages("tidyverse")
install.packages("nleqslv")

library("ivreg")
library(tidyverse)
library(nleqslv) # to solve non linear first order conditions

rm(list = ls())


# Set seed for reproducibility
set.seed(123)

# Number of products
J <- 10
# Number of markets
T <- 1000
#Number of observations
N <- T*J

# N x 1 object with market ID integers: (1,...,1,2,...,2,...)
t <- rep(1:T, each=J)

# Let products be assigned randomly to two alternative groups
groupid <- sample(c(1, 2), size = N, replace = TRUE)

# Give a unique numeric to each unique market-group interaction
code_tg <- as.numeric(interaction(groupid, t, drop = TRUE))

# observed quality for each observation: standard uniform
xj  <- 0.5*runif(N)
# marginal cost: increasing in quality 
mc <-  xj + 0.5*runif(N)
# unobserved quality for each observation: standard Normal
xi <- 0.5*rnorm(N)  

# Create a data frame with exogenous data
exog.data <- data.frame(t, groupid,  xj, mc, xi)


# Demand parameters
alpha <- 1  # price
beta1 <- -2   # const
beta2 <- 1   # xj
sigma <- 0.5   # sigma
theta <- c(alpha,beta1,beta2,sigma)

#Load in functions

#===============================================================================

source("mkt_share_nlogit.R")
source("foc_nlogit.R")

#===============================================================================


#Test share function on mkt 1

mkt_ind <- t == 1
price0_t <- 1.25*mc[ mkt_ind ] # startingguesses: price_pre_i
exog.data_t = exog.data[mkt_ind ,]
s_t <- mkt_share(price0_t, exog.data_t,  theta)
sj_t = s_t$sj
sj_g = s_t$sj_g
sum(s_t$sj)

solution <- nleqslv(price0_t, foc_nlogit, jac=NULL, exog.data_t, theta)

isconv_t <- solution$termcd #should be 1 (can fail even when code correct)
isconv_t

price_t <- solution$x
price_t

markups <- price_t-mc[ mkt_ind ] #should be >0
markups
mc[ mkt_ind ]
rm(mkt_ind,price0_t,exog.data_t,sj_t,solution,markups, isconv_t, price_t) # Remove checking variables


#=======================================================================
# Get endogenous data


pj_all = c()
sj_all = c()
sj_g_all = c()
s0_all = c()
isconv_sum = 0

for (i in 1:T) {

  mkt_ind <- t == i
  price0_t <- 1.25*mc[ mkt_ind ] # startingguesses: price_pre_i
  exog.data_t = exog.data[mkt_ind ,]
  solution <- nleqslv(price0_t, foc_nlogit, jac=NULL, exog.data_t, theta)
  isconv_t <- solution$termcd
  price_t <- solution$x
  s_t <- mkt_share(price_t, exog.data_t,  theta)
  
  sj_t = s_t$sj
  sj_g_t = s_t$sj_g
  
  s0_t <- 1-sum(sj_t)
  s0_t <- rep(s0_t,J)

  pj_all = c(pj_all, price_t)
  sj_all = c(sj_all, sj_t)
  sj_g_all = c(sj_g_all, sj_g_t)
  s0_all = c(s0_all, s0_t)
  isconv_sum = isconv_sum + isconv_t

}


# frac converge
print(isconv_sum/T)

#summarise

summary(pj_all)
summary(sj_all)
summary(sj_g_all)
summary(s0_all)

# Invert market share
delta = log(sj_all)-log(s0_all)

#gen ln sj_g  variable
ln_sj_g = log(sj_g_all)

#OLS
model_ols <- lm(delta~xj+ln_sj_g+pj_all)
summary(model_ols) 
#Note positive bias of parameters on each endog variable

#IV with cost data
# Run simple logit, instrumenting for prices not nesting term
cor(pj_all,mc)
cor(xi,mc)

model_iv <- ivreg(delta~xj+ln_sj_g|pj_all|mc)
summary(model_iv)
#Note poor price parameter

# Sum/count first column for rows with the same value in the second column
temp <- as.numeric(tapply(code_tg, code_tg, length))
sum_tg <- temp[code_tg]


cor(ln_sj_g,sum_tg)


model_iv <- ivreg(delta~xj|pj_all+ln_sj_g|mc+sum_tg)
summary(model_iv)
#Better.


#===============================================================================

#Elasticities

#cross 

a_hat = -model_iv$coefficients[2] #change sign to positive (as in notation)
sig_hat = model_iv$coefficients[3] #change sign to positive (as in notation)

mkt_ind <- t == 1
sj_t = sj_all[mkt_ind]
sj_g_t = sj_g_all[mkt_ind]
pj_t = pj_all[mkt_ind]
g_t = groupid[mkt_ind]

# J x J symmetric indicator: row and column belong to same group
gg_t = (t(matrix(g_t,nrow=10,ncol=10))==matrix(g_t,nrow=10,ncol=10))
gg_t <- apply(gg_t, c(1,2), as.numeric)

# Get J x J matrix of cross-elasticities (a cell is row's share wrt column's price)
# Matrix creates column-wise so transpose to get products j=1,...,J on columns:

sp_mkt <- a_hat*pj_t*sj_t # (J x 1) j is produce whose price changes
CrossElas1 = matrix(sp_mkt, nrow=J,ncol=J)
CrossElas1 = t(CrossElas1)

sg_t = sj_t/sj_g_t
thing = (1+ (1/sg_t)*sig_hat/(1-sig_hat))
thing = matrix(thing, nrow=J,ncol=J) # (J x 1) j is produce whose price change

CrossElas =  gg_t*CrossElas1*thing + (1-gg_t)*CrossElas1

#Get own-price elasticity vector and put it on diagonal:
ownElas <- -a_hat*pj_t*(1-sig_hat*sj_g_t-(1-sig_hat)*sj_t)/(1-sig_hat)
diag(CrossElas) <- ownElas

#Take a look (round to 3 dec. places)
round(CrossElas,3)
round(pj_t,3)

plot(pj_t,ownElas)
#Comment on flexibility relative to logit
#Comment on magnitude of own-price elasticities


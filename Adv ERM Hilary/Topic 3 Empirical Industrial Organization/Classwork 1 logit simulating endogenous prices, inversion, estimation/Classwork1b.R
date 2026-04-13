# Classwork 1b: Generating data from a logit oligopoly model

install.packages("ivreg")
install.packages("tidyverse")
install.packages("nleqslv")

#Clear and add libraries
rm(list = ls())


library("ivreg")
library(tidyverse)
library(nleqslv) # to solve non linear first order conditions


#LOAD IN FUNCTIONS (enter as a block)===========================================

source("mkt_share.R")
source("foc_logit.R")

# Set seed for reproducibility
set.seed(123)

# SET EXOGENOUS VARIABLES ======================================================

# Number of products,  markets, observations
J <- 10
T <- 1000
N <- T*J

# Let product ownership be assigned randomly to three alternative firms
f <- sample(c(1, 2, 3), size = N, replace = TRUE)

# N x 1 object with market ID integers: (1,...,1,2,...,2,...)
t <- rep(1:T, each=J)

# Give a unique numeric to each unique market-firm interaction
tf <- as.numeric(interaction(f, t, drop = TRUE))


# observed characteristics, unobserved characteristics and cost
xj  <- runif(N)
xi <- rnorm(N)  
mc <- xj + rnorm(N)

#a data frame with exogenous data
exog.data <- data.frame(t, f,  xj, mc, xi)


# Demand parameters
alpha <- 1  # price
beta1 <- -1   # const
beta2 <- 1   # xj
theta <- c(alpha,beta1,beta2)


#Checking: market share and Multi-product Nash pricing for market 1=============

mkt_ind <- t == 1
p0_t <- mc[ mkt_ind ] # guess
exog.data_t = exog.data[mkt_ind ,]
s_t<-mkt_share(p0_t,exog.data_t,  theta) 

sum(s_t) #should sum to <1

solution <- nleqslv(p0_t, foc_logit, jac=NULL, exog.data_t,theta )

isconv_t <- solution$termcd #should be 1 (can fail even when code correct)
isconv_t

price_t <- solution$x
price_t

markups <- price_t-mc[ mkt_ind ] #should be >0
markups

rm(mkt_ind,p0_t,exog.data_t,s_t,markups,isconv_t, price_t, solution) # Remove checking variables

# GET ENDOGENOUS DATA ALL MARKETS===============================================

#Empty objects for concat
pj_all = c()
sj_all = c()
s0_all = c()
isconv_sum = 0

for (k in 1:T) {
  
  mkt_ind <- t == k
  p0_t <- 1*mc[ mkt_ind ] # starting guess
  exog.data_t = exog.data[mkt_ind ,]
  
  Y <- nleqslv(p0_t, foc_logit, jac=NULL, exog.data_t,theta )
  
  isconv_t <- Y$termcd
  p1_t <- Y$x
  sj_t <- mkt_share(p1_t, exog.data_t,  theta)
  s0_t <- 1-sum(sj_t)
  s0_t <- rep(s0_t,J)
  
  pj_all = c(pj_all, p1_t)
  sj_all = c(sj_all, sj_t)
  s0_all = c(s0_all, s0_t)
  isconv_sum = isconv_sum + isconv_t
  
}

# frac converge
print(isconv_sum/T)

# INVERT AND ESTIMATE===========================================================

# Invert market share
delta = log(sj_all)-log(s0_all)


#OLS
model_ols <- lm(delta~xj+pj_all)
summary(model_ols) 

cor(pj_all,xi)
cor(pj_all,mc)
cor(xi,mc)

#IV: cost-shifters

model_iv <- ivreg(delta~xj|pj_all|mc)
summary(model_iv)


# Build markup-shifter (BLP) instruments: 
# sum of co-owned product characteristics


# Sum/count first column for rows with the same value in the second column
temp <- as.numeric(tapply(xj, tf, sum))
BLP_sumx <- temp[tf]

#IV: markup-shifters

model_iv <- ivreg(delta~xj|pj_all|BLP_sumx)
summary(model_iv)

cor(pj_all,BLP_sumx)

#2SLS with mc and BLP iv
model_iv <- ivreg(delta~xj|pj_all|mc+BLP_sumx)
summary(model_iv)


#===============================================================================

#MATRIX OF OWN AND CROSS-PRICE ELASTICITIES FOR MARKET t=1

a_hat = -model_iv$coefficients[2] #change sign to positive (as in notation)
mkt_ind <- t == 1
sj_t = sj_all[mkt_ind]
pj_t = pj_all[mkt_ind]


# Get JxJ matrix of cross-elasticities 
  # (a cell is elasticity of row's share wrt column's price)


sp_mkt <- a_hat*pj_t*sj_t # (J x 1) j is produce whose price changes
CrossElas = matrix(sp_mkt, nrow=J,ncol=J) # fills column-wise

# Matrix creates column-wise so transpose to get products j=1,...,J on columns:
CrossElas = t(CrossElas)

#Get own-price elasticity vector and put it on diagonal:
ownElas <- -a_hat*pj_t*(1-sj_t)
diag(CrossElas) <- ownElas

#Take a look (round to 3 dec. places)
round(CrossElas,3)

#Check that the high share products give high cross elasticities:
round(sj_t,3)



#Check relationship between price and ownElas
plot(pj_t, ownElas, main="Elast and Price", xlab="Price ", ylab="Elast ", pch=19)


#END============================================================================



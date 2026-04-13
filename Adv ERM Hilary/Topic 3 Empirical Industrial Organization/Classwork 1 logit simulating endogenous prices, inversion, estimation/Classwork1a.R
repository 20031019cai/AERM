rm(list = ls())
install.packages("ivreg")
install.packages("tidyverse")
install.packages("nleqslv")

library("ivreg")
library(tidyverse)
library(nleqslv) # to solve non linear first order conditions


# Set seed for reproducibility
set.seed(123)

# Number of products
J <- 10
# Number of markets
T <- 1000
#Number of observations
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

#Enter market share function before proceeding

#===============================================================================
#FUNCTIONS 
#===============================================================================

#Market share function==========================================================

mkt_share <- function(price_t, exog.data_t,  theta) {
  
  
  xj = exog.data_t$xj
  xi = exog.data_t$xi 
  
  u = theta[2] + theta[3]*xj - theta[1]*price_t + xi
  numerator <- exp(u)
  
  # This sums the columns with the same market
  denom  <- 1+sum(numerator)
  
  s <- numerator / denom
  return(s)
}

#===============================================================================

#Empty objects for concat
pj_all = c()
sj_all = c()
s0_all = c()

#Marginal cost pricing
for (k in 1:T) {
  
  mkt_ind <- t == k
  p1_t <- mc[ mkt_ind ] # price set to marginal cost
  exog.data_t = exog.data[mkt_ind ,]
  

  sj_i <- mkt_share(p1_t, exog.data_t,  theta)
  s0_i <- 1-sum(sj_i)
  s0_i <- rep(s0_i,J)
  
  pj_all = c(pj_all, p1_t)
  sj_all = c(sj_all, sj_i)
  s0_all = c(s0_all, s0_i)
  
  
}

#Remove redundant
rm(s0_i,sj_i,p1_t)

# Invert market share
delta = log(sj_all)-log(s0_all)


#OLS
model_ols <- lm(delta~xj+pj_all)
summary(model_ols) 

cor(pj_all,xi)
cor(pj_all,mc)
cor(xi,mc)

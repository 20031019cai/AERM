# Workscript for Classwork 2 (autos). 
#
# Plan:
# Get original BLP car data 
# OLS estimates: clustered and iid SEs
# Postestimation: own-price elasticities, markups and marginal costs; cross-price elasticity matrix
# Create BLP instruments from product characteristics
# IV estimates: clustered and iid SEs
# Postestimation: own-price elasticities, markups and marginal costs; cross-price elasticity matrix

install.packages("tidyverse")
install.packages("lmtest")
install.packages("sandwich")
install.packages("modelsummary")
install.packages("ivreg")

rm(list = ls())

library(tidyverse)
library(lmtest)
library(sandwich)
library(modelsummary)
library("ivreg")

# ==============================================================================
# Load functions

source("elasticities_logit.R")
source("make_BLP_IVs.R")
source("make_firm_shares.R")

# ==============================================================================
# Load original BLP data.
# Do this from File>Import Dataset > from text (base) >...
  
  

autoData <- read_csv("autoData.csv")

# I obtained this file from Genzkow and Shapiro's replication zip file. It's the original data file used in BLP.
# at https://web.stanford.edu/~gentzkow/research/blp_replication.zip  
 
#=============================================================================== 
# First have a look at the variables in spreadsheet form:
  
  View(autoData)

# Summarize the data. Compare with BLP's descriptive stats (Tables 1 & 2)

  summary(autoData)

  # In the datafile, the variable "logit_depvar" is ln(s_jt)-ln(s_0t) but it seems wrong so we
  
# Create delta_jt = ln (s_jt) - ln (s_0t)
  
  autoData$share_agg <- ave(autoData$share, autoData$year, FUN = sum)
  autoData$share_0 <- 1-autoData$share_agg
  
  autoData$delta = log(autoData$share) - log(autoData$share_0)

  mean(autoData$share_0)
  
  mean(autoData$delta)

# ==============================================================================
# Run the logit model. 
 

  model_ols <- lm(delta~hpwt+air+mpd+space+price, data=autoData)
  summary(model_ols)
 
# Cluster standard errors by car model variable "model_id"
  
  
# vcov tells it to do clustered SEs
# Cluster tells it to cluster on model 
# HC1 tells it to allow heteroskedasticity  
  
  model_ols_cl <- coeftest(model_ols,  vcov = vcovCL, type = "HC1", cluster = ~model_id)
  
# Compare SEs; compare with Table 3 in BLP
  
  modelsummary(list("OLS1"=model_ols, "OLS2"=model_ols_cl))

# Assign price parameter, price and share.
  
  alpha <- model_ols$coefficients[6]
  p = autoData$price
  s = autoData$share
  
# ==============================================================================  
# Calculate  own-price elasticities (put in value of alpha at start)
  
  elast_own_ols = alpha*p*(1-s)
  plot(density(elast_own_ols))
  summary(elast_own_ols)
  
# Fraction that are inelastic/implausible
  
  
  mean(elast_own_ols>-1)
  
# Look at how own-price elasticities vary with price and market share (flexibility)
  
  plot(p, elast_own_ols, main="Elast and Price", xlab="Price ", ylab="Elast ", pch=19)
  #plot(s, elast_own_ols, main="Elast and Share", xlab="Share ", ylab="Elast ", pch=19)
  
  
# Lerner indices, markups, and marginal costs assuming single-product Nash equilibrium
  
  lerner = -1/elast_own_ols
  plot(density(lerner))

  markup = p*lerner
  summary(markup)

  
  mc = p- markup
  summary(mc)
  plot(density(mc))
  
#Cross and own-price elasticity matrix in a specific market (year=1971)
#Column price, row product
  
#use the function "elasticities_logit.R" in the folder, putting in the price parameter for alpha
  Elast_ols = elasticities_logit(alpha)
  round(Elast_ols,4)
  
# ==============================================================================  
#Create instrumental variables for price using sums of characteristics
  
  #use the function "make_BLP_IVs.R" in the logit folder
  # this gives you 10 IVs in the data frame BLP_IVs
  #IV1-5 are "own" and IV6-10 are rival
  
  BLP_IVs <- make_BLP_IVs()
  View(BLP_IVs)
# ==============================================================================  
  
# Run simple logit, instrumenting for prices
  

  model_iv <- ivreg(autoData$delta~hpwt+air+mpd+space|price|
                      BLP_IVs$IV1+BLP_IVs$IV2+BLP_IVs$IV3+BLP_IVs$IV4+BLP_IVs$IV5+
                      BLP_IVs$IV6+BLP_IVs$IV7+BLP_IVs$IV8+BLP_IVs$IV9+BLP_IVs$IV10, 
                    data=autoData)
  # Look at parameters
  summary(model_iv)
  
# Get 
  
  model_iv_cl <- coeftest(model_iv,
                                        vcov = vcovCL,
                                        type = "HC1",
                                        cluster = ~autoData$model_id)
  

# Equivalent to BLP's Table 3:
  
  modelsummary(list("OLS"=model_ols_cl, "IV"=model_iv_cl))
  
# ==============================================================================
# Calculate (absolute) own-price elasticities (put in value of alpha at start)
  
  # for some reason, price parameter 
  model_iv$coefficients
  alpha2 <- model_iv$coefficients[2]
  elast_own_iv = alpha2*p*(1-s)
  summary(elast_own_iv)
  
# Fraction that are inelastic:
  

  mean(elast_own_iv>-1)
  
  
# Multi-product Nash equilibrium pricing: marginal costs etc.
# Need firm-level market shares.  
# Use make_firm_shares.R
  
  sf<-make_firm_shares()
  
  markup = -(1/alpha2)*(1/(1-sf))
  summary(markup)
  
  lerner = markup/p
  summary(lerner)

  mc = p - markup 
  summary(mc) 
  mean(mc<0)

  elast_own_iv = -1/lerner
  plot(p, elast_own_iv, main="Elast and Price", xlab="Price ", ylab="Elast ", pch=19)
  plot(s, elast_own_iv, main="Elast and Share", xlab="Share ", ylab="Elast ", pch=19)
  
  
    
# Own- and cross-price elasticity matrix (row share wrt column price )
# Note that it's inflexible: cross-elasticity same for both products in any column
# Products are
  # HD60071 min price,  (Honda, $3.4k, share=0.00014)
  # CVBISC71 med price, (Chevy, $7.6k,  share=0.00038)
  # LCMRKI71 max price, (Lexus, $21.8k,  share=0.00045)
  
  Elast_iv = elasticities_logit(alpha2)
  round(Elast_iv,4)  
  
  
  
  
#end
  
# ==============================================================================
  

  

  

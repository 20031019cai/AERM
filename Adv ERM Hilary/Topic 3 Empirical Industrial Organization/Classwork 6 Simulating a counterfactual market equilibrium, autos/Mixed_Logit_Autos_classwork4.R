# Classwork 4 (2026)
# This code estimates a random coefficients logit model on the autos data
# We use BLP instruments, calculate elasticities, and assess flexibility

# Note that in BLPestimator's vignette (linked below) they give an auto example
# which specifies the consumer heterogeneity differently from this code.
# Here, We do it as in the lecture notes' specification of BLP, in which income
# is interacted with price but no other random term. I do this differently from
# the vignette to give you another (hopefully useful) example of how to code BLP.
# This is actually closer to BLP (1995) than the example in the vignette

# Also different from the vignette: 
   # I load in the data rather than use data in library(BLPestimatoR)

# Vignette:
# (https://cran.r-project.org/web/packages/BLPestimatoR/vignettes/blp_intro.html)

#Install
install.packages("tidyverse")
install.packages("modelsummary")
install.packages("BLPestimatoR")
install.packages("nleqslv")

rm(list = ls())
set.seed(130) 

library(BLPestimatoR)
library(tidyverse)
library(modelsummary)
library(nleqslv) # to solve non linear first order conditions



#===============================================================================

#There are two auto datasets: 
  #1. autoData.csv (product chars)
  #2. dummies_cars.csv (ownership)


autoData <- read_csv("autoData.csv")
own_pre <- read_csv("dummies_cars.csv",)
#ownership_cars: JT x F ownership dummy (0,1) (F=26=#firms)

#Rename some variables in autoData:
colnames(autoData)[colnames(autoData) == "firm_id"] <- "firmid"
colnames(autoData)[colnames(autoData) == "year"] <- "cdid"
colnames(autoData)[colnames(autoData) == "car_id"] <- "id"
#cdid market identifier: t
#id observation identifier: jt (N=2217 unique values)

#This package likes the market number to start counting from 1:
autoData$cdid = autoData$cdid - 70 


# Add owner data matrix to productData

colnames(own_pre) <- paste0("company", 1:26) #column names
autoData <- cbind(autoData, own_pre)


#===============================================================================
#Make the BLP ownership instruments

nobs <- nrow(autoData)

#Make a data frame with the 5 characteristics
X <- data.frame(autoData$const, autoData$hpwt,
                autoData$air, autoData$mpg, autoData$space)

#Prepare two nobs x 5 matrices
sum_other <- matrix(NA, nobs, ncol(X))
sum_rival <- matrix(NA, nobs, ncol(X))

#Loop over observations adding characteristics within the same year:

for (i in 1:nobs) {
  other_ind <- autoData$firmid == autoData$firmid[i] &
    autoData$cdid == autoData$cdid[i] &
    autoData$id != autoData$id[i]
  rival_ind <- autoData$firmid != autoData$firmid[i] &
    autoData$cdid == autoData$cdid[i]
  sum_other[i, ] <- colSums(X[other_ind == 1, ])
  sum_rival[i, ] <- colSums(X[rival_ind == 1, ])
}

#Give them column names:
colnames(sum_other) <- paste0("IV", 1:5)
colnames(sum_rival) <- paste0("IV", 6:10)

#Add them to the data frame
autoData <- cbind(autoData, sum_other, sum_rival)

head(autoData)
#===============================================================================
# Get demographic draws: inverse income, for interaction w/ price a la BLP(1995)
# Like BLP we draw from log-normal distribution with the following parameters.

# Standard deviation of log incomes, assuming empirically given in BLP(1995):
sigma_v = 1.72

# Log income means for years 1971 - 1990:
incomeMeans = c(2.01156, 2.06526, 2.07843, 2.05775, 2.02915, 2.05346, 2.06745,
                2.09805, 2.10404, 2.07208, 2.06019, 2.06561, 2.07672, 2.10437, 
                2.12608, 2.16426, 2.18071, 2.18856, 2.21250, 2.18377)

# 500 draws per market, standard normal (can reduce if your code runs slow):
ndraws = 500
nu_0 <-matrix( rnorm(20*ndraws,mean=0,sd=1), nrow=20, ncol=ndraws) 

# Q: should you use same 500 draws in different markets? A: Don't have to, won't.
head(nu_0[,1:10]) #t x ndraws

y_means <-matrix( incomeMeans, nrow=20, ncol=ndraws) #t x ndraws of means 

# income draws are log normal; convert it from logs to levels 
y_it = exp(y_means + sigma_v*nu_0)

hist(log(y_it),500)

# We want the inverse a la BLP (1995)
inv_inc = 1/y_it

#take a look (out of curiosity)
hist(inv_inc[inv_inc<1],500)

# BLPEstimatoR needs us to label columns like this:
colnames(inv_inc) <- paste0("draw_", 1:ndraws)
inv_inc[1:5,1:5] #take a look to check

# BLPEstimatoR needs market identifier (1-20) in first column (& must be labelled such)
inv_inc <- cbind(c(1:20), inv_inc)
colnames(inv_inc)[1] <- "cdid"

# Take a look to check headings, markets, etc are as required:
inv_inc[1:5,1:5]


# BLPEstimatoR wants demographics as as a "list" in R:
demog_auto <- list(invinc = inv_inc)

#===============================================================================
# Set up the model using as.formula()

#Set up the model and assign to blps_model
  #The market share variable  is first and precedes a "~": share~. 
  #Then the explanatory variables are in four potentially-overlapping groups:

  #share ~ linear variables in utility |  exogenous variables in utility | 
  #         variables with random coefficients | further instruments

  # Note that instruments automatically include "exogenous variables in utility" 
  # and 'further instruments' refer to instruments extra to these. 

  # when a "0" is included to stop BLPestimatoR adding a constant (its default).
  #         ... not needed as we have put one in explicitly.

  #Note that price does not have a linear parameter in the BLP specification.

blps_model <- as.formula("share ~ 0 + const +  hpwt + air + mpd + space |
                          0 + const + hpwt + air + mpd + space |
                          0 + price + const + hpwt + air + mpd + space |
                          0 + IV1 + IV2 + IV3 + IV4 + IV5 + IV6 + IV7 + IV8 + IV9 + IV10")


#===============================================================================
# Create a data object, using BLP_data(). 

# Specify the product, market, cluster, productData, demographicdraws, etc.


car_data <- BLP_data(
  model = blps_model,
  market_identifier = "cdid",             #market: t
  product_identifier = "id",              #unit of observation: (jt)
               
  productData = autoData,
  additional_variables = paste0("company", 1:26), # we use this in merger (later)
  group_structure = "model_id",           #clustering variable
  blp_inner_tol = 1e-9,                   #tolerance of contraction
  blp_inner_maxit = 5000,                 #max iterations of contraction
  integration_method = "MLHS",            #Modified Latin Hypercube Sampling
  demographic_draws = demog_auto,
  integration_accuracy = ndraws,          #ndraws draws for integral
  integration_seed = 48
  
)


#===============================================================================
# Create theta_guesses: the starting values for the sigma/nonlinear parameters 
# This also specifies heterogeneity: NA means model will not estimate a parameter.


# We want a (K+1=6) x 2 matrix.
# Arranged so we get zeros where we don't want a sigma parameter
# Like BLP: I want unobserved heterogeneity on the 5 characteristics not price
# ... and I want demographics on price not any other chars:

#theta_guesses <- matrix(c(0,-2.8,3.1,1.2,-1.33,-1.78,-23.77,0,0,0,0,0), nrow=6, ncol=2)

theta_guesses <- matrix(c(0,-14.12, 5.32, -3.32, -3.17, -0.62,-25.3,0,0,0,0,0), nrow=6, ncol=2)
#These theta guesses get it to a lower GMM than the first set.

# Naming the rows and columns is needed.
rownames(theta_guesses) <- c("price", "const", "hpwt", "air", "mpd", "space")
colnames(theta_guesses) <- c("unobs_sd", "invinc")

# Replace the zeros with NAs: BLPestimator will interpret this as no parameter
theta_guesses[theta_guesses == 0] <- NA

#Take a look to make sure
theta_guesses


#===============================================================================
# Now we can run our first BLP estimator:
#  Need to feed it the data, theta, and the minimization routine details

car_est <- estimateBLP(
  blp_data = car_data,
  par_theta2 = theta_guesses,
  solver_method = "BFGS", solver_maxit = 1000, solver_reltol = 1e-10,
  standardError = "cluster",
  extremumCheck = FALSE, 
  printLevel = 1 # this outputs interesting information on the steps
)

# Take a look and compare with BLP(1995):

summary(car_est)

# In trials I found  I got numerical problems in my optimization
# So I had to try different starting values.
# It should be robust to small changes in starting values (s.t. them not causing numerical problems)
# I found that there were local minima and the lowest (hopefuly global) minimum required several attempts.

delta_est <- car_est$delta

## update mean utility in data ( always use update_BLP_data() to update data object to maintain consistent data )
delta_data <- data.frame(
  "id" = car_data$parameters$product_id,
  "cdid" = car_data$parameters$market_id,
  "delta" = delta_est
)
car_data_updated <- update_BLP_data(
  data_update = delta_data,
  blp_data = car_data
)

car_data <- car_data_updated

#===============================================================================
#Calculate consumer-draw level market shares, s_ij, as needed for elasticities

# Need to create a labeled nonlinear parameters matrix from the estimates

# When using demographics, the colnames must match the names of provided demographics 
# (as in demographic_draws) and "unobs_sd". Row names of par_theta2 must
# match random coefficients as specified in model. 



theta2_price <- car_est$theta_rc["invinc*price"]
theta2_price

theta2_all <- matrix(car_est$theta_rc)
theta2_all

theta2 <- matrix(NA, nrow = 6, ncol = 2)

rownames(theta2) <- c( "const", "hpwt", "air", "mpd","space", "price")
colnames(theta2) <- c("unobs_sd", "invinc")

theta2[1,1] = car_est$theta_rc["unobs_sd*const"]
theta2[2,1] = car_est$theta_rc["unobs_sd*hpwt"]
theta2[3,1] = car_est$theta_rc["unobs_sd*air"]
theta2[4,1] = car_est$theta_rc["unobs_sd*mpd"]
theta2[5,1] = car_est$theta_rc["unobs_sd*space"]
theta2[6,2] = car_est$theta_rc["invinc*price"]

theta2


shareObj <- getShareInfo(
  blp_data = car_data,
  par_theta2 = theta2,
  printLevel = 2
)

#have a look at s_ij
shareObj$sij[1:5,1:5]


#===============================================================================
#Get elasticiticies

# Calculate elasticities  (for market t=20 i.e. year 1990)
Elas_autos<-get_elasticities(
  blp_data = car_data,
  share_info = shareObj,
  theta_lin = 0,
  variable = "price",
  market = "20"
)

# Check own-price elasticities to see if implausible:
elast_own_BLP = diag(Elas_autos)
summary(elast_own_BLP)
mean(elast_own_BLP>-1)

#there are perhaps some inelastic in market 16 and 17 but only a few.


#===============================================================================

#Now look at cross-elasticity and own-elasticity patterns


# 5589 is a Yugo selling at $3,393
# 5478 is a Ford Fiesta at $4,834
# 5582 is a WV Passat selling at $11,301
# 5501 is a Lexus at $16,105
# 5428 is a BMW 3 at $25,558
# 5434 is a BMW 7 series at $37.490436
# 5447 is a Cadillac at  $40,589
# 5590 is a Porsche 911 selling at $44,759


#Calculate elasticities for selected products
# Order them from cheapest to most expensive:

Elas_autos_v2<-get_elasticities(
  blp_data = car_data,
  share_info = shareObj,
  theta_lin = 0,
  variable = "price",
  products = c("5589", "5478", "5582", "5501","5428", "5434", "5447", "5590"),
  market = "20"
)


#Label columns for ease of reading
colnames(Elas_autos_v2) <- c("Yugo", "Fiesta", "Passat", "Lexus", "BMW3", 
                             "BMW7", "Cadillac", "Porsche")

rownames(Elas_autos_v2) <- c("Yugo", "Fiesta", "Passat", "Lexus", "BMW3", 
                             "BMW7", "Cadillac", "Porsche")

# effect of change in column price on row share
round(Elas_autos_v2,4)



#===============================================================================

# Counterfactual merger

# Here (more or less) we follow steps at the bottom of:
# https://cran.r-project.org/web/packages/BLPestimatoR/vignettes/blp_intro.html

#===============================================================================


## Pre-Merger data

#company dummies: nobs x 26
own_pre <- as.matrix(car_data$data$additional_data[, paste0("company", 1:26)])


theta2_price <- car_est$theta_rc["invinc*price"]
theta2_all <- matrix(car_est$theta_rc)



#===============================================================================

## back out marginal costs (assuming multi-product Nash)

market_id <- car_data$parameters$market_id
nmkt <- length(unique(market_id))     # number of markets
markups <- numeric(length(market_id)) # number of obs

sh <- shareObj$shares
prices_pre <- car_data$data$X_rand[, "price"]


for (i in 1:nmkt) {
  mkt_ind <- market_id == i
  share_i <- sh[ mkt_ind ]
  price_pre_i <- prices_pre[ mkt_ind ]
  scalar_i <- matrix(1 / share_i) %*% matrix(price_pre_i, nrow = 1)
  elasticities_i <- get_elasticities(
    blp_data = car_data,
    share_info = shareObj,
    theta_lin = 0, #theta1_price,
    variable = "price",
    market = i,
    printLevel = 0
  )
  
  derivatives_i <- elasticities_i / scalar_i # partial derivatives of shares wrt price
  own_pre_i <- own_pre[ mkt_ind, ]
  own_prod_pre_i <- own_pre_i %*% t(own_pre_i) # if element (i,j) equals 1, that means that prod i and j are produced by same firm

  #uses solve(A) to get inverse of A (see slide p126 of lectures)
  markups[mkt_ind] <- c(-solve(t(derivatives_i) * own_prod_pre_i) %*% share_i)
}
marg_cost <- prices_pre - markups

#Look at markups/marginal costs implied by demand estimates & multi-product Nash
summary(markups)
summary(marg_cost)
summary(marg_cost<0)

#===============================================================================

# Merger between company 16 and 19 (i.e. GM and Chrysler)
prices_post <- numeric(2217)

#Inspect product ownership
colSums(own_pre[market_id == 1,])

#New ownership matrix

own_post <- cbind(
  own_pre[, 1:15],
  own_pre[, 16] + own_pre[, 19],
  own_pre[, 17:18],
  own_pre[, 20:26]
)

theta_lin=0
source("foc_bertrand_mkt.R")



#Pick a market from 1-20
i=1

    mkt_ind <- market_id == i
  own_post_i <- own_post[ mkt_ind, ]
  own_prod_post_i <- own_post_i %*% t(own_post_i)
  price_pre_i <- prices_pre[ mkt_ind ]
  
  solution <- nleqslv(
    x = price_pre_i, foc_bertrand_mkt, # startingguesses: price_pre_i
    own_prod = own_prod_post_i,
    blp_data = car_data,
    mkt = i,
    marg_cost = marg_cost,
    theta_rc = theta2
  )
  
#Solved post-merger prices
prices_post[ market_id == i ] <- solution$x


#Percentage change in prices from the merger 
pch = 100*(prices_post[ market_id == i ] - prices_pre[ market_id == i ])/prices_pre[ market_id == i ]

#Have a look at percentage price changes
summary(round(pch,2))


# This has same length as pch. Is TRUE if product j is merger insider (FALSE if not) 
insider <- own_pre[market_id==i, 16]==1 | own_pre[market_id==i, 19]==1

#Calculate mean effects for three groups and make a nice table:

mean_effects = c(mean(pch), mean(pch[insider==TRUE]), mean(pch[insider==FALSE]))
mean_effectsM <- matrix(mean_effects, ncol = 3)
colnames(mean_effectsM) <- paste0(c("all", "insider","outsider")) #column names

round(mean_effectsM,4)


#This is quite a dramatic merger.
#Average outsider prices increase.
#Prices fall for some outsider products. 
  #This is surprising-aren't prices strategic complements?
  #Possible reason: merging firm loses its most price-sensitive consumers who other firms compete for by reducing price.

#It seems common for the equation solver not to work (the chosen market and merger should work).
#One has to be patient as it can be slow.
#There are more reliable alternative ways to solve for equilibrum prices than the Newton method.
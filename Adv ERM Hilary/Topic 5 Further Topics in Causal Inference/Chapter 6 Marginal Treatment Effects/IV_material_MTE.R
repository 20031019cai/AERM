## Chapter 4 - MTE Model 
# This is meant as a simple example of estimation using MTE.
# 


library(dplyr)
library(ggplot2)
library(stats)

## ---- An example of representation for ATE, ATT, ATNT, MTE
#Generating 5000 observations random sample
N<- 5000
sigma1<- 1.15^2; sigma0<- 0.5^2; sigmaV<- 1;
alpha0<- 0.02; phi<- 0.2;
gamma0<- 0.2; gamma1<- 0.51; gamma2<- 0.1; gamma3 <- -0.25; rho <- 0.5


set.seed(180618)
eps<- rnorm(N)
Z1<- rnorm(N,-1,1.5)
Z2<- rnorm(N,1,1.5)
Z3<- rbernoulli(N,p=0.3)
U1<- sigma1*eps
U0<- rho*U1 + sigma0*rnorm(N)
V<- rnorm(N,0,sigmaV)

# Generate vector of ones (5000x1)
ones<- rep(1,5000)

# Generate Z=(1,Z1,Z2)
Z<- cbind(ones,Z1,Z2,Z3) 

# Generate Gamma Vector-> Gamma=(Gamma0,..,Gamma2)
Gamma <- rbind(gamma0,gamma1,gamma2,gamma3) 



# Generating Potential Outcomes
Y1<- alpha0 + phi + U1
Y0<- alpha0 + U0

# Get latent Variable I= Z%*%Gamma + V
# Note the selection on gains (Y_1 - Y_0)
I <-  Y1-Y0 + (Z%*%Gamma) + V 
U <- -(U1-U0 + V)
F_U <- ecdf(U)

# Get D= Indicator[I>=0]
D<- ifelse(I>=0, 1, 0)

# Generate Y= DY1 + (1-D)Y0
Y<- D*Y1 + (1-D)*Y0

# Save in a data frame
Df<- data.frame('Delta' = Y1 - Y0, 'Y' = Y,
                'Y1' = Y1, 'Y0' = Y0,
                'DecVar' = D, 'Latent' = I, 'UnObs' = V,
                'Z1' = Z1, 'Z2' = Z2, 'Z3' = Z3)


#Plot
(ATE<- Df %>% pull(Delta) %>% mean())
(ATT<- Df %>% filter(DecVar == 1) %>% pull(Delta) %>% mean())
(ATNT<- Df %>% filter(DecVar == 0) %>% pull(Delta) %>% mean())
(MTE <- Df %>% filter(abs(phi+(Z%*%Gamma)) <= 0.005) %>% pull(Delta) %>%
    mean())

##---- Check the support of the propensity score
model<- glm(DecVar~Z1+Z2,family = binomial(link="probit"),data=Df)
Df$y_pred = predict(model, Df, type="response")

ggplot(Df, aes(x = y_pred, y = after_stat(density))) +
  geom_density(alpha = 0.4, fill = 'azure') 



##---- Estimate the ATE using the package ivmte
library(ivmte)

mte_Estimate <-ivmte(
  data = Df,
  outcome = Y,
  target = "ate",
  m0 = ~ u, # defines m0 = E(Y_0 | X,u)
  m1 = ~ u, # defines m0 = E(Y_1 | X,u)
  propensity =DecVar ~ Z1 + Z2 + Z3,
  bootstraps = 100
)
summary(mte_Estimate)
# This is works well!

## ---- Represent the moment response function
# the Moment response functions are stored in the subvector mtr.coef

# Extract the coefficients
mte_Estimate$mtr.coef
m0<- function(u) mte_Estimate$mtr.coef[1] + u*mte_Estimate$mtr.coef[2]
m1<- function(u) mte_Estimate$mtr.coef[3] + u*mte_Estimate$mtr.coef[4]
u <- seq(0,1,0.05)
par(mfrow = c(1,2))

#Plot the moment response functions

# Plot m_0(u)
plot(u,m0(u),type="l", col="green", ylim=range(-0.5,1))
# Plot m_1(u)
lines(u,m1(u),type="l", col = 'dodgerblue')
legend(0,-0.2,c("m_0(u)","m_1(u)"), lwd=c(2,2), col=c("green","dodgerblue"))
abline(0,0)

# Plot gains

# True gains
# Using loess in stats package
r1 <- loess(Y1-Y0~F_U(U),span = 0.75,deg = 1)
d1 <- predict(r1, u, se = TRUE,na.omit)

# Plot m_1(u)-m0(u) against E(Y_1 - Y_0 | u)
plot(F_U(U),Y1-Y0, ylab = "m_1(u) - m_0(u)")
lines(u,m1(u) - m0(u),type="l", lwd = 5, col="firebrick", ylim=range(-0.5,1))
lines(u,d1$fit, lwd = 5, col="green")

# The average gains decrease with the cost of participation u. This was expected given our generated data.
# In fact, our estimate of the MTE matches reasonably well the relation between Y_1 - Y_0 and F_U(U) 

## - Estimation with the pakage localIV reveals a problem of this package.
#library(localIV)
# Estimate the ATE using the package mte
#mod <- mte(selection = DecVar ~ Z1 + Z2 , outcome = Y ~ DecVar, data = Df,method = "localIV")
#ate <- ace(mod, "ate")
#mte_vals <- mte_at(u = seq(0.05, 0.95, 0.1), model = mod)
# the results are wrong - The estimation of the MTE seems to be biased (see below)

#library(ggplot2)

#ggplot(mte_vals, aes(x = u, y = value)) +
#  geom_line(size = 1) +
#  xlab("Latent Resistance U") +
#  ylab("Estimates of MTE at Average values of X") +
#  theme_minimal(base_size = 14)
# The estimate MTE is almost always negative, whereas 



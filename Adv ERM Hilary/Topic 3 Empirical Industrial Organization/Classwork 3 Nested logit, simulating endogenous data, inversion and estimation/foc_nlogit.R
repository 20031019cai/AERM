foc_nlogit <- function(par, exog.data_t,  theta) {
  # argument par: candidate for unknown prices
  
  a = theta[1]
  b1 = theta[2]
  b2 = theta[3]
  sigma = theta[4]
  
  xj = exog.data_t$xj
  xi = exog.data_t$xi 
  mc = exog.data_t$mc
  firmid = exog.data_t$firmid
  
  s<-mkt_share(par,exog.data_t,  theta)
  
  sj<-s$sj     # j's share
  sj_g<-s$sj_g # j's share of group g
  
  #Single-product Nash; case of alpha=1:
  
  markup <- (1-sigma)/(1-sigma*sj_g - (1-sigma)*sj)
  
  differences <- par - mc - markup
  
  return(differences)
}

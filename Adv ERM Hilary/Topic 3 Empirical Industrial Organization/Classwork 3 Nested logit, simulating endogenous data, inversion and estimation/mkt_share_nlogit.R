mkt_share <- function(price_t, exog.data_t,  theta) {
  
  
  a = theta[1]
  b1 = theta[2]
  b2 = theta[3]
  sigma = theta[4]
  
  xj = exog.data_t$xj
  xi = exog.data_t$xi 
  mc = exog.data_t$mc
  groupid = exog.data_t$groupid
  
  u = b1 + b2*xj - a*price_t + xi
  numj <- exp(u/(1-sigma))
  
  #sum over products of the same group
  y <- match(groupid,unique(groupid)) #this generates a local groupid with no gaps
  denomg <- as.numeric(tapply(numj, y, sum)) #G x 1, G is 1 or 2
  denomj_g <- denomg[y]                        #J x 1
  
  sj_g <- numj / denomj_g  #share within group g
  
  numg <- exp( (1-sigma)*log(denomj_g) )
  
  sg <- numg / (1 + sum(exp((1-sigma)*log(denomg)) ))
  
  sj = sj_g*sg
  
  
  
  s = data.frame(sj,sj_g)
  
  return(s)
}

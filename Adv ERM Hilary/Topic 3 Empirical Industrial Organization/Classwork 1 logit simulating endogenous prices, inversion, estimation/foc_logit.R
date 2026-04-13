#System of Nash Prices Function=================================================


foc_logit <- function(price, exog.data_t,  theta) {
  
  s<-mkt_share(price,exog.data_t,  theta) 
  
  mc = exog.data_t$mc
  f = exog.data_t$f
  
  #Gen t-specific index (some t don't have all three firms)
  y <- match(f,unique(f)) 
  
  agg <- as.numeric(tapply(s, y, sum))
  sf <- agg[y]
  
  #markup <- 1/(1-s) #Single-product Nash
  markup <- 1/(1-sf) #Multi-product Nash
  
  differences <- price - (mc + markup)
  
  return(differences)
}

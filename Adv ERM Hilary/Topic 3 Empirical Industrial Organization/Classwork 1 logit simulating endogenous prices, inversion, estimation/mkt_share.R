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

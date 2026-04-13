elasticities_logit <- function(alpha)
{
  
  #alpha is the price parameter and should be negative 
  
  #this function uses the car data assigned as follows:
  #autoData <- read_csv("C:...original_data_GS.csv")
  
  #pick market 1: logical indicators
  mkt_ind <- autoData$year == 71
  
  s_mkt <- autoData$share[mkt_ind ]
  p_mkt <- autoData$price[mkt_ind ]
  
  #number of products J_mkt
  J_mkt <- length(s_mkt)
  
  
  # Get own- and cross- price elasticity matrix
  
  # Cross elasticity of j share wrt k price = alpha*s_k*p_k (simple logit model)
  # Elasticity matrix form: j's on rows, k's on cols
  
  
  # Multiply share and price, for each product k:
  sp_mkt <- -s_mkt*p_mkt
  
  # Get J x J matrix:
  Elas = matrix(sp_mkt, nrow=J_mkt,ncol=J_mkt)
  
  # Since matrix creates column-wise, transpose to get products on columns:
  Elas = t(Elas)
  
  #Get own-price elasticity vector and put it on diagonal:
  ownElas <- p_mkt*(1-s_mkt)
  diag(Elas) <- ownElas
  
  #Mult by abs value of price coeff:
  Elas = alpha*Elas
  
  # product_selector <- [logical]  
  # Return some products using product_selector (which is logical form):
  # summary(autoData$price[autoData$year==71])
  # HD60071 min price (Honda, $3.4k)
  # CVBISC71 med price (Chevy, $7.6k)
  # LCMRKI71 max price (Lexus, $21.8k)
  
  m1 = which(autoData$model_id == "HD60071")
  m2 = which(autoData$model_id == "CVBISC71")
  m3 = which(autoData$model_id == "LCMRKI71")
  
  Elas = Elas[c(m1,m2,m3), c(m1,m2,m3)]
  
  colnames(Elas) <- c("Honda", "Chevy", "Lexus") 
  rownames(Elas) <- c("Honda", "Chevy", "Lexus") 
  
  return(Elas)  
  
}

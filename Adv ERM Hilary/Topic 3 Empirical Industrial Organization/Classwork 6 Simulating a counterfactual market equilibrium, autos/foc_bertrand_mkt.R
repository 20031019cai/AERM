foc_bertrand_mkt <- function(par, own_prod, blp_data, mkt, marg_cost, theta_rc) {
  # argument par: candidate for post merger prices
  # arguments own_prod, blp_data, mkt, marg_cost,  theta_rc: see previous code blocks
  
  # post merger updates: update the BLP_data object for market i
  market_ind <- blp_data$parameters$market_id == mkt
  
  
  tmp <- data.frame(
    "id" = blp_data$parameters$product_id,
    "cdid" = blp_data$parameters$market_id,
    "price" = blp_data$data$X_rand[, "price"]
  )
  
  tmp$price[ market_ind ] <- par
  
  
  new_blp_data <- update_BLP_data(
    blp_data = blp_data,
    data_update = tmp
  )
  
  ShareObj <- getShareInfo(
    blp_data = new_blp_data,
    par_theta2 = theta_rc,
    printLevel = 0
  )
  
  implied_shares <- as.matrix(ShareObj$shares[market_ind])
  
  elasticities_post_mkt <- get_elasticities(
    blp_data = new_blp_data,
    share_info = ShareObj,
    theta_lin = 0,
    variable = "price",
    market = mkt,
    printLevel = 0
  )
  
  scalar_mkt <- matrix(1 / implied_shares) %*% matrix(par, nrow = 1)
  derivatives_mkt <- elasticities_post_mkt / scalar_mkt
  
  markups_post <- c(-solve(t(derivatives_mkt) * own_prod) %*% implied_shares)
  differences <- par - marg_cost[market_ind] - markups_post
  
  return(differences)
}

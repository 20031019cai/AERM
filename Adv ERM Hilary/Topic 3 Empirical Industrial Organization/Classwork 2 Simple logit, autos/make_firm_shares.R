make_firm_shares <- function(){


# construct instruments
nobs <- nrow(autoData)
Z <- autoData$share 

sum_samef <- matrix(NA, nobs, 1)

for (i in 1:nobs) {
  same_ind <- autoData$firm_id == autoData$firm_id[i] &
    autoData$year == autoData$year[i]
  
  sum_samef[i,1] <- sum(Z[same_ind == 1])
  
}

colnames(sum_samef) <- paste0("sharef", 1)

return(sum_samef)

}
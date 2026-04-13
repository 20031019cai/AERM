make_BLP_IVs <- function(){


#this function uses the car data assigned as follows:
#autoData <- read_csv("C:...original_data_GS.csv")  
  
# construct instruments
nobs <- nrow(autoData)
Z <- data.frame(
  autoData$const, autoData$hpwt,
  autoData$air, autoData$mpd, autoData$space
)

sum_other <- matrix(NA, nobs, ncol(Z))
sum_rival <- matrix(NA, nobs, ncol(Z))

for (i in 1:nobs) {
  other_ind <- autoData$firm_id == autoData$firm_id[i] &
    autoData$year == autoData$year[i] &
    autoData$model_id != autoData$model_id[i]
  rival_ind <- autoData$firm_id != autoData$firm_id[i] &
    autoData$year == autoData$year[i]
  
  
  sum_other[i, ] <- colSums(Z[other_ind == 1, ])
  sum_rival[i, ] <- colSums(Z[rival_ind == 1, ])
  
}

colnames(sum_other) <- paste0("IV", 1:5)
colnames(sum_rival) <- paste0("IV", 6:10)
BLP_IVs <- data.frame(sum_other, sum_rival)
return(BLP_IVs)

}
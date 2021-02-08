library(dplyr)
library(readxl)
library(stats4)

########Calculate negative log likelihood for 6 parameters
Simple_Lik_func <- function (paramlist)
{
  #paramlist <- c(0.023686729, 0.006109094, 0.136253986, 0.006796158, 0.099881378,  0.123039868)
  
  I_1<- paramlist[1];
  I_2<- paramlist[2];
  I_3<- paramlist[3];
  A_1<- paramlist[4];
  B_1<- paramlist[5];
  G_1<- paramlist[6];
  
  U_1 <- (I_1 + A_1*Deter$P1 + B_1*Deter$D1 + G_1*Deter$F1)
  U_2 <- (I_2 + A_1*Deter$P2 + B_1*Deter$D2 + G_1*Deter$F2)
  U_3 <- (I_3 + A_1*Deter$P3 + B_1*Deter$D3 + G_1*Deter$F3)
  U_4 <- (      A_1*Deter$P4 + B_1*Deter$D4 + G_1*Deter$F4)
  
  Tot_Denom <- (exp(U_1) +  exp(U_2) +  exp(U_3) + exp(U_4))
  
  p_1 <-exp(U_1)/Tot_Denom
  p_2 <-exp(U_2)/Tot_Denom
  p_3 <-exp(U_3)/Tot_Denom
  p_4 <-exp(U_4)/Tot_Denom
  
  prob <- (p_1 * Deter$BC1) + (p_2 * Deter$BC2) + (p_3 * Deter$BC3) + (p_4 * Deter$BC4)
  
  return(-sum(log(prob)))
}

#Create a Q matrix of 20x6 from Standard Normal
Q_matrix <- t(replicate(50, rnorm(n = 6, mean = 0,sd = 1)))

### Compute Continuous Heterogentiy
Continuous_Likelihood <- function(param_list)
{
  #param_list <- (1:27)/100
  meancnt <- 6
  mean_matrix <- matrix(param_list[1:meancnt], nrow = 1)
  
  ## Create Gamma upper triangular matrix
  gamma_utm <- matrix(0L, nrow = meancnt, ncol = meancnt)
  gamma_utm[upper.tri(gamma_utm, diag = TRUE)] <- tail(param_list, -meancnt)
  
  ### Update Variance 
  MVNQ <- Q_matrix %*% gamma_utm
  ### update Mean
  MVNQ <-  t(apply(MVNQ, 1, function(x){x+ mean_matrix}))
  lllist <- apply(MVNQ, 1, Simple_Lik_func)
  return(-mean(lllist)) 
}


##########otimrsltbckup <- otimrslt
startlist <- c(-0.2, 0.1, 1.09, 0.1, 0.7,0.1,-1.69, 0.1, 0.27, 0.1, 0.37, 0.1);
otimrslt <- optimx::optimr(par=startlist,fn = Continuous_Likelihood, 
                           method = "BFGS",
                           hessian = TRUE)
View(otimrslt)
##saveRDS(otimrslt, "Optimresult")
### Extract standard errors for parameters
### Select the 6x6 matrix from the 27x27 hessian matrix
### Invert the 6x6 Hessian matrix to get Variance Covariance matrix
###
coefficient_hessian <- (otimrslt$hessian[1:6, 1:6])
std_error_mat <- solve(coefficient_hessian)
std_error_vals <- diag(std_error_mat)

otimrslt$par[1:6]/ std_error_vals


testmat <- matrix(0L, nrow = 6, ncol =  6)
testmat[upper.tri(testmat, diag = TRUE)]<- c(
  -6.149811e-07,-1.147272e-06, -5.157877e-07 , 1.969842e-07,
  2.878527e-07, -1.788305e-08, 1.150546e-07,  1.020549e-07 , 9.616033e-07,
  4.238154e-07,  1.903003e-06, 2.754066e-06, -1.896241e-06 , 8.086099e-07,
 -1.555130e-06, -7.525633e-07,-2.103802e-06,  7.333393e-07 ,-2.485881e-06,
  1.627960e-06, -2.031646e-07) 
View(testmat)
#############
SavedOptimresult <- readRDS("Optimresult.rds")
View(SavedOptimresult)

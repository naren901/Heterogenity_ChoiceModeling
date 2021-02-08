library(dplyr)
library(stats4)
library(maxLik)
library(readxl)

## Load detergent data
detergent_file <- "Detergent Data set.xlsx"
Deter <- read_excel(detergent_file, col_names = TRUE)

## Load number of household data
nh_ds = read.table("hhobs.asc")
nhh  = nrow(nh_ds);
rperhh = 125;
nparam = 6;
# ## Generate (492*125) x 6 
rand_values <- matrix(0L, nrow = nhh*rperhh, ncol = 6)
rnromvalues  <- replicate(nhh*rperhh, rnorm(nparam, mean = 0,sd =1));
rand_matrix <- matrix(rnromvalues, nrow = nhh*rperhh, ncol = nparam, byrow =TRUE)
 
##Global for LL function
uniquehhids <- unique(Deter$Household_id)
startvalue <- c(-0.2, 0.1, 1.09, 0.1, 0.7,0.1,-1.69, 0.1, 0.27, 0.1, 0.37, 0.1);

### LL part begins
### For each household record, generate 125 probabilities 
currentrand <- rand_matrix[1:125,]
Compute_Neg_LL <- function(params)
{
  #params <- c(-0.2, 0.1, 1.09, 0.1, 0.7,0.1,-1.69, 0.1, 0.27, 0.1, 0.37, 0.1);
  hhmeanprob <- matrix(0L,nrow = nhh, ncol = 1)
  print (params)
 
  for (hh in 1:length(uniquehhids)) 
  {
    ### Pick next 125 random vectors
   # currentrand <- rand_matrix[(rperhh*(hh-1)+1):(rperhh*hh),] ###125x6
   # currentrand <- rand_matrix[1:125,] ###125x6
    hhincidents <- subset(Deter, Deter$Household_id == uniquehhids[hh])
    incidcnts <- nrow(hhincidents) ## No of transactions
    
    prob_hh <- matrix(1L, nrow = rperhh, ncol = 1)
  
    I_02 <- params[1] + params[2] * currentrand[,1] #Intercept
    I_03 <- params[3] + params[4] * currentrand[,2] #Intercept
    I_04 <- params[5] + params[6] * currentrand[,3] #Intercept
    P_01 <- params[7] + params[8] * currentrand[,4] #Price coef
    D_01 <- params[9] + params[10] * currentrand[,5] # Display coef
    F_01 <- params[11] + params[12] * currentrand[,6] ##Feature coef
  
    ### Calculate probabilites for each incident of a household
    for(hhi in  1:incidcnts)
    {
      hhrecord <-  hhincidents[hhi,]
    
      U_1 <-  exp(       hhrecord$P1*P_01 +hhrecord$D1*D_01 + hhrecord$F1*F_01 );
      U_2 <-  exp(I_02 + hhrecord$P2*P_01 +hhrecord$D2*D_01 + hhrecord$F2*F_01 );
      U_3 <-  exp(I_03 + hhrecord$P3*P_01 +hhrecord$D3*D_01 + hhrecord$F3*F_01 );
      U_4 <-  exp(I_04 + hhrecord$P4*P_01 +hhrecord$D4*D_01 + hhrecord$F4*F_01 );
    
      prob_local <- U_1*hhrecord$BC1 + U_2*hhrecord$BC2 +U_3*hhrecord$BC3 + U_4*hhrecord$BC4
      prob_local <- (prob_local/(U_1 + U_2 + U_3 + U_4))
      
      prob_hh[,1] <- prob_hh[,1] * prob_local
    }

    hhmeanprob[hh,] <- mean(prob_hh)
  }

      hhmeanprob[,1] <- log(hhmeanprob[,1]);
      llvalue <- sum(hhmeanprob[,1]);
      print(llvalue);
      return(-llvalue)
      #return( hhmeanprob[,1]);
}
 
optimresult <- optimx::optimr(par = startvalue, fn = Compute_Neg_LL,
                              method = "BFGS", hessian =  FALSE,
                              control=list(maximize = FALSE) )

View(optimresult)
saveRDS(optimresult, "KT_Implentation_Result_WO_Hessian.rds")

optimresult_hess <- optimx::optimr(par = startvalue, fn = Compute_Neg_LL,
                              method = "BFGS", hessian =  TRUE,
                              control=list(maximize = FALSE) )
saveRDS(optimresult_hess, "KT_Implentation_Result_Hessian.rds")


#maxlikresult <- maxLik(Compute_Neg_LL, start = startvalue, method = "BFGS")


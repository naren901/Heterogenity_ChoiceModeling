library(dplyr)

Null_model_loglikelihood <- function(Dataset)
{
 
    Prob1 <- mean(Dataset$BC1)
    Prob2 <- mean(Dataset$BC2)
    Prob3 <- mean(Dataset$BC3)
    Prob4 <- mean(Dataset$BC4)
    
    Lik <- Dataset$BC1 *Prob1 + Dataset$BC2*Prob2 + Dataset$BC3*Prob3 +Dataset$BC4 *Prob4
    View(Lik)
    
     neglikelihood <- 2*sum (log(Lik))
     
     return (neglikelihood)
}

########Calculate negative log likelihood for 6 parameters
Log_lik_func <- function (I_1, I_2, I_3, A_1, B_1, G_1)
{
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

### Negative log likelihood for Brand specific coefficients
Log_lik_func_xp <- function (I_1, I_2, I_3, A_1, A_2, A_3, A_4, B_1, G_1, d_1)
{
  U_1 <- (I_1 + A_1*Deter$P1 + B_1*Deter$D1 + G_1*Deter$F1+d_1*Deter$Household_id)
  U_2 <- (I_2 + A_2*Deter$P2 + B_1*Deter$D2 + G_1*Deter$F2+d_1*Deter$Household_id)
  U_3 <- (I_3 + A_3*Deter$P3 + B_1*Deter$D3 + G_1*Deter$F3+d_1*Deter$Household_id)
  U_4 <- (      A_4*Deter$P4 + B_1*Deter$D4 + G_1*Deter$F4+d_1*Deter$Household_id)
  
  Tot_Denom <- sum(exp(U_1), exp(U_2), exp(U_3), exp(U_4))
  p_1 <-exp(U_1)/Tot_Denom
  p_2 <-exp(U_2)/Tot_Denom
  p_3 <-exp(U_3)/Tot_Denom
  p_4 <-exp(U_4)/Tot_Denom
  
  prob <- (p_1 * Deter$BC1) + (p_2 * Deter$BC2) + (p_3 * Deter$BC3) + (p_4 * Deter$BC4)
  return(-sum(log(prob)))
}

#### Compute simple loyalty values using averaging
Compute_Simple_Loyalty <-function(DeterData)
{
  Loyalty_df <- DeterData[c("Household_id", "BC1","BC2","BC3",  "BC4" )]
  Loyalty_df$BC1_cs <- 0
  Loyalty_df$BC2_cs <- 0
  Loyalty_df$BC3_cs <- 0
  Loyalty_df$BC4_cs <- 0
  
  currentid <- Loyalty_df$Household_id[1]
  for(i in  (2:nrow(Loyalty_df)))
  {
    if (currentid == Loyalty_df[i,]$Household_id )
    {
      Loyalty_df[i, ]$BC1_cs <- Loyalty_df[i-1, ]$BC1_cs + Loyalty_df[i-1,]$BC1
      Loyalty_df[i, ]$BC2_cs <- Loyalty_df[i-1, ]$BC2_cs + Loyalty_df[i-1,]$BC2
      Loyalty_df[i, ]$BC3_cs <- Loyalty_df[i-1, ]$BC3_cs + Loyalty_df[i-1,]$BC3
      Loyalty_df[i, ]$BC4_cs <- Loyalty_df[i-1, ]$BC4_cs + Loyalty_df[i-1,]$BC4
      
    }else {
      currentid <- Loyalty_df[i,]$Household_id
      Loyalty_df[i, ]$BC1_cs <- 0
      Loyalty_df[i, ]$BC2_cs <- 0
      Loyalty_df[i, ]$BC3_cs <- 0
      Loyalty_df[i, ]$BC4_cs <- 0
    } 
  }
  Loyalty_df_temp <- Loyalty_df
  
  csTotal <- (Loyalty_df$BC1_cs + Loyalty_df$BC2_cs + Loyalty_df$BC3_cs + Loyalty_df$BC4_cs)
 
  Loyalty_df_temp$BC1_cs <- Loyalty_df_temp$BC1_cs/csTotal
  Loyalty_df_temp$BC2_cs <- Loyalty_df_temp$BC2_cs/csTotal
  Loyalty_df_temp$BC3_cs <- Loyalty_df_temp$BC3_cs/csTotal
  Loyalty_df_temp$BC4_cs <- Loyalty_df_temp$BC4_cs/csTotal
  
  Loyalty_df_temp$BC1_cs[is.nan(Loyalty_df_temp$BC1_cs)] <- 0
  Loyalty_df_temp$BC2_cs[is.nan(Loyalty_df_temp$BC2_cs)] <- 0
  Loyalty_df_temp$BC3_cs[is.nan(Loyalty_df_temp$BC3_cs)] <- 0
  Loyalty_df_temp$BC4_cs[is.nan(Loyalty_df_temp$BC4_cs)] <- 0
  
  DeterData$SL1 <- Loyalty_df_temp$BC1_cs
  DeterData$SL2 <- Loyalty_df_temp$BC2_cs
  DeterData$SL3 <- Loyalty_df_temp$BC3_cs
  DeterData$SL4 <- Loyalty_df_temp$BC4_cs
  
  return(DeterData)
}

##Negative Log likelihood with Loyalties
Log_lik_func_Simple_Loyalty <- function(I_1, I_2, I_3, A_1, B_1, G_1, L1, L2, L3, L4)
{
  U_1 <- (I_1 + A_1*Loyalty_data$P1 + B_1*Loyalty_data$D1 + G_1*Loyalty_data$F1 + L1*Loyalty_data$SL1)
  U_2 <- (I_2 + A_1*Loyalty_data$P2 + B_1*Loyalty_data$D2 + G_1*Loyalty_data$F2 + L2*Loyalty_data$SL2)
  U_3 <- (I_3 + A_1*Loyalty_data$P3 + B_1*Loyalty_data$D3 + G_1*Loyalty_data$F3 + L3*Loyalty_data$SL3)
  U_4 <- (      A_1*Loyalty_data$P4 + B_1*Loyalty_data$D4 + G_1*Loyalty_data$F4 + L4*Loyalty_data$SL4)
  
  Tot_Denom <- (exp(U_1) +  exp(U_2) +  exp(U_3) + exp(U_4))
  
  p_1 <-exp(U_1)/Tot_Denom
  p_2 <-exp(U_2)/Tot_Denom
  p_3 <-exp(U_3)/Tot_Denom
  p_4 <-exp(U_4)/Tot_Denom
  
  prob <- (p_1 * Deter$BC1) + (p_2 * Deter$BC2) + (p_3 * Deter$BC3) + (p_4 * Deter$BC4)
  return(-sum(log(prob)))
  
}
 
 
#Assigning initial value of (alpha) as 0.6 for the brand chosen, and (1-alpha)/3 for the other brands
## Generate Loyalties wherein the recent repeat purchases are weighted more
Computee_Weighted_Loyalty <- function(Detert_data, alpha =0.6)
{
  LD <- Detert_data[c("Household_id","Brand_chosen","BC1","BC2","BC3","BC4")]
  LD$Wl1 <- 0
  LD$Wl2 <- 0
  LD$Wl3 <- 0
  LD$Wl4 <- 0
  currentid <- LD$Household_id[1]
  
  ## Initialize loyalties for first instance of purchase 
  LD$Wl1[1] <- (alpha*LD$BC1[1] + ((1-alpha)*(1- LD$BC1[1]))/3)
  LD$Wl2[1] <- (alpha*LD$BC2[1] + ((1-alpha)*(1- LD$BC2[1]))/3)
  LD$Wl3[1] <- (alpha*LD$BC3[1] + ((1-alpha)*(1- LD$BC3[1]))/3)
  LD$Wl4[1] <- (alpha*LD$BC4[1] + ((1-alpha)*(1- LD$BC4[1]))/3)

  # Create cumulative purchase occasions until previous instance
  for(i in  (2:nrow(LD)))
  {
    if (currentid == LD[i,]$Household_id )
    {
      LD[i, ]$Wl1 <- LD[i-1, ]$Wl1*(alpha) + LD[i-1,]$BC1*(1 - alpha)
      LD[i, ]$Wl2 <- LD[i-1, ]$Wl2*(alpha) + LD[i-1,]$BC2*(1 - alpha)
      LD[i, ]$Wl3 <- LD[i-1, ]$Wl3*(alpha) + LD[i-1,]$BC3*(1 - alpha)
      LD[i, ]$Wl4 <- LD[i-1, ]$Wl4*(alpha) + LD[i-1,]$BC4*(1 - alpha)
    }
    else 
    {
      ### When consumer has changed, reset loyalties for first record again
      currentid <- LD[i,]$Household_id
      LD$Wl1[i] <- (alpha*LD$BC1[i] + ((1-alpha)*(1- LD$BC1[i]))/3)
      LD$Wl2[i] <- (alpha*LD$BC2[i] + ((1-alpha)*(1- LD$BC2[i]))/3)
      LD$Wl3[i] <- (alpha*LD$BC3[i] + ((1-alpha)*(1- LD$BC3[i]))/3)
      LD$Wl4[i] <- (alpha*LD$BC4[i] + ((1-alpha)*(1- LD$BC4[i]))/3)
    } 
  } 
  
  Detert_data$WL1 <- LD$Wl1
  Detert_data$WL2 <- LD$Wl2
  Detert_data$WL3 <- LD$Wl3
  Detert_data$WL4 <- LD$Wl4
  
  return(Detert_data)
}


##Negative Log likelihood with Loyalties
Neg_Log_lik_func_Weighted_Loyalty <- function(I_1, I_2, I_3, A_1, B_1, G_1, L1, L2, L3, L4)
{
  U_1 <- (I_1 + A_1*Weighted_Loyalty_data$P1 + B_1*Weighted_Loyalty_data$D1 + G_1*Weighted_Loyalty_data$F1 + L1*Weighted_Loyalty_data$WL1)
  U_2 <- (I_2 + A_1*Weighted_Loyalty_data$P2 + B_1*Weighted_Loyalty_data$D2 + G_1*Weighted_Loyalty_data$F2 + L2*Weighted_Loyalty_data$WL2)
  U_3 <- (I_3 + A_1*Weighted_Loyalty_data$P3 + B_1*Weighted_Loyalty_data$D3 + G_1*Weighted_Loyalty_data$F3 + L3*Weighted_Loyalty_data$WL3)
  U_4 <- (      A_1*Weighted_Loyalty_data$P4 + B_1*Weighted_Loyalty_data$D4 + G_1*Weighted_Loyalty_data$F4 + L4*Weighted_Loyalty_data$WL4)
  
  Tot_Denom <- (exp(U_1) +  exp(U_2) +  exp(U_3) + exp(U_4))
  
  p_1 <-exp(U_1)/Tot_Denom
  p_2 <-exp(U_2)/Tot_Denom
  p_3 <-exp(U_3)/Tot_Denom
  p_4 <-exp(U_4)/Tot_Denom
  
  prob <- (p_1 * Deter$BC1) + (p_2 * Deter$BC2) + (p_3 * Deter$BC3) + (p_4 * Deter$BC4)
  return(-sum(log(prob)))
}

##########################################################################################
### Segmentation likelihood
### This code is not completely modular to handle changes in parameter length dynamically
### Because the segmentation part is not made moduler yet.
###########################################################################################
Segmentation_likelihood <- function(params)
{
  nop <- 6;
  loopcnt <- floor(length(params)/nop)
  start <-1; end <- nop;
  uniqueids <- unique (Deter$Household_id)
  Probresults <- matrix(0L, nrow = length(uniqueids), ncol = loopcnt)
  
  for (sn in 1:loopcnt)
  {
    parvec <- params[start:end]
    print (parvec)
    U_1 <- (parvec[1]+  parvec[4]*Deter$P1 + parvec[5]*Deter$D1 + parvec[6]*Deter$F1)
    U_2 <- (parvec[2] + parvec[4]*Deter$P2 + parvec[5]*Deter$D2 + parvec[6]*Deter$F2)
    U_3 <- (parvec[3] + parvec[4]*Deter$P3 + parvec[5]*Deter$D3 + parvec[6]*Deter$F3)
    U_4 <- (            parvec[4]*Deter$P4 + parvec[5]*Deter$D4 + parvec[6]*Deter$F4)
 
    Tot_Denom <- (exp(U_1) +  exp(U_2) +  exp(U_3) + exp(U_4))
    
    p_1 <-exp(U_1)/Tot_Denom
    p_2 <-exp(U_2)/Tot_Denom
    p_3 <-exp(U_3)/Tot_Denom
    p_4 <-exp(U_4)/Tot_Denom
    
    prob <- (p_1 * Deter$BC1) + (p_2 * Deter$BC2) + (p_3 * Deter$BC3) + (p_4 * Deter$BC4)
    
    ###Compute Household wise likelihoods
    probrslt <- data.frame(Deter$Household_id, prob)
    colnames(probrslt) <- c("id", "prob")
    prodresult <- probrslt %>% aggregate.data.frame(by= list(probrslt$id), FUN = prod)
    Probresults[,sn] <- prodresult$prob
    
    ## Update the param indexes
    start <- start + end
    end <- nop + end
  }
  
  
  #### Handles only 2 segments
  segparam <- tail(x = params,n = 1)
  seg1=exp(segparam)/(1+exp(segparam)); 
  seg2=1-seg1;
  segvec <- c(seg1,seg2);
  
  ### Multiply each column by segment prob values
  Probresults_Re <- t( t(Probresults)*segvec)
  ProbresultsSum <- rowSums(Probresults_Re)
  #ProbresultsSum <-  Probresults[,1] * segvec[1] + Probresults[,2] * segvec[2]
  
  return(-sum (log (ProbresultsSum)))
}

######################

## Generate Loyalties wherein the recent repeat purchases are weighted more
NonLinear_loyalties_func <- function(Detergent_data, alpha)
{
  
  LD <- Detert_data[c("Household_id","Brand_chosen","BC1","BC2","BC3","BC4")]
  LD$Wl1 <- 0
  LD$Wl2 <- 0
  LD$Wl3 <- 0
  LD$Wl4 <- 0
  LD$Wl_NLP1 <- 0
  LD$Wl_NLP2 <- 0
  LD$Wl_NLP3 <- 0
  LD$Wl_NLP4 <- 0
  
  currentid <- LD$Household_id[1]
  
  ## Initialize loyalties for first instance of purchase 
  LD$Wl1[1] <- (alpha*LD$BC1[1] + ((1-alpha)*(1- LD$BC1[1]))/3)
  LD$Wl2[1] <- (alpha*LD$BC2[1] + ((1-alpha)*(1- LD$BC2[1]))/3)
  LD$Wl3[1] <- (alpha*LD$BC3[1] + ((1-alpha)*(1- LD$BC3[1]))/3)
  LD$Wl4[1] <- (alpha*LD$BC4[1] + ((1-alpha)*(1- LD$BC4[1]))/3)
  
  
  
  
  currentid <- Recent_Loyalty_dataframe$Household_id[1]
  current_brand_chosen <- Recent_Loyalty_dataframe$Brand_chosen[1]
  
  #Assigning initial value of (alpha) as 0.6 for the brand chosen, and (1-alpha)/3 for the other brands
  #alpha = 0.6
  if(current_brand_chosen == 1)
  {
    LD$Wl_NLP1<- 1
    LD$Wl_NLP2<- (-1/3)
    LD$Wl_NLP3<- (-1/3)
    LD$Wl_NLP4<- (-1/3)
  }
  else 
  {if(current_brand_chosen == 2)
  {
    LD$Wl_NLP1<- (-1/3)
    LD$Wl_NLP2<- 1
    LD$Wl_NLP3<- (-1/3)
    LD$Wl_NLP4<- (-1/3)
  }
    else
    {if(current_brand_chosen == 3)
    {
      LD$Wl_NLP1<- (-1/3)
      LD$Wl_NLP2<- (-1/3)
      LD$Wl_NLP3<- 1
      LD$Wl_NLP4<- (-1/3)
    }
      else
      {
        LD$Wl_NLP1<- (-1/3)
        LD$Wl_NLP2<- (-1/3)
        LD$Wl_NLP3<- (-1/3)
        LD$Wl_NLP4<- 1
      }
    }
  }  
  
  # Create cumulative purchase occasions until previous instance
  for(i in  (2:nrow(Recent_Loyalty_dataframe)))
  {
    if (currentid == Recent_Loyalty_dataframe[i,]$Household_id )
    {
      Recent_Loyalty_dataframe[i, ]$B1_loyalty <- Recent_Loyalty_dataframe[i-1, ]$B1_loyalty * (alpha) + Recent_Loyalty_dataframe[i-1,]$BC1 * (1 - alpha)
      Recent_Loyalty_dataframe[i, ]$B2_loyalty <- Recent_Loyalty_dataframe[i-1, ]$B2_loyalty * (alpha) + Recent_Loyalty_dataframe[i-1,]$BC2 * (1 - alpha)
      Recent_Loyalty_dataframe[i, ]$B3_loyalty <- Recent_Loyalty_dataframe[i-1, ]$B3_loyalty * (alpha) + Recent_Loyalty_dataframe[i-1,]$BC3 * (1 - alpha)
      Recent_Loyalty_dataframe[i, ]$B4_loyalty <- Recent_Loyalty_dataframe[i-1, ]$B4_loyalty * (alpha) + Recent_Loyalty_dataframe[i-1,]$BC4 * (1 - alpha)
      
      Recent_Loyalty_dataframe[i, ]$B1_loyalty_NLP <- Recent_Loyalty_dataframe[i-1, ]$B1_loyalty_NLP * (alpha) + Recent_Loyalty_dataframe[i-1, ]$B1_loyalty - Recent_Loyalty_dataframe[i-1,]$BC1 
      Recent_Loyalty_dataframe[i, ]$B2_loyalty_NLP <- Recent_Loyalty_dataframe[i-1, ]$B2_loyalty_NLP * (alpha) + Recent_Loyalty_dataframe[i-1, ]$B2_loyalty - Recent_Loyalty_dataframe[i-1,]$BC2 
      Recent_Loyalty_dataframe[i, ]$B3_loyalty_NLP <- Recent_Loyalty_dataframe[i-1, ]$B3_loyalty_NLP * (alpha) + Recent_Loyalty_dataframe[i-1, ]$B3_loyalty - Recent_Loyalty_dataframe[i-1,]$BC3 
      Recent_Loyalty_dataframe[i, ]$B4_loyalty_NLP <- Recent_Loyalty_dataframe[i-1, ]$B4_loyalty_NLP * (alpha) + Recent_Loyalty_dataframe[i-1, ]$B4_loyalty - Recent_Loyalty_dataframe[i-1,]$BC4 
    }
    else 
    {
      currentid <- Recent_Loyalty_dataframe[i,]$Household_id
      current_brand_chosen <- Recent_Loyalty_dataframe$Brand_chosen[i]
      
      #Assigning initial value of (alpha) as 0.6 for the brand chosen, and (1-alpha)/3 for the other brands
      if(current_brand_chosen == 1)
      {
        Recent_Loyalty_dataframe$B1_loyalty[i] <- alpha
        Recent_Loyalty_dataframe$B2_loyalty[i] <- (1-alpha)/3
        Recent_Loyalty_dataframe$B3_loyalty[i] <- (1-alpha)/3
        Recent_Loyalty_dataframe$B4_loyalty[i] <- (1-alpha)/3
        
        Recent_Loyalty_dataframe$B1_loyalty_NLP[i] <- 1
        Recent_Loyalty_dataframe$B2_loyalty_NLP[i] <- (-1/3)
        Recent_Loyalty_dataframe$B3_loyalty_NLP[i] <- (-1/3)
        Recent_Loyalty_dataframe$B4_loyalty_NLP[i] <- (-1/3)
      }
      else 
      {if(current_brand_chosen == 2)
      {
        Recent_Loyalty_dataframe$B2_loyalty[i] <- alpha
        Recent_Loyalty_dataframe$B1_loyalty[i] <- (1-alpha)/3
        Recent_Loyalty_dataframe$B3_loyalty[i] <- (1-alpha)/3
        Recent_Loyalty_dataframe$B4_loyalty[i] <- (1-alpha)/3
        
        Recent_Loyalty_dataframe$B2_loyalty_NLP[i] <- 1
        Recent_Loyalty_dataframe$B1_loyalty_NLP[i] <- (-1/3)
        Recent_Loyalty_dataframe$B3_loyalty_NLP[i] <- (-1/3)
        Recent_Loyalty_dataframe$B4_loyalty_NLP[i] <- (-1/3)
      }
        else
        {if(current_brand_chosen == 3)
        {
          Recent_Loyalty_dataframe$B3_loyalty[i] <- alpha
          Recent_Loyalty_dataframe$B1_loyalty[i] <- (1-alpha)/3
          Recent_Loyalty_dataframe$B2_loyalty[i] <- (1-alpha)/3
          Recent_Loyalty_dataframe$B4_loyalty[i] <- (1-alpha)/3
          
          Recent_Loyalty_dataframe$B3_loyalty_NLP[i] <- 1
          Recent_Loyalty_dataframe$B1_loyalty_NLP[i] <- (-1/3)
          Recent_Loyalty_dataframe$B2_loyalty_NLP[i] <- (-1/3)
          Recent_Loyalty_dataframe$B4_loyalty_NLP[i] <- (-1/3)
        }
          else
          {
            Recent_Loyalty_dataframe$B4_loyalty[i] <- alpha
            Recent_Loyalty_dataframe$B1_loyalty[i] <- (1-alpha)/3
            Recent_Loyalty_dataframe$B2_loyalty[i] <- (1-alpha)/3
            Recent_Loyalty_dataframe$B3_loyalty[i] <- (1-alpha)/3
            
            Recent_Loyalty_dataframe$B4_loyalty_NLP[i] <-  1
            Recent_Loyalty_dataframe$B1_loyalty_NLP[i] <-  (-1/3)
            Recent_Loyalty_dataframe$B2_loyalty_NLP[i] <-  (-1/3)
            Recent_Loyalty_dataframe$B3_loyalty_NLP[i] <-  (-1/3)
          }
        }
      }   #End of initializing loyalty values for alpha, as the household id is a new one
      
    } #End of check for household and flag 
    
  } #End of for loop
  
  return(Recent_Loyalty_dataframe)
}




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
Q_matrix <- t(replicate(20, rnorm(n = 6, mean = 0,sd = 1)))
 
 
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
  
  return(mean(lllist)) 
}
 

##########otimrsltbckup <- otimrslt
startlist <- (1:27)/100
otimrslt <- optimx::optimr(par=startlist,fn = Continuous_Likelihood, method = "BFGS", hessian = TRUE)
##saveRDS(otimrslt, "Optimresult")

coefficient_hessian <- (otimrslt$hessian[1:6, 1:6])
std_error_mat <- solve(coefficient_hessian)
std_error_vals <- diag(std_error_mat)

otimrslt$par[1:6]/ std_error_vals


 

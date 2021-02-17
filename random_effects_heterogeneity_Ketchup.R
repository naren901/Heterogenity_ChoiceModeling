# Random effects model

# Liquid detergents dataset
ld = read.table("detxmat6.asc")
# The data is in long format - we reshape it to a wide format
reshape_data = function(x){
  start_ind = (x[[1]][1]-1)*6+1
  end_ind = (x[[1]][1]-1)*6+6
  return(ld[start_ind:end_ind,])
}
res1= lapply(as.list(c(1:6849)), reshape_data)
reshape_data2 = function(x){
  return (c(x[4,],x[5,],x[6,]))
}
res2 = lapply(res1, reshape_data2)
res2 = lapply(res2, unlist)
data_in = do.call(rbind.data.frame,res2)
colnames(data_in) = c(paste(c("Wisk","Tide","Era","Surf"),rep(".price",4),sep = ""),
                      paste(c("Wisk","Tide","Era","Surf"),rep(".disp",4),sep = ""),
                      paste(c("Wisk","Tide","Era","Surf"),rep(".feat",4),sep = ""))

Choice = read.table("D://Personal//gauss_programs//prakhyasir_rep//detymat.asc")

# observations per household
nh = read.table("D://Personal//gauss_programs//prakhyasir_rep//hhobs.asc")

N = dim(data_in)[1]

# Total number of households
HH = dim(nh)[1]
# number of alternatives
nc = length(unique(Choice$V1))
# number of simulation draws per observation
R = 200

# Function to create draws of a certain mean and sd
fun_create_draws = function(n,mean,sd){
  return (rnorm(n,mean,sd))
}

repeat_fn = function(x,n){
  return(rep(x,n))
}
# Error draws -- these draws are fixed 
# Error draws are taken for each alternative specific shock to utility
# and for as many draws in each likelihood evaluation
e = mat.or.vec(N,nc*R)
f1 = function(sd)
{ return(unlist(lapply(as.list(nh$V1), fun_create_draws, mean=0,sd=sd)))
}
e[,1:(nc*R)] = sapply(as.list(rep(1,nc*R)),FUN = f1)
dim(e)

# Draws for household level parameters
# these draws are fixed and are not meant to be varying during optimization run
# draws for intercepts
sd=1
f2 = function(sd){
  return (unlist(mapply(rep,fun_create_draws(n = HH,mean = 0,sd = sd),as.list(nh$V1))))
}
b01draw = mat.or.vec(N,R); b02draw = mat.or.vec(N,R); b03draw = mat.or.vec(N,R)
b1draw = mat.or.vec(N,R); b2draw = mat.or.vec(N,R); b3draw = mat.or.vec(N,R)

b01draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)
b02draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)
b03draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)

# draws for coefficients on price, display and feature -  attribute effects
b1draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)
b2draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)
b3draw[,1:R] =sapply(as.list(rep(1,R)),FUN=f2)

# Choices are available as single categorical variable
# are converted to a set of dummy variables
D = fastDummies::dummy_columns(Choice$V1)
colnames(D) = c("Choice","D3","D4","D1","D2")

elementprod = function(x,B){
  return (x*B)
}

# Baseline model 
MNL_logit = function(betas){
  bprice = betas[1]
  bdisp = betas[2]
  bfeat = betas[3]
  b01 = 0
  b02 = betas[4]
  b03 = betas[5]
  b04 = betas[6]
  
  xb1 = exp(b01+bprice*data_in[,1]+bdisp*data_in[,5]+bfeat*data_in[,9])
  xb2 = exp(b02+bprice*data_in[,2]+bdisp*data_in[,6]+bfeat*data_in[,10])
  xb3 = exp(b03+bprice*data_in[,3]+bdisp*data_in[,7]+bfeat*data_in[,11])
  xb4 = exp(b04+bprice*data_in[,4]+bdisp*data_in[,8]+bfeat*data_in[,12])
  DNR = xb1+xb2+xb3+xb4
  prob1 = xb1/DNR; prob2 = xb2/DNR; prob3 = xb3/DNR; prob4 = xb4/DNR
  ll = prob1*D[,"D1"]+prob2*D[,"D2"]+prob3*D[,"D3"]+prob4*D[,"D4"]
  return(log(ll))
}
betas0 = c(-2,0.6,.44,2.20,1.8,0.65)
out_MNL = maxLik(MNL_logit,start = betas0,method = "BFGS")
out_MNL2 = maxLik(MNL_logit,start = betas0,method = "BFGS")


# Simulated likelihood function
# with logit errors
random_eff_ll <- function(params){
 
  b_price = params[1]; sb_price = params[2]
  b_display = params[3]; sb_display = params[4]
  b_feature = params[5]; sb_feature = params[6]
  b02 = params[7]; sb02 = params[8]
  b03 = params[9]; sb03 = params[10]
  b04 = params[11];sb04 = params[12]
  b01 = 0         
  
  # intercepts
  b02 = b02 + sb02*b01draw;   b03 = b03 + sb03*b02draw;   b04 = b04 + sb04*b03draw;
  # attributes
  b_price   = b_price   + sb_price*b1draw
  b_display = b_display + sb_display*b2draw
  b_feature = b_feature + sb_feature*b3draw
  
  xb1 = exp(b01+apply(b_price,2,elementprod,B=data_in[,1])+apply(b_display,2,elementprod,B=data_in[,5])+apply(b_feature,2,elementprod,B=data_in[,9]))
  xb2 = exp(b02+apply(b_price,2,elementprod,B=data_in[,2])+apply(b_display,2,elementprod,B=data_in[,6])+apply(b_feature,2,elementprod,B=data_in[,10]))
  xb3 = exp(b03+apply(b_price,2,elementprod,B=data_in[,3])+apply(b_display,2,elementprod,B=data_in[,7])+apply(b_feature,2,elementprod,B=data_in[,11]))
  xb4 = exp(b04+apply(b_price,2,elementprod,B=data_in[,4])+apply(b_display,2,elementprod,B=data_in[,8])+apply(b_feature,2,elementprod,B=data_in[,12]))
  DNR = xb1+xb2+xb3+xb4
  prob1 = apply(xb1/DNR,1,mean)
  prob2 = apply(xb2/DNR,1,mean)
  prob3 = apply(xb3/DNR,1,mean)
  prob4 = apply(xb4/DNR,1,mean)
  ll = prob1*D[,"D1"]+prob2*D[,"D2"]+prob3*D[,"D3"]+prob4*D[,"D4"]
  #print(sum(ll))
  return (log(ll))
}

params0 = c(-0.2,.1,1.09,.1,.7,.1,-1.69,.1,.27,.1,.37,.1) 
params1 = c(-1.69,.1,.4,.1,.65,.1,-0.2,.1,1.09,.1,.7,.1)
params2 = c(-1.59,.1,.17,.1,.37,.1,1.25,.1,.91,.1,.25,.1)

out1 = maxLik(random_eff_ll,start = params0,method = "BFGS")
out2 = maxLik(random_eff_ll,start = params1,method = "BFGS")
out3 = maxLik(random_eff_ll,start = params2,method = "BFGS")

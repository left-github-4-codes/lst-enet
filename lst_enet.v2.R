#Y. Zuo on June and July 10, 2023

rm(list=ls())
#library(glmnet)
library(lars)
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

##################################################
## Main body of the program
##################################################
#Simulating Data 
## for gaussian

lambda=c(0.1, 1, 10, 50, 200)  # could be longer and with smaller steps, tuning
alpha=seq(0, 1, length=5)       #5 points, could be more, tuning
set.seed(86)
n <- 100; p <- 25 # number of observations and variables, tuning as wanted
beta <- rep(0,p); beta[1:6] <- 1 # 24% nonzero coefficients
N_lars=100; alpha0=3; gam=1; n.folds=5 #allow LARS perform 100 times
sigma <- 0.5 # controls signal-to-noise ratio

eps <- 0.1 # contamination level, tuning as wanted
m <- ceiling(eps*n) # observations to be contaminated


#Simulation design I
x <- matrix(rnorm(n*(p-1), sigma),nrow=n)
x<-cbind(rep(1,n), x) #n by p matrix
e <- rnorm(n,0,1) # error terms

# #Simulation design II
# sigma0=matrix(0, nrow=p, ncol=p)
# sigma1=matrix(0, nrow=p1, ncol=p1); sigma2=matrix(0, nrow=p2, ncol=p2)
# rho1=0.95; rho2=0.05
# for (i in 1:p1){ for(j in 1:p1) { sigma1[i,j]=rho1^{abs(i-j)} } }
# for (i in 1:p2){ for(j in 1:p2) { sigma2[i,j]=rho2^{abs(i-j)} } }
# sigma0[1:p1, 1:p1]<-sigma1; sigma0[(p1+1):p,(p1+1):p]<-sigma2
# 
# x <- matrix(rmvnorm(n*p, mean=rep(0,p), sigma=sigma0), nrow=n, ncol=p )
# #x<-cbind(rep(1,n), x) #n by p matrix, this should be the design of the paper
# e <- rnorm(n,0,1);  # error terms

#Contamination scheme one:
eout <- e; eout[1:m] <- eout[1:m] + 10 # vertical outliers
yout <- c(x %*% beta + sigma * eout) # response
xout <- x; xout[1:m,] <- xout[1:m,] + 10 # bad leverage points
x=xout; y=yout

#contamination scheme two (based on the proof of Theorem 4.2)
#ind<-sample(n, m); Del=10^4;Kap=10^6
#yout <- c(x %*% beta + sigma * eout) # response
#xout[ind,]<-c(Del,rep(0,(p-1))); yout[ind]<-Del*Kap

#################################################################
lst_enet <- function(xout, yout, n.folds)
{#Train-test Split 70%-30% split
train_rows <- sample(1:n, .7*n, replace = F)

x.train <- xout[train_rows,]
y.train <- yout[train_rows]

x.test <- xout[-train_rows,]
y.test <- yout[-train_rows]

#using CV to get lambda and alpha
al_lam_vec<-CV_lst_enet(x.train, y.train, alpha, lambda, alpha0, N_lars, gam, n.folds)

alpha <- al_lam_vec[1]; lambda <- al_lam_vec[2]

#input the alpha and lambda to get lst_beta
beta_lst=get_lst_beta(x.train, y.train, lambda, alpha, alpha0, N_lars, gam)

rmse<- sqrt(mean((y.test-x.test%*%beta_lst)^2))
return(list(beta=beta_lst, RMSE=rmse))
}

#############################################################
### Cross-Validation (n.folds, typically 5 to 10)
#############################################################

CV_lst_enet1<- function(x, y, alpha, lambda, alpha0, N_lars, gam, n.folds)
{
  folds <- list() # flexible object for storing folds
  fold.size <- floor(nrow(x)/n.folds)
  remain <- 1:nrow(x) # all obs are in #sequence of original row numbers 
  
  for (i in 1:n.folds)
  {
    select <- sample(remain, fold.size, replace = FALSE)
    #randomly sample "fold_size" from the 'remaining observations' (via indices)
    
    folds[[i]] <- select # store indices
    #write a special statement for the last fold - if there are 'leftover points'
    
    if (i == n.folds){
      folds[[i]] <- remain # the last fold containing all left points
    }
    
    #update remaining indices to reflect what was taken out
    remain <- setdiff(remain, select)  #remain<-remain[-select]
  }
  
  results <- matrix(0, nrow=n.folds, ncol=3)
  for (i in 1:n.folds)
  { 
    #retrieve the ith fold
    indis <- folds[[i]] #unpack into a vector
    train_x <- x[-indis, ] #split into train and test sets
    train_y <- y[-indis]   #hold-one-cross-validation
    test_x <- x[indis,]; test_y <- y[indis] 
   
    al_len=length(alpha); lam_len=length(lambda)
    al_lam_rmse=matrix(0, ncol=3, nrow=lam_len)
    # a matrix to store alpha, lambda, and rmse
    
    #the following could be optimized by utilizing the mapply or outer
    for(ii in 1:lam_len)
    {
      RMSE=rep(0,al_len)
      for(jj in 1:al_len)
      { beta_new<-get_lst_beta(train_x, train_y, lambda[ii], 
                               alpha[jj], alpha0, N_lars, gam)
      
      wi=get_weight_wi(test_x, test_y, beta_new, alpha0) 
      indx_bt=which(wi==1); K =sum(wi)  
      pred<-test_x%*%beta_new; rn=as.vector(test_y-pred)
      mse=sum ((rn)^2*indx_bt)/K; RMSE[jj]<- sqrt(mse)
      }
      
      id=which.min(RMSE); min_rmse=RMSE[id]
      al_lam_rmse[ii,]<-c(alpha[id],lambda[ii],min_rmse)
    }  
    id_min=which.min(al_lam_rmse[,3])
    results[i,]<-al_lam_rmse[id_min,]
  }
  id_minn=which.min(results[,3])
  alpha=results[id_minn,1]; lambda=results[id_minn,2]
  return(c(alpha, lambda))
}
CV_lst_enet=cmpfun(CV_lst_enet1)

#################################################################

get_weight_wi1=function(x, y, beta, alpha)
{ # Input: x (n by p), y and beta (p by 1) vector, alpha a constant
  # Output: weight w_i=indicator(|r_i-med(r_i)|/mad(r_i)\leq alpha)  
  
  n=dim(x)[1]; p=dim(x)[2]
  wn=rn=index=rep(0, n); sd=rep(1, n)
  rn=y-x%*%matrix(beta, nrow=p)
  
  med=median(rn); madd=mad(rn)
  sd=abs(rn-med)/madd 
  #generalized deviations from median 
  index=which(sd<=alpha) # which returns indices meeting the condition
  #create a weight vector based on index
  wn[index] <- 1 # change entries to one based on the which returned indices
  return(c(wn))   # a weight vector (1 by n) with zero or one entries
}
get_weight_wi=cmpfun(get_weight_wi1)
####################################################################
#
get_lst_beta1=function(x, y, lambda, alpha, alpha0, N_lars, gam)
{
  # default value is 1 for alpha0, gamma, N_ls,
  # lambda and alpha are scalar.
  # alpha0 is used in the definition of LST
  # x is a n by (p-1) (or p) independent co-variate, y is a n by 1 response vector
  # N_lars is total the number to allow lars to run
  ################
  #initialization
  ###############
  n = dim(x)[1];  p = dim(x)[2]
  N1 = n * (n - 1) / 2           #total number of pick two distinct points
  N2 = choose(n, floor((n + 1) / 2)) #maximum number of K sequence in the paper
  N=min(N1, N2, N_lars)
  M_combi = matrix(0, nrow = 2, ncol = N1) #combination matrix
  M_bet = matrix(0, nrow = p, ncol = 2 * p)
  bet_0 = P_zero= matrix(0, nrow = p, ncol = 1)
  OBV_min=10 ^ 8 
  T_lars = CR= 0 
  bet_final = rep(0, p)
  delta=0.5
  lam<-lambda #for simplicity
  a<-alpha    #for simplicity
  
  M_combi = combn(1:n, 2)   # all combinations of two district indices
  ##################################################
  # big for loop considering all possible two points
  ##################################################
   for (k in 1:N1)  
  {
    i = M_combi[1, k]; j = M_combi[2, k]  #pick two points (i, j)
    xi = x[i,]; xj = x[j, ]; yi = y[i];  yj = y[j]
    
    dif = xi - xj       # 1 by p vector
    temp = which(dif != 0)
    
    l<-temp[1]      # just using the first index in temp1, reducing to one loop
    bet_0[l] <- (yi - yj) / (xi[l] - xj[l])  #no intercept any more (cf AA1 of ZZ22)
    
    ############################
    # create the 2p-beta matrix
    ###########################
    for (j in 1:p)
    { v0 = v00 = bet_0; v0[j, ] <- v0[j, ] + delta
      M_bet[, j] <- v0; v00[j, ] <- v00[j, ] - delta
      M_bet[, j + p] <- v00
    }
    ################
    #initialization
    ###############
    rn = I_bet = wi = rep(0, n); sd = rep(1, n) 
    for (kk in 1:(2*p)) #for loop
    { bt = M_bet[, kk]
      
      #calculate local Lars and update the beta_lars and OBV
      rn = y-x%*%matrix(bt, nrow= p); med = median(rn); madd = mad(rn)
      abdv = abs(rn - med);  sd = abdv / madd
      # 
      indx_bt=which(sd<=alpha0)
      wi[indx_bt]<-1 # change entries to one based on the which returned indices
      K =length(indx_bt); D_seqn <- sd[indx_bt]
      
      if (length(unique(D_seqn)) != K) {break} #could using differentiability for 
      else                                   #the computation otherwise
      { CR = CR + 1    #CR is the counter for the all D-sequences
        xx<-x[indx_bt,]; yy<-y[indx_bt]   #sub-dataset based on bt
      
        fit1=lars(xx, yy, type="lar", max.steps=900, use.Gram=F) #call lars to get a new beta
        fitted<-predict(fit1, newx=xx, s=10, type="fit")$fit
        bt<-predict(fit1, s=10, type="coef")$coef  #bt=bt_lars$coef; 
        T_lars = T_lars + 1
        
        indx_bt=get_weight_wi(xx, yy, bt, alpha0)
      
        rn=as.vector(yy-fitted); SS_lars =sum ((rn)^2*indx_bt)/sum(indx_bt)
        OBV_new=SS_lars + lam*( (1-a)*sum((abs(bt))^(gam))+a*sum(bt^2) )
        
        if (OBV_new < OBV_min){OBV_min = OBV_new; bet_final=bt}
      }# end of else
      
      if (CR > N) {break}
    } #end of 2nd loop on kk
    
    if (T_lars > N_lars) {break }
  }#big 1st loop on k end
  
  return(bet_final)
} # function end
get_lst_beta=cmpfun(get_lst_beta1)
#######################################################################
#

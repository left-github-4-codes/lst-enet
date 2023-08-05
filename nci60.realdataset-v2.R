#compare four-plus procedures by Y. Zuo on 07/17/23
#
rm(list=ls())
install.packages("enetLTS")
library(enetLTS)
install.packages("elasticnet")
library(elasticnet)
install.packages("robustHD")
library(robustHD)
library(glmnet)
install.packages("elasticnet")
install.packages("elasticnet")
library(elasticnet)
library(lars)
library(enetLTS)
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

##################################################
## Main body of the program
##################################################
#Simulating Data 
## for gaussian
library(robustHD)

data("nci60")
madvec<-apply(protein, 2, mad) #92 --> max of madvec
indexs<- which(madvec==median(madvec)) #(75, 85)
col_number=indexs[1]   #75
yout<-protein[,col_number]
correlations<-apply(gene, 2, corHuber,yout)
k1=100; k2=1000-k1
keep1<-partialOrder(abs(correlations), k1, decreasing=T)
keep2<-partialOrder(abs(correlations), k2, decreasing=F)
xout<-gene[,c(keep1,keep2)]

n=dim(xout)[1]; p=dim(xout)[2]
#set.seed(86)
# n <- 100; p <- 50 # number of observations and variables, tuning as wanted
# 
# #there are three choices for the true beta0
# 
# beta0 <- rep(0,p); 
# p1=floor(0.06*p); p2=p-p1
# beta0[1:p1] <- 3 # 24% nonzero coefficients
# #other choices of beta0
# #beta <- rep(1, p) # all are zeros
# #beta[6:10]<-2; beta[16:20]<-2 # zeros and 2 are alternating
# 
# sigma <- 0.5 # controls signal-to-noise ratio
 R=50; alph0=6
# 
# eps <- 0.2; alpha0=1 # contamination level, tuning as wanted
# m <- ceiling(eps*n) # observations to be contaminated
# 
p1.m=p2.m=p3.m=p4.m=p5.m=matrix(0, nrow=R, ncol=p)
rmse.lst_enet=rmse.lasso=rmse.lars=rmse.enetLTS=rmse.enet=rep(0, R)

# measure.lst_enet=measure.lasso=measure.lars=matrix(0,nrow=R, ncol=4)
# measure.enetLTS=measure.enet=matrix(0,nrow=R, ncol=4)
# #first component is L2-error, followed by TSDR, FSDR, the last one is the RMSE

# sigma0=matrix(0, nrow=p, ncol=p)
# sigma1=matrix(0, nrow=p1, ncol=p1); sigma2=matrix(0, nrow=p2, ncol=p2)
# rho1=0.95; rho2=0.05
# for (i in 1:p1){ for(j in 1:p1) { sigma1[i,j]=rho1^{abs(i-j)} } }
# for (i in 1:p2){ for(j in 1:p2) { sigma2[i,j]=rho2^{abs(i-j)} } }
# sigma0[1:p1, 1:p1]<-sigma1; sigma0[(p1+1):p,(p1+1):p]<-sigma2

#i=1
for(i in 1:R){
  
  # x <- matrix(rnorm(n*(p), sigma),nrow=n, ncol=p)
  # #x<-cbind(rep(1,n), x) #n by p matrix, this should be the design of the paper
  # e <- rnorm(n,0,1);  # error terms
  # #eout[1:m] <- eout[1:m] + 10 # vertical outliers
  
  
  # x <- matrix(rmvnorm(n*p, mean=rep(0,p), sigma=sigma0), nrow=n, ncol=p )
  # #x<-cbind(rep(1,n), x) #n by p matrix, this should be the design of the paper
  # e <- rnorm(n,0,1);  # error terms
  
  
  #Contamination scheme one:
  # eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
  # yout <- c(x %*% beta0 + sigma * eout) # response
  # xout <- x; xout[1:m,] <- xout[1:m,] + 20 # bad leverage points
  
  #Contamination scheme two:
  # eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
  # ind<-sample(n, m)
  # yout <- c(x %*% beta0 + sigma * eout) # response
  # xout <- x; xout[ind,] <- xout[ind,] + 20 # bad leverage points
  
  # #contamination scheme three (based on the proof of Theorem 4.2)
  # eout <- e; eout[1:m] <- eout[1:m] + 20 # vertical outliers
  # ind<-sample(n, m); Del=10^4;Kap=10^6
  # yout <- c(x %*% beta0 + sigma * eout) # response
  # xout <- x; xout[ind,]<-c(Del,rep(0,(p-1))); yout[ind]<-Del*Kap
  # 
  
  #################################################################
  #lst_enet <- function(xout, yout, n.folds)
  #{#Train-test Split 70%-30% split
  train_rows <- sample(1:n, .7*n, replace = F)
  
  x.train <- xout[train_rows,]
  y.train <- yout[train_rows]
  
  x.test <- xout[-train_rows,]
  y.test <- yout[-train_rows]
  
  #########################
  #lst_enet
  #########################
  lambda=seq(0.01, 10, by=20/20)  # could be longer and with smaller steps, tuning
  alpha=seq(0, 1, by=0.2)       #5 points, could be more, tuning
  N_lars=100;  gam=1; n.folds=5 #allow LARS perform 100 times
  
  al_lam_vec<-CV_lst_enet(x.train, y.train, alpha, lambda, alpha0, N_lars, gam, n.folds)
  
  al <- al_lam_vec[1]; lam<- al_lam_vec[2]
  
  #input the alpha and lambda to get lst_beta
  beta_lst=get_lst_beta(x.train, y.train, lam, al, alpha0, N_lars, gam)
  
   # wi=get_weight_wi(x.test, y.test, beta_lst, alpha0) 
   # indx_bt=which(wi==1); K =sum(wi)  
   # pred<-x.test%*%beta_lst; rn=as.vector(y.test-pred)
   # mse=sum ((rn)^2*indx_bt)/K; lst_enet.rmse<- sqrt(mse)
   # 
  print("t1")
  lst_enet.rmse<- sqrt(mean((y.test-x.test%*%beta_lst)^2))
  
  p1.m[i,]=beta_lst
  rmse.lst_enet[i]=lst_enet.rmse
  # L2.error_lst=sum((beta0-beta_lst)^2) 
  # TSDR_lst=TSDR(beta0, beta_lst); FSDR_lst=FSDR(beta0,beta_lst)
  # 
  
  # measure.lst_enet[i,1] <-L2.error_lst; measure.lst_enet[i,2] <-TSDR_lst
  # measure.lst_enet[i,3] <-FSDR_lst; measure.lst_enet[i,4]<-lst_enet.rmse
  
  #print(measure.lst_enet)
  #print(beta_lst); print(lst_enet.rmse)
  ###############################
  #LASSO
  ##############################
  # lasso<- cv.glmnet(x=x.train, y=y.train, type.measure="mse", alpha=1,
  #                   family="gaussian", nlambda=100)
  # 
  # lasso.predicted <- predict(lasso, lasso$lambda.1se, new=x.test)
  # 
  # beta_lasso<-coef(lasso, s="lambda.min")
  # 
  # lasso.rmse <- sqrt(mean((y.test - lasso.predicted)^2))
  
  fit2<- lars(x.train, y.train, type="lasso", max.steps=200, use.Gram=F)
  # 
  lasso.predicted <-predict(fit2, newx=x.test, s=10, type="fit")$fit
  # 
  beta_lasso<-predict(fit2,  s=10, type="coefficients")$coef
  # 
  # beta_lasso=bt_lasso$coef
  # 
  # fit3<-enet(x.train, y.train, lambda=0, max.steps=200)
  # 
  # beta_lasso<-predict.enet(fit3, s=0.99, type="coef", mode="fraction")$coef
  
  print("t2")
  lasso.rmse <- sqrt( mean((y.test-lasso.predicted)^2) )
  #lasso.rmse <- sqrt( mean((y.test-x.test%*% beta_lasso)^2) )
  p2.m[i,]=beta_lasso
  rmse.lasso[i]=lasso.rmse
  # 
  # L2.error_lasso=sum((beta0-beta_lasso)^2) 
  # TSDR_lasso=TSDR(beta0, beta_lasso); FSDR_lasso=FSDR(beta0,beta_lasso)
  # 
  # 
  # measure.lasso[i,1] <-L2.error_lasso; measure.lasso[i,2] <-TSDR_lasso
  # measure.lasso[i,3] <-FSDR_lasso; measure.lasso[i,4]<-lasso.rmse
  
  #print(measure.lasso[1,])
  #print(as.vector(beta_lasso[-1])); print(lasso.rmse)
  ###############################
  #LARS
  ###############################
  fit<- lars(x.train, y.train, type="lar", max.steps=100, use.Gram=F)
  
  lars.predicted <-predict(fit, newx=x.test, s=10, type="fit")$fit
  
  bt_lars<-predict(fit,  s=10, type="coefficients")
  
  beta_lars=bt_lars$coef
  
  print("t3")
  lars.rmse <- sqrt( mean((y.test-lars.predicted)^2) )
  
  p3.m[i,]=beta_lars
  rmse.lars[i]=lars.rmse
  # L2.error_lars=sum((beta0-beta_lars)^2) 
  # TSDR_lars=TSDR(beta0, beta_lars); FSDR_lars=FSDR(beta0,beta_lars)
  # 
  # 
  # measure.lars[i,1] <-L2.error_lars; measure.lars[i,2] <-TSDR_lars
  # measure.lars[i,3] <-FSDR_lars; measure.lars[i,4]<-lars.rmse
  # 
  #print(measure.lars)
  #print(as.vector(beta_lars)); print(lars.rmse)
  ###############################
  #enetLTS
  ###############################
  #install.packages("enetLTS")
  #library(enetLTS)
  # 
    fit1 <- enetLTS(x.train, y.train) #enetLTS can not handle xout with a 
   # #column that has the same element, e.g. all are 1 or c.
   # 
    beta_enetLTS<-as.vector(coef(fit1))
   # #this beta containing the beta0 which lead to its length p+1
   # 
   # #enetLTS.predicted <-predict(fit1, x.test, vers="raw", "response")
   # 
    print("t4")
  # # 
    enetLTS.rmse<- sqrt(mean((y.test-cbind(1,x.test)%*%beta_enetLTS)^2))
  # # 
    p4.m[i,]=beta_enetLTS[-1]
    rmse.enetLTS[i]=enetLTS.rmse
  #  L2.error_lts=sum((beta0-beta_enetLTS[-1])^2) 
  #  TSDR_lts=TSDR(beta0, beta_enetLTS[-1]); FSDR_lts=FSDR(beta0,beta_enetLTS[-1])
  # # 
  # # 
  #  measure.enetLTS[i,1] <-L2.error_lts; measure.enetLTS[i,2] <-TSDR_lts
  #  measure.enetLTS[i,3] <-FSDR_lts; measure.enetLTS[i,4]<-enetLTS.rmse
  
  #print(measure.enetLTS)
  #print(beta_enetLTS[-1]); print(enetLTS.rmse)
  ###############################
  #elasticnet
  ###############################
  #install.packages("elasticnet")
  #library(elasticnet)
  
  fit2<-enet(x.train, y.train, lambda=0.5, max.steps=1400)
  
  beta_enet<-predict.enet(fit2, s=150, type="coef", mode="step")$coef
  
  print("t5")
  
  enet.rmse<-sqrt(mean((y.test-x.test%*%beta_enet)^2))
  
  p5.m[i,]=beta_enet
  rmse.enet[i]=enet.rmse
  
  # L2.error_enet=sum((beta0-beta_enet)^2) 
  # TSDR_enet=TSDR(beta0, beta_enet); FSDR_enet=FSDR(beta0,beta_enet)
  # 
  # 
  # measure.enet[i,1] <-L2.error_enet; measure.enet[i,2] <-TSDR_enet
  # measure.enet[i,3] <-FSDR_enet; measure.enet[i,4]<-enet.rmse
  # 
  #print(measure.enet[1,])
  #print(beta_enet); print(enet.rmse)
  print("i"); print(i)
}

#print("contamination level epsilon"); print(ep)
m.lst=rmse.lst_enet
m.lasso=rmse.lasso
m.lars=rmse.lars
m.enetLTS=rmse.enetLTS
m.enet=rmse.enet

p1.mean=colMeans(p1.m)
p2.mean=colMeans(p2.m)
p3.mean=colMeans(p3.m)
p4.mean=colMeans(p4.m)
p5.mean=colMeans(p5.m)

p1.dev=p1.m-matrix(p1.mean, byrow=T, nrow=R, ncol=p)
p2.dev=p2.m-matrix(p2.mean, byrow=T, nrow=R, ncol=p)
p3.dev=p3.m-matrix(p3.mean, byrow=T, nrow=R, ncol=p)
p4.dev=p4.m-matrix(p4.mean, byrow=T, nrow=R, ncol=p)
p5.dev=p5.m-matrix(p5.mean, byrow=T, nrow=R, ncol=p)

p1.emse=sum(p1.dev*p1.dev)/R
p2.emse=sum(p2.dev*p2.dev)/R
p3.emse=sum(p3.dev*p3.dev)/R
p4.emse=sum(p4.dev*p4.dev)/R
p5.emse=sum(p2.dev*p2.dev)/R

dev.new(width=5, height=5, unit="cm")

dev.new(width=10, height=5, unit="in")

dev.new(width=100, height=50, unit="px")

par(mfrow = c(1,2))
boxplot(cbind(p1.emse, p2.emse, p3.emse, p4.emse, p5.emse), 
          names=c("P1", "P2","P3", "P4", "P5"),
          xlab="nci60 (n=59, p=1000)",
          ylab="EMSE",
          border="green"
  )

boxplot(cbind(m.lst, m.lasso, m.lars, m.enetLTS,  m.enet), 
        names=c("P1", "P2","P3", "p4", "P5"),
        xlab="nci60 (n=59, p=1000)",
        ylab="RMSE",
        border="blue"
)

 # boxplot(cbind(m.lst[,2], m.lasso[,2], m.lars[,2],   m.enet[,2]), 
 #         names=c("P1", "P2","P3",  "P5"),
 #         xlab="20% contamination",
 #         ylab="TSDR",
 #         border="red"
 # )
 # boxplot(cbind(m.lst[,3], m.lasso[,3], m.lars[,3],   m.enet[,3]), 
 #         names=c("P1", "P2","P3", "P5"),
 #         xlab="20% contamination",
 #         ylab="FSDR",
 #        border="orange"
 # )




#####################################################

TSDR<-function(beta_0, beta_hat)
{id=which(beta_0==0); den=length(id); neu=length(which(beta_hat[id]==0))
if (den>0){ return(neu/den)} else{return("no zero component"); break}
}

FSDR<- function(beta_0, beta_hat)
{ id=which(beta_0!=0); den=length(id); neu=length(which(beta_hat[id]==0))
return(neu/den)
}

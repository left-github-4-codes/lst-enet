#readme.first.txt by Y.Zuo on 08/04/23

(i) Install all needed packages and load those packages by library function

install.packages("enetLTS")
install.packages("elasticnet")
install.packages("robustHD")
install.packages("lars")
install.packages("glmnet")

library(lars)
library(enetLTS)
library(mvtnorm)
library(robustHD)
library(glmnet)
library(lars)
library(mvtnorm)
library(robustbase)
library(compiler)
enableJIT(1)

(ii) run functions (they are included in lst_enet.v2.R)
   get_weight_wi1
   get_lst_beta1
   CV_lst_enet1
   TSDR (at the end of compare.R)
   FSDR (at the end of compare.R)
   
(iii) now you can run  any of following 
   compare.R
   compare1.R
   comparev2.R
   comparev3.R
to get all the boxplots in Figures 2-7

(iv) run 
     nci60.realdataset.R, or
     nci60.realdataset-v2.R
to get Figure 8 in Example 7.4

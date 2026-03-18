# rm(list=ls())
library(MASS)
library(glmnet)
library(HDeconometrics)
library(tictoc)

n = 100
p = 300
cov = matrix(numeric(p*p),ncol=p)
for(i in 1:p){
  for(j in 1:p){
    cov[i,j] = 0.5^abs(i-j)
  }
}
num = 50
seed = sample(1:1e6, num)
# seed = 1
res_GCV_type1 = matrix(numeric(5*num),ncol=5); colnames(res_GCV_type1) = c("Ave", "TN", "TP", "MSE", "AUC")
res_AIC_type1 = matrix(numeric(5*num),ncol=5); colnames(res_AIC_type1) = c("Ave", "TN", "TP", "MSE", "AUC")
res_BIC_type1 = matrix(numeric(5*num),ncol=5); colnames(res_AIC_type1) = c("Ave", "TN", "TP", "MSE", "AUC")
res_AICC_type1= matrix(numeric(5*num),ncol=5); colnames(res_AICC_type1)= c("Ave", "TN", "TP", "MSE", "AUC")
res_HQC_type1 = matrix(numeric(5*num),ncol=5); colnames(res_HQC_type1) = c("Ave", "TN", "TP", "MSE", "AUC")

L <-1
E <-3
H<-3
T.s <- seq(from=0.1,to=0.5,length=L)
T.alpha <- seq(from=0.1,to=0.9,length=E)
T.H <-seq(from=0.1,to=3,length=H)

AUC <- function(score, truth) {
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  r <- rank(score, ties.method = "average")
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
ts <- Sys.time()

####################################################################################################################
##############################       GCV            #################################################################
####################################################################################################################

total_start = Sys.time()
for(w in 1:num){
  cat("Simulation(GCV):", w, "start \n")
  tic()
  
  set.seed(seed[w])
  # p=500
  t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov)  ##100  60
  
  ## beta function
  m = runif(n, min = 0, max = 1)
  type = 1
  n.b = n
  NZ.p = p*0.2
  oo = sample(c(1:p), NZ.p, replace=FALSE)
  
  ch.BETAf = matrix(numeric(p*n), ncol=p)
  ch.BETA = matrix(numeric(NZ.p*n), ncol=NZ.p)
  for (j in 1:NZ.p){
    ch.BETA[((n-n.b)+1):n,j]<-seq(from=1, to=3, length=n.b)
  }
  
  if(type==1){
    ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA    ######## M1
  }else if(type==2){
    ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    ######## M2
  }
  
  t_Y = matrix(numeric(n*1),ncol=1)
  e = rnorm(n,0,1)
  for (i in 1:n){
    t_Y[i] = t_X[i,]%*%ch.BETAf[i,]+e[i]
  }
  colnames(t_X)<-paste("X",c(1:p))
  
  
  sam <- sample(1:n, n*0.9)
  sam20 <- setdiff(1:n, sam)
  rge<-scale(t_X[sam,],center=TRUE,scale=TRUE)
  tge<-scale(t_Y[sam],center=TRUE,scale=FALSE)
  me <- m[sam]
  M.obs<-1:length(me)
  
  # p=500
  # variance = apply(t_X, 2, function(x) var(x))
  # rge = rge[, order(variance, decreasing = TRUE)[1:p]]
  # ch.BETAf = ch.BETAf[, order(variance, decreasing = TRUE)[1:p]]
  
  #################### NP start
  parameterMatrix <- NULL
  for(i in 1:length(T.s)){
    for(j in 1:length(T.alpha)){
      for(l in 1:length(T.H)){
        parameterMatrix <-rbind(parameterMatrix,c(T.s[i],T.alpha[j],T.H[l]))
      }
    }
  }
  PARA<-data.frame(matrix(numeric((ncol(parameterMatrix)+length(me))*nrow(parameterMatrix)),nrow=nrow(parameterMatrix)))
  PARA[,1:3]<-parameterMatrix
  
  for (tp in 1:nrow(PARA)){
    
    CV.ERROR.al<-matrix(numeric(1*length(me)),ncol=length(me))
    
    ### opne i
    for (i in 1:length(me)) {
      weight  <- exp(-(me - me[i])^2 / PARA[tp, 3])
      sw <- sqrt(weight)
      
      K.cv_rge <- rge * sw
      K.cv_tge <- tge * sw
      
      fit <- glmnet(K.cv_rge, K.cv_tge, alpha=PARA[tp,2], lambda=PARA[tp,1])
      prediction <- predict(fit, newx = K.cv_rge, s = PARA[tp,1])
      beta <- fit$beta
      
      sigma = PARA[tp,1] * (1 - PARA[tp,2] + PARA[tp,2] / pmax(abs(beta), 1e-12))
      A <- crossprod(K.cv_rge) 
      diag(A) <- diag(A) + sigma 
      H <- K.cv_rge %*% solve(A, t(K.cv_rge)) 
      
      CV.ERROR.al[i] <- sum((K.cv_tge - prediction)^2) / (length(K.cv_tge) - sum(diag(H)))^2
    } 
    PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
  }
  
  
  opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
  sct.TP <- PARA[opt_idx, 1:3, drop = FALSE]
  colnames(sct.TP) <- c("lambda", "alpha", "h")
  
  RESULT<-matrix(numeric(length(me)*(4)),ncol=(4))
  BETA<-matrix(numeric(p*length(me)),nrow=length(me))
  
  
  for (i in 1:length(me)){
    weight  <- exp(-(me - me[i])^2 / sct.TP[i,3])
    sw <- sqrt(weight)
    K_rge <- rge * sw
    K_tge <- tge * sw
    
    fit <- glmnet(K_rge, K_tge, alpha=sct.TP[i,2], lambda = sct.TP[i,1])
    BETA[i,] <- coef(fit)[-1]
    
    RESULT[i,1]<-sum(BETA[i,ch.BETAf[sam[i],]==0]==0)/sum(ch.BETAf[sam[i],]==0) #T.N
    RESULT[i,2]<-sum(BETA[i,ch.BETAf[sam[i],]!=0]!=0)/sum(ch.BETAf[sam[i],]!=0) #T.P
    RESULT[i,3]<-t(tge[i]-rge[i,]%*%(BETA[i,]))%*%(tge[i]-rge[i,]%*%(BETA[i,]))
    truth <- as.integer(ch.BETAf[sam[i], ] != 0)  
    score <- abs(BETA[i,])
    RESULT[i,4] <- AUC(score, truth) 
  } 
  
  ######################
  res_GCV_type1[w,1] = mean((RESULT[is.na(RESULT[,2])==FALSE,1]+RESULT[is.na(RESULT[,2])==FALSE,2])/2)
  res_GCV_type1[w,2] = mean(RESULT[,1], na.rm=TRUE) #T.N
  res_GCV_type1[w,3] = mean(RESULT[,2], na.rm=TRUE) #T.P
  res_GCV_type1[w,4] = mean(RESULT[,3], na.rm=TRUE) #MSE
  res_GCV_type1[w,5] = mean(RESULT[,4], na.rm=TRUE) #AUC
  toc()
}

total_end = Sys.time()
total = total_end - total_start
cat("total time: ", total, "\n")

sct_GCV = sct.TP
####################################################################################################################
##############################       AIC            #################################################################
####################################################################################################################

total_start = Sys.time()
for(w in 1:num){
  cat("Simulation(AIC):", w, "start \n")
  
  tic()
  set.seed(seed[w])
  # p=500
  t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov)  ##100  60
  
  ## beta function
  m = runif(n, min = 0, max = 1)
  type = 1
  n.b = n
  NZ.p = p*0.2
  oo = sample(c(1:p), NZ.p, replace=FALSE)
  
  ch.BETAf = matrix(numeric(p*n), ncol=p)
  ch.BETA = matrix(numeric(NZ.p*n), ncol=NZ.p)
  for (j in 1:NZ.p){
    ch.BETA[((n-n.b)+1):n,j]<-seq(from=1, to=3, length=n.b)
  }
  
  if(type==1){
    ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA    ######## M1
  }else if(type==2){
    ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    ######## M2
  }
  
  # y = seq(1,1,length=length(m))
  # plot(m,y,ylim=c(0,0.5),ylab ="",xlab ="", yaxt='n',main="M1")
  # abline(v=m)
  # par(new=TRUE)
  # plot(m,ch.BETAf[,oo[2]])
  
  t_Y = matrix(numeric(n*1),ncol=1)
  e = rnorm(n,0,1)
  for (i in 1:n){
    t_Y[i] = t_X[i,]%*%ch.BETAf[i,]+e[i]
  }
  colnames(t_X)<-paste("X",c(1:p))
  
  
  sam <- sample(1:n, n*0.9)
  sam20 <- setdiff(1:n, sam)
  rge<-scale(t_X[sam,],center=TRUE,scale=TRUE)
  tge<-scale(t_Y[sam],center=TRUE,scale=FALSE)
  me <- m[sam]
  M.obs<-1:length(me)
  
  #################### NP start
  parameterMatrix <- NULL
  for(i in 1:length(T.s)){
    for(j in 1:length(T.alpha)){
      for(l in 1:length(T.H)){
        parameterMatrix <-rbind(parameterMatrix,c(T.s[i],T.alpha[j],T.H[l]))
      }
    }
  }
  PARA<-data.frame(matrix(numeric((ncol(parameterMatrix)+length(me))*nrow(parameterMatrix)),nrow=nrow(parameterMatrix)))
  PARA[,1:3]<-parameterMatrix
  
  for (tp in 1:nrow(PARA)){
    
    CV.ERROR.al<-matrix(numeric(1*length(me)),ncol=length(me))
    
    ### opne i
    for (i in 1:length(me)){
      weight  <- exp(-(me - me[i])^2 / PARA[tp, 3])
      sw <- sqrt(weight)
      
      K.cv_rge <- rge * sw
      K.cv_tge <- tge * sw
      
      fit <- ic.glmnet(K.cv_rge, K.cv_tge, alpha=PARA[tp,2], lambda=PARA[tp,1], crit = "aic")
      CV.ERROR.al[i] <- fit$ic[2]
    } 
    PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
  }
  
  
  opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
  sct.TP <- PARA[opt_idx, 1:3, drop = FALSE]
  colnames(sct.TP) <- c("lambda", "alpha", "h")
  
  RESULT<-matrix(numeric(length(me)*(4)),ncol=(4))
  BETA<-matrix(numeric(p*length(me)),nrow=length(me))
  
  
  for (i in 1:length(me)){
    weight  <- exp(-(me - me[i])^2 / sct.TP[i,3])
    sw <- sqrt(weight)
    K_rge <- rge * sw
    K_tge <- tge * sw
    
    fit <- glmnet(K_rge, K_tge, alpha=sct.TP[i,2], lambda = sct.TP[i,1])
    BETA[i,]<-coef(fit)[-1]
    
    RESULT[i,1]<-sum(BETA[i,ch.BETAf[sam[i],]==0]==0)/sum(ch.BETAf[sam[i],]==0) #T.N
    RESULT[i,2]<-sum(BETA[i,ch.BETAf[sam[i],]!=0]!=0)/sum(ch.BETAf[sam[i],]!=0) #T.P
    RESULT[i,3]<-t(tge[i]-rge[i,]%*%(BETA[i,]))%*%(tge[i]-rge[i,]%*%(BETA[i,]))
    truth <- as.integer(ch.BETAf[sam[i], ] != 0)  
    score <- abs(BETA[i,])
    RESULT[i,4] <- AUC(score, truth)
  } 
  
  ######################
  res_AIC_type1[w,1] = mean((RESULT[is.na(RESULT[,2])==FALSE,1]+RESULT[is.na(RESULT[,2])==FALSE,2])/2)
  res_AIC_type1[w,2] = mean(RESULT[,1], na.rm=TRUE) #T.N
  res_AIC_type1[w,3] = mean(RESULT[,2], na.rm=TRUE) #T.P
  res_AIC_type1[w,4] = mean(RESULT[,3], na.rm=TRUE) #MSE
  res_AIC_type1[w,5] = mean(RESULT[,4], na.rm=TRUE) #AUC
  toc()
}

total_end = Sys.time()
total = total_end - total_start
cat("total time: ", total, "\n")

####################################################################################################################
#############################       BIC            ##############################################################
####################################################################################################################

total_start = Sys.time()
for(w in 1:num){
  cat("Simulation(BIC):", w, "start \n")
  tic()
  
  set.seed(seed[w])
  # p=500
  t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov)  ##100  60
  
  ## beta function
  m = runif(n, min = 0, max = 1)
  type = 1
  n.b = n
  NZ.p = p*0.2
  oo = sample(c(1:p), NZ.p, replace=FALSE)
  
  ch.BETAf = matrix(numeric(p*n), ncol=p)
  ch.BETA = matrix(numeric(NZ.p*n), ncol=NZ.p)
  for (j in 1:NZ.p){
    ch.BETA[((n-n.b)+1):n,j]<-seq(from=1, to=3, length=n.b)
  }
  
  if(type==1){
    ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA    ######## M1
  }else if(type==2){
    ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    ######## M2
  }
  
  # y = seq(1,1,length=length(m))
  # plot(m,y,ylim=c(0,0.5),ylab ="",xlab ="", yaxt='n',main="M1")
  # abline(v=m)
  # par(new=TRUE)
  # plot(m,ch.BETAf[,oo[2]])
  
  t_Y = matrix(numeric(n*1),ncol=1)
  e = rnorm(n,0,1)
  for (i in 1:n){
    t_Y[i] = t_X[i,]%*%ch.BETAf[i,]+e[i]
  }
  colnames(t_X)<-paste("X",c(1:p))
  
  
  sam <- sample(1:n, n*0.9)
  sam20 <- setdiff(1:n, sam)
  rge<-scale(t_X[sam,],center=TRUE,scale=TRUE)
  tge<-scale(t_Y[sam],center=TRUE,scale=FALSE)
  me <- m[sam]
  M.obs<-1:length(me)
  
  #################### NP start
  parameterMatrix <- NULL
  for(i in 1:length(T.s)){
    for(j in 1:length(T.alpha)){
      for(l in 1:length(T.H)){
        parameterMatrix <-rbind(parameterMatrix,c(T.s[i],T.alpha[j],T.H[l]))
      }
    }
  }
  PARA<-data.frame(matrix(numeric((ncol(parameterMatrix)+length(me))*nrow(parameterMatrix)),nrow=nrow(parameterMatrix)))
  PARA[,1:3]<-parameterMatrix
  
  for (tp in 1:nrow(PARA)){
    
    CV.ERROR.al<-matrix(numeric(1*length(me)),ncol=length(me))
    
    ### opne i
    for (i in 1:length(me)){
      weight  <- exp(-(me - me[i])^2 / PARA[tp, 3])
      sw <- sqrt(weight)
      
      K.cv_rge <- rge * sw
      K.cv_tge <- tge * sw
      
      fit <- ic.glmnet(K.cv_rge, K.cv_tge, alpha=PARA[tp,2], lambda=PARA[tp,1], crit = "bic")
      CV.ERROR.al[i] <- fit$ic[1]
    } 
    PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
  }
  
  
  opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
  sct.TP <- PARA[opt_idx, 1:3, drop = FALSE]
  colnames(sct.TP) <- c("lambda", "alpha", "h")
  
  RESULT<-matrix(numeric(length(me)*(4)),ncol=(4))
  BETA<-matrix(numeric(p*length(me)),nrow=length(me))
  
  
  for (i in 1:length(me)){
    weight  <- exp(-(me - me[i])^2 / sct.TP[i,3])
    sw <- sqrt(weight)
    K_rge <- rge * sw
    K_tge <- tge * sw
    
    fit <- glmnet(K_rge, K_tge, alpha=sct.TP[i,2], lambda = sct.TP[i,1])
    BETA[i,]<-coef(fit)[-1]
    
    RESULT[i,1]<-sum(BETA[i,ch.BETAf[sam[i],]==0]==0)/sum(ch.BETAf[sam[i],]==0) #T.N
    RESULT[i,2]<-sum(BETA[i,ch.BETAf[sam[i],]!=0]!=0)/sum(ch.BETAf[sam[i],]!=0) #T.P
    RESULT[i,3]<-t(tge[i]-rge[i,]%*%(BETA[i,]))%*%(tge[i]-rge[i,]%*%(BETA[i,]))
    truth <- as.integer(ch.BETAf[sam[i], ] != 0)  
    score <- abs(BETA[i,])
    RESULT[i,4] <- AUC(score, truth)
  } 
  
  ######################
  res_BIC_type1[w,1] = mean((RESULT[is.na(RESULT[,2])==FALSE,1]+RESULT[is.na(RESULT[,2])==FALSE,2])/2)
  res_BIC_type1[w,2] = mean(RESULT[,1], na.rm=TRUE) #T.N
  res_BIC_type1[w,3] = mean(RESULT[,2], na.rm=TRUE) #T.P
  res_BIC_type1[w,4] = mean(RESULT[,3], na.rm=TRUE) #MSE
  res_BIC_type1[w,5] = mean(RESULT[,4], na.rm=TRUE) #AUC
  toc()
}

total_end = Sys.time()
total = total_end - total_start
cat("total time: ", total, "\n")

####################################################################################################################
##############################       AICC            ############################################################
####################################################################################################################

total_start = Sys.time()
for(w in 1:num){
  cat("Simulation(AICC):", w, "start \n")
  tic()
  
  set.seed(seed[w])
  # p=500
  t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov)  ##100  60
  
  ## beta function
  m = runif(n, min = 0, max = 1)
  type = 1
  n.b = n
  NZ.p = p*0.2
  oo = sample(c(1:p), NZ.p, replace=FALSE)
  
  ch.BETAf = matrix(numeric(p*n), ncol=p)
  ch.BETA = matrix(numeric(NZ.p*n), ncol=NZ.p)
  for (j in 1:NZ.p){
    ch.BETA[((n-n.b)+1):n,j]<-seq(from=1, to=3, length=n.b)
  }
  
  if(type==1){
    ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA    ######## M1
  }else if(type==2){
    ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    ######## M2
  }
  
  # y = seq(1,1,length=length(m))
  # plot(m,y,ylim=c(0,0.5),ylab ="",xlab ="", yaxt='n',main="M1")
  # abline(v=m)
  # par(new=TRUE)
  # plot(m,ch.BETAf[,oo[2]])
  
  t_Y = matrix(numeric(n*1),ncol=1)
  e = rnorm(n,0,1)
  for (i in 1:n){
    t_Y[i] = t_X[i,]%*%ch.BETAf[i,]+e[i]
  }
  colnames(t_X)<-paste("X",c(1:p))
  
  
  sam <- sample(1:n, n*0.9)
  sam20 <- setdiff(1:n, sam)
  rge<-scale(t_X[sam,],center=TRUE,scale=TRUE)
  tge<-scale(t_Y[sam],center=TRUE,scale=FALSE)
  me <- m[sam]
  M.obs<-1:length(me)
  
  #################### NP start
  parameterMatrix <- NULL
  for(i in 1:length(T.s)){
    for(j in 1:length(T.alpha)){
      for(l in 1:length(T.H)){
        parameterMatrix <-rbind(parameterMatrix,c(T.s[i],T.alpha[j],T.H[l]))
      }
    }
  }
  PARA<-data.frame(matrix(numeric((ncol(parameterMatrix)+length(me))*nrow(parameterMatrix)),nrow=nrow(parameterMatrix)))
  PARA[,1:3]<-parameterMatrix
  
  for (tp in 1:nrow(PARA)){
    
    CV.ERROR.al<-matrix(numeric(1*length(me)),ncol=length(me))
    
    ### opne i
    for (i in 1:length(me)){
      weight  <- exp(-(me - me[i])^2 / PARA[tp, 3])
      sw <- sqrt(weight)
      
      K.cv_rge <- rge * sw
      K.cv_tge <- tge * sw
      
      fit <- ic.glmnet(K.cv_rge, K.cv_tge, alpha=PARA[tp,2], lambda=PARA[tp,1], crit = "aicc")
      CV.ERROR.al[i] <- fit$ic[3]
    } 
    PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
  }
  
  
  opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
  sct.TP <- PARA[opt_idx, 1:3, drop = FALSE]
  colnames(sct.TP) <- c("lambda", "alpha", "h")
  
  RESULT<-matrix(numeric(length(me)*(4)),ncol=(4))
  BETA<-matrix(numeric(p*length(me)),nrow=length(me))
  
  
  for (i in 1:length(me)){
    weight  <- exp(-(me - me[i])^2 / sct.TP[i,3])
    sw <- sqrt(weight)
    K_rge <- rge * sw
    K_tge <- tge * sw
    
    fit <- glmnet(K_rge, K_tge, alpha=sct.TP[i,2], lambda = sct.TP[i,1])
    BETA[i,]<-coef(fit)[-1]
    
    RESULT[i,1]<-sum(BETA[i,ch.BETAf[sam[i],]==0]==0)/sum(ch.BETAf[sam[i],]==0) #T.N
    RESULT[i,2]<-sum(BETA[i,ch.BETAf[sam[i],]!=0]!=0)/sum(ch.BETAf[sam[i],]!=0) #T.P
    RESULT[i,3]<-t(tge[i]-rge[i,]%*%(BETA[i,]))%*%(tge[i]-rge[i,]%*%(BETA[i,]))
    truth <- as.integer(ch.BETAf[sam[i], ] != 0)  
    score <- abs(BETA[i,])
    RESULT[i,4] <- AUC(score, truth)
  } 
  
  ######################
  res_AICC_type1[w,1] = mean((RESULT[is.na(RESULT[,2])==FALSE,1]+RESULT[is.na(RESULT[,2])==FALSE,2])/2)
  res_AICC_type1[w,2] = mean(RESULT[,1], na.rm=TRUE) #T.N
  res_AICC_type1[w,3] = mean(RESULT[,2], na.rm=TRUE) #T.P
  res_AICC_type1[w,4] = mean(RESULT[,3], na.rm=TRUE) #MSE
  res_AICC_type1[w,5] = mean(RESULT[,4], na.rm=TRUE) #AUC
  toc()
}

total_end = Sys.time()
total = total_end - total_start
cat("total time: ", total, "\n")

####################################################################################################################
###############################       HQC             ######################################################
####################################################################################################################

total_start = Sys.time()
for(w in 1:num){
  
  cat("Simulation(HQC):", w, "start \n")
  tic()
  
  set.seed(seed[w])
  # p=500
  t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov)  ##100  60
  
  ## beta function
  m = runif(n, min = 0, max = 1)
  type = 1
  n.b = n
  NZ.p = p*0.2
  oo = sample(c(1:p), NZ.p, replace=FALSE)
  
  ch.BETAf = matrix(numeric(p*n), ncol=p)
  ch.BETA = matrix(numeric(NZ.p*n), ncol=NZ.p)
  for (j in 1:NZ.p){
    ch.BETA[((n-n.b)+1):n,j]<-seq(from=1, to=3, length=n.b)
  }
  
  if(type==1){
    ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA    ######## M1
  }else if(type==2){
    ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    ######## M2
  }
  
  # y = seq(1,1,length=length(m))
  # plot(m,y,ylim=c(0,0.5),ylab ="",xlab ="", yaxt='n',main="M1")
  # abline(v=m)
  # par(new=TRUE)
  # plot(m,ch.BETAf[,oo[2]])
  
  t_Y = matrix(numeric(n*1),ncol=1)
  e = rnorm(n,0,1)
  for (i in 1:n){
    t_Y[i] = t_X[i,]%*%ch.BETAf[i,]+e[i]
  }
  colnames(t_X)<-paste("X",c(1:p))
  
  
  sam <- sample(1:n, n*0.9)
  sam20 <- setdiff(1:n, sam)
  rge<-scale(t_X[sam,],center=TRUE,scale=TRUE)
  tge<-scale(t_Y[sam],center=TRUE,scale=FALSE)
  me <- m[sam]
  M.obs<-1:length(me)
  
  #################### NP start
  parameterMatrix <- NULL
  for(i in 1:length(T.s)){
    for(j in 1:length(T.alpha)){
      for(l in 1:length(T.H)){
        parameterMatrix <-rbind(parameterMatrix,c(T.s[i],T.alpha[j],T.H[l]))
      }
    }
  }
  PARA<-data.frame(matrix(numeric((ncol(parameterMatrix)+length(me))*nrow(parameterMatrix)),nrow=nrow(parameterMatrix)))
  PARA[,1:3]<-parameterMatrix
  
  for (tp in 1:nrow(PARA)){
    
    CV.ERROR.al<-matrix(numeric(1*length(me)),ncol=length(me))
    
    ### opne i
    for (i in 1:length(me)){
      weight  <- exp(-(me - me[i])^2 / PARA[tp, 3])
      sw <- sqrt(weight)
      
      K.cv_rge <- rge * sw
      K.cv_tge <- tge * sw
      
      fit <- ic.glmnet(K.cv_rge, K.cv_tge, alpha=PARA[tp,2], lambda=PARA[tp,1], crit = "hqc")
      CV.ERROR.al[i] <- fit$ic[4]
    } 
    PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
  }
  
  
  opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
  sct.TP <- PARA[opt_idx, 1:3, drop = FALSE]
  colnames(sct.TP) <- c("lambda", "alpha", "h")
  
  RESULT<-matrix(numeric(length(me)*(4)),ncol=(4))
  BETA<-matrix(numeric(p*length(me)),nrow=length(me))
  
  
  for (i in 1:length(me)){
    weight  <- exp(-(me - me[i])^2 / sct.TP[i,3])
    sw <- sqrt(weight)
    K_rge <- rge * sw
    K_tge <- tge * sw
    
    fit <- glmnet(K_rge, K_tge, alpha=sct.TP[i,2], lambda = sct.TP[i,1])
    BETA[i,]<-coef(fit)[-1]
    
    RESULT[i,1]<-sum(BETA[i,ch.BETAf[sam[i],]==0]==0)/sum(ch.BETAf[sam[i],]==0) #T.N
    RESULT[i,2]<-sum(BETA[i,ch.BETAf[sam[i],]!=0]!=0)/sum(ch.BETAf[sam[i],]!=0) #T.P
    RESULT[i,3]<-t(tge[i]-rge[i,]%*%(BETA[i,]))%*%(tge[i]-rge[i,]%*%(BETA[i,]))
    truth <- as.integer(ch.BETAf[sam[i], ] != 0)  
    score <- abs(BETA[i,])
    RESULT[i,4] <- AUC(score, truth)
  } 
  
  ######################
  res_HQC_type1[w,1] = mean((RESULT[is.na(RESULT[,2])==FALSE,1]+RESULT[is.na(RESULT[,2])==FALSE,2])/2)
  res_HQC_type1[w,2] = mean(RESULT[,1], na.rm=TRUE) #T.N
  res_HQC_type1[w,3] = mean(RESULT[,2], na.rm=TRUE) #T.P
  res_HQC_type1[w,4] = mean(RESULT[,3], na.rm=TRUE) #MSE
  res_HQC_type1[w,5] = mean(RESULT[,4], na.rm=TRUE) #AUC
  toc()
}

total_end = Sys.time()
total = total_end - total_start
cat("total time: ", total, "\n")
te <- Sys.time()
t <- te - ts
cat("total time: ", t, "\n")

round(colMeans(res_GCV_type1), 3)
round(colMeans(res_AIC_type1), 3)
round(colMeans(res_BIC_type1), 3)
round(colMeans(res_AICC_type1), 3)
round(colMeans(res_HQC_type1), 3)




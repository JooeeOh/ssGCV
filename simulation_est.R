rm(list=ls())
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

res_GCV_type1 = matrix(numeric(5*num),ncol=5); colnames(res_GCV_type1) = c("Ave", "TN", "TP", "MSE", "AUC")


L <-3
E <-3
H <-3
T.s <- seq(from=0.1,to=0.5,length=L)
T.alpha <- seq(from=0.1,to=0.9,length=E)
T.H <-seq(from=0.1,to=3,length=H)


ts <- Sys.time()

## generate data
t_X = mvrnorm(n, mu = numeric(1*p), Sigma = cov) 
m = runif(n, min = 0, max = 1)

## beta function (2 type)
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
  ch.BETAf[order(m, decreasing=TRUE),oo]<-ch.BETA 
}else if(type==2){
  ch.BETAf[order(m, decreasing=FALSE),oo]<-ch.BETA    
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

## NetworkProfiler
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
  for (i in 1:length(me)) {
    CV.ERROR.al[i] <- ssGCV(X = rge,
                            y = tge,
                            mod_vec = me,
                            target = i,
                            lambda = PARA[tp, 1],
                            alpha = PARA[tp, 2],
                            bandwidth = PARA[tp, 3])
  } 
  PARA[tp,c(4:(4+nrow(rge)-1))]<-CV.ERROR.al
}


opt_idx <- apply(PARA[, 4:(3 + length(me)), drop = FALSE], 2, which.min)  
opt_mx <- PARA[opt_idx, 1:3, drop = FALSE]
colnames(opt_mx) <- c("lambda", "alpha", "h"); rownames(opt_mx) <- NULL


BETA <- matrix(numeric(p*length(me)),nrow=length(me))
for (i in 1:length(me)){
  weight  <- exp(-(me - me[i])^2 / opt_mx[i,3])
  sw <- sqrt(weight)
  K_rge <- rge * sw
  K_tge <- tge * sw
  
  fit <- glmnet(K_rge, K_tge, alpha=opt_mx[i,2], lambda = opt_mx[i,1])
  BETA[i,] <- coef(fit)[-1]
} 

BETA
opt_mx







###########################################################
# Set up
###########################################################
library(causaldrf)
library(GGally)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(CAM)
# Gaussian kernel
K = function(x,y,s){
  m = matrix(rep(0,length(x)*length(y)),nrow=length(x))
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m[i,j] = exp(-(x[i]-y[j])^2/(2*s^2))
    }
  }
  return(m)
}
# Median of interpoint distance
dis = function(){
  disx = disz = rep(0,n*(n-1)/2)
  id = 1
  for(i in c(1:(n-1))){
    for(j in c((i+1):n)){
      disx[id] = abs(X[i]-X[j])
      disz[id] = abs(Z[i]-Z[j])
      id = id + 1
    }
  }
  sx = median(disx); sz = median(disz)
  return(c(sx,sz))
}

###########################################################
# EDA
###########################################################
data("nmes_data")
nmes_data = nmes_data[nmes_data$TOTALEXP>0,] # Focus on smokers
dim(nmes_data)
summary(nmes_data)
cor(nmes_data[,c(1,2,3,4,12)])

# Vairbles
X = log(nmes_data$packyears)
Z = nmes_data$LASTAGE
Y = log(nmes_data$TOTALEXP)

# Sample allocation
set.seed(7654)
n = 1000
n1 = n2 = 300; m = n - (n1+n2)
id = c(1:nrow(nmes_data))
id = sample(id,n)
X = X[id]; Z = Z[id]; Y = Y[id]
id = c(1:n)
id1 = sample(id,n1); id = setdiff(id,id1)
id2 = sample(id,n1); id = setdiff(id,id2)
id3 = id
id12 = c(id1,id2); n12 = n1 + n2

summary(cbind(X,Y,Z))

s = dis()
sx = s[1]; sz = s[2]

xx = seq(min(X),max(X),length.out=n+2)[-c(1,n+2)]
gamma = 2.9
alpha1 = 1/sqrt(n1)
alpha2 = 1/sqrt(n2)
alpha = 1/sqrt(n12)
xi = 1/sqrt(m)

# KAR
aY = t(Y[id3]) + (sqrt(gamma)-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
K1ZZ = K(Z[id1],Z[id1],sz); K13ZZ = K(Z[id1],Z[id3],sz)
K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id1],sx); K1XX = K(X[id1],X[id1],sx)
K34XX = K(X[id3],xx,sx); K14XX = K(X[id1],xx,sx)
KaXX = K3XX + 
  (sqrt(gamma)-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
  (sqrt(gamma)-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
y_KAR = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
y_KAR = t(y_KAR)

# KPA
aY = t(Y[id3]) + (-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
KaXX = K3XX + 
  (-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
  (-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
y_KPA = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
y_KPA = t(y_KPA)

# KAR-2
aY = t(Y[id3]) + (sqrt(gamma)-1)*t(Y[id12])%*%solve(K(Z[id12],Z[id12],sz)+n12*alpha*diag(n12))%*%K(Z[id12],Z[id3],sz)
K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
KaXX = K3XX + 
  (sqrt(gamma)-1)^2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K1XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX) +
  (sqrt(gamma)-1)*2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX)
Ka4 = K34XX + K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K14XX
y_KAR2 = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
y_KAR2 = t(y_KAR2)

# KIV
W = K1XX%*%solve(K1ZZ + n12*alpha*diag(n12))%*%K13ZZ
y_KIV = t(ginv(W%*%t(W)+m*xi*K1XX)%*%W%*%Y[id3])%*%K14XX
y_KIV = t(y_KIV)

#KernelReg
KXX = K(X,X,sx)
KX4 = K(X,xx,sx)
y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
y_KReg = t(y_KReg)

###########################################################
# Linear estimate
###########################################################
# AR
aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_AR = c(b_AR)*(xx-mean(X)) + mean(Y)

# IV
aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_IV = c(b_IV)*(xx-mean(X)) + mean(Y)

# PA
aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_PA = c(b_PA)*(xx-mean(X)) + mean(Y)

# Basic
b = solve(t(X)%*%X)%*%t(aX)%*%Y
y_b = c(b)*(xx-mean(X)) + mean(Y)


###########################################################
# Visualisation
###########################################################
Mycol = brewer.pal(10, 'Paired') 
Mycol = c(Mycol[c(6,5)],Mycol[-c(5,6)])
dat2 = rbind(data.frame('X'=X,'Y'=Y,'Method'=rep("Observations",n)),
             data.frame('X'=xx,'Y'=y_KAR,'Method'=rep('KAR',n)),
             data.frame('X'=xx,'Y'=y_KAR2,'Method'=rep('KAR.2',n)),
             data.frame('X'=xx,'Y'=y_KIV,'Method'=rep('KIV',n)),
             data.frame('X'=xx,'Y'=y_KPA,'Method'=rep('KPA',n)),
             data.frame('X'=xx,'Y'=y_KReg,'Method'=rep('KReg',n)),
             data.frame('X'=xx,'Y'=y_AR,'Method'=rep('AR',n)),
             data.frame('X'=xx,'Y'=y_IV,'Method'=rep('IV',n)),
             data.frame('X'=xx,'Y'=y_PA,'Method'=rep('PA',n)),
             data.frame('X'=xx,'Y'=y_b,'Method'=rep('OLS',n)))
gs1 = ggplot(dat2)+
  geom_point(size=0.45,alpha=0.5,aes(x=X,y=Y,color=Method))+
  theme_classic()+
  ylab("Medical Expenditures (log)")+
  xlab("Smoking Amount (log)")+
  theme(text = element_text(size = 16))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6],
                              'AR'=Mycol[7],'IV'=Mycol[8],'PA'=Mycol[9],'OLS'=Mycol[10]))+
  guides(color = guide_legend(override.aes = list(size = 4)))

# ICA-based results
X = log(nmes_data$packyears)
Z = nmes_data$LASTAGE
Y = log(nmes_data$TOTALEXP)
CAM(cbind(X,Y,Z))$Adj

###########################################################
# Different groups for training and testing
###########################################################
n = 1000
n1 = n2 = 300; m = n - (n1+n2)
gamma = 2.9
# Testing samples
te = nmes_data[nmes_data$MALE==0,]
xx = log(te$packyears)
zz = te$LASTAGE
yy = log(te$TOTALEXP)

set.seed(7654)
n_repeat = 50
pre3 = matrix(rep(0,n_repeat*9),nrow=n_repeat)

for(k in 1:n_repeat){
  # Training samples
  tr = nmes_data[nmes_data$MALE==1,]
  X = log(tr$packyears)
  Z = tr$LASTAGE
  Y = log(tr$TOTALEXP)
  
  id = c(1:nrow(tr))
  id = sample(id,n)
  X = X[id]; Z = Z[id]; Y = Y[id]
  id = c(1:n)
  id1 = sample(id,n1); id = setdiff(id,id1)
  id2 = sample(id,n1); id = setdiff(id,id2)
  id3 = id
  id12 = c(id1,id2); n12 = n1 + n2
  s = dis()
  sx = s[1]; sz = s[2]
  
  ###########################################################
  # Kernel
  ###########################################################
  # KAR
  aY = t(Y[id3]) + (sqrt(gamma)-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
  K1ZZ = K(Z[id1],Z[id1],sz); K13ZZ = K(Z[id1],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id1],sx); K1XX = K(X[id1],X[id1],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id1],xx,sx)
  KaXX = K3XX + 
    (sqrt(gamma)-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
    (sqrt(gamma)-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
  y_KAR = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KAR = t(y_KAR)
  # Error
  pre3[k,1] = mean((yy-y_KAR)^2)
  
  # KAR-2
  aY = t(Y[id3]) + (sqrt(gamma)-1)*t(Y[id12])%*%solve(K(Z[id12],Z[id12],sz)+n12*alpha*diag(n12))%*%K(Z[id12],Z[id3],sz)
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  KaXX = K3XX + 
    (sqrt(gamma)-1)^2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K1XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX) +
    (sqrt(gamma)-1)*2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX)
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K14XX
  y_KAR2 = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KAR2 = t(y_KAR2)
  # Error
  pre3[k,2] = mean((yy-y_KAR2)^2)
  
  # KIV
  aY = Y[id3]
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  KaXX = K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K1XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX) 
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K14XX
  y_KIV = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KIV = t(y_KIV)
  # Error
  pre3[k,3] = mean((yy-y_KIV)^2)
  
  # KPA
  aY = t(Y[id3]) + (-1)*t(Y[id12])%*%solve(K(Z[id12],Z[id12],sz)+n12*alpha*diag(n12))%*%K(Z[id12],Z[id3],sz)
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  KaXX = K3XX + 
    (-1)^2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K1XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX) +
    (-1)*2*K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%t(K31XX)
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n12*alpha*diag(n12))%*%K14XX
  y_KPA = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KPA = t(y_KPA)
  # Error
  pre3[k,4] = mean((yy-y_KPA)^2)
  
  #KernelReg
  KXX = K(X,X,sx)
  KX4 = K(X,xx,sx)
  y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
  y_KReg = t(y_KReg)
  # Error
  pre3[k,5] = mean((yy-y_KReg)^2)
  
  ###########################################################
  # Linear estimate
  ###########################################################
  # AR
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
  b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_AR = c(b_AR)*(xx-0.5)
  pre3[k,6] = mean((yy-y_AR)^2)
  
  # IV
  aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_IV = c(b_IV)*(xx-0.5)
  pre3[k,7] = mean((yy-y_IV)^2)
  
  # PA
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_PA = c(b_PA)*(xx-0.5)
  pre3[k,8] = mean((yy-y_PA)^2)
  
  # Basic
  b = solve(t(X)%*%X)%*%t(aX)%*%Y
  y_b = c(b)*(xx-mean(X)) + mean(Y)
  pre3[k,9] = mean((yy-y_b)^2)
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lpre3 = log(pre3,base=10)
colMeans(lpre3)
res = data.frame('MSE'=c(lpre3[,1],lpre3[,2],lpre3[,3],lpre3[,4],
                         lpre3[,5],lpre3[,6],lpre3[,7],lpre3[,8],lpre3[,9]),
                 'Method'=factor(rep(c('KAR','KAR.2','KIV','KPA','AR','KReg','IV','PA','OLS'),each=n_repeat),
                                 levels = c('KAR','KAR.2','KIV','KPA','AR','KReg','IV','PA','OLS')))
gs2 = ggplot(res,aes(y=MSE))+
  ylab("Prediction Error (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],
                              'AR'=Mycol[6],'IV'=Mycol[7],'PA'=Mycol[8],'OLS'=Mycol[9]))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

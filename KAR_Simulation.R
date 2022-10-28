###########################################################
# Set up 
###########################################################
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(CAM)
library(cowplot)
library(Rmisc)

###########################################################
# Build the population in Case 1
###########################################################
set.seed(123)
n1 = 2.5e2; n2 = 2.5e2; m = 2e2
n12 = sum(n1,n2)
n = sum(c(n1,n2,m))
id1 = c(1:n1); id2 = c((n1+1):(n1+n2)); id3 = c((n1+n2+1):n)
id12 = c(id1,id2)

S = matrix(c(1, 0.3, 0.2,
             0.3, 1, 0,
             0.2, 0, 1), nrow = 3, byrow = T)
tmp = mvrnorm(n,rep(0,3),S)
C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
X = pnorm((V+W)/sqrt(2))
Z = pnorm(W)-0.5
Y = 0.75*C - 0.25*Z + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)

xx = seq(0,1,length.out=n+2)[-c(1,n+2)]
ty = log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)

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
s = dis()
sx = s[1]; sz = s[2]

# kernel methods
gamma = 2
alpha1 = 1.5/sqrt(n1)
alpha2 = 1.5/sqrt(n2)
alpha = 1.5/sqrt(n12)
xi = 1.5/sqrt(m)
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
K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
W = K1XX%*%solve(K1ZZ + n12*alpha*diag(n12))%*%K13ZZ
y_KIV = t(ginv(W%*%t(W)+m*xi*K1XX)%*%W%*%Y[id3])%*%K14XX
y_KIV = t(y_KIV)
# KPA
aY = t(Y[id3]) + (-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
K1ZZ = K(Z[id1],Z[id1],sz); K13ZZ = K(Z[id1],Z[id3],sz)
K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id1],sx); K1XX = K(X[id1],X[id1],sx)
K34XX = K(X[id3],xx,sx); K14XX = K(X[id1],xx,sx)
KaXX = K3XX + 
  (-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
  (-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
y_KPA = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
y_KPA = t(y_KPA)
#KernelReg
KXX = K(X,X,sx)
KX4 = K(X,xx,sx)
y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
y_KReg = t(y_KReg)

# linear methods
# AR
aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_AR = c(b_AR)*(xx-mean(X))+mean(Y)

# IV
aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_IV = c(b_IV)*(xx-mean(X))+mean(Y)

# PA
aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
y_PA = c(b_PA)*(xx-mean(X))+mean(Y)

# Basic
b = solve(t(X)%*%X)%*%t(aX)%*%Y
y_b = c(b)*(xx-mean(X))+mean(Y)

###########################################################
# Visualisation
###########################################################
dat = data.frame(X,Y,'Method'=rep('Observations',n))
dat = rbind(dat,data.frame('X'=xx,'Y'=ty,'Method'=rep('True Model',n)))
dat = rbind(dat,data.frame('X'=xx,'Y'=y_KAR,'Method'=rep('KAR',n)),
             data.frame('X'=xx,'Y'=y_KAR2,'Method'=rep('KAR.2',n)),
             data.frame('X'=xx,'Y'=y_KIV,'Method'=rep('KIV',n)),
             data.frame('X'=xx,'Y'=y_KPA,'Method'=rep('KPA',n)),
            data.frame('X'=xx,'Y'=y_KReg,'Method'=rep('KReg',n)))
Mycol = brewer.pal(10, 'Paired') 
Mycol = c(Mycol[c(6,5)],Mycol[-c(5,6)])

g1 = ggplot(dat)+
  geom_point(size=0.45,alpha=0.5,aes(x=X,y=Y,color=Method))+
  scale_color_manual(values=c('Observations'='grey','True Model'=Mycol[1],
                              'KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6]),
                     name="Method")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 4)))

dat2 = data.frame(X,Y,'Method'=rep('Observations',n))
dat2 = rbind(dat2,data.frame('X'=xx,'Y'=ty,'Method'=rep('True Model',n)))
dat2 = rbind(dat2,data.frame('X'=xx,'Y'=y_KAR,'Method'=rep('KAR',n)),
             data.frame('X'=xx,'Y'=y_AR,'Method'=rep('AR',n)),
             data.frame('X'=xx,'Y'=y_IV,'Method'=rep('IV',n)),
             data.frame('X'=xx,'Y'=y_PA,'Method'=rep('PA',n)),
             data.frame('X'=xx,'Y'=y_b,'Method'=rep('OLS',n)))
g2 = ggplot(dat2)+
  geom_point(size=0.45,alpha=0.5,aes(x=X,y=Y,color=Method))+
  scale_color_manual(values=c('Observations'='grey','True Model'=Mycol[1],'KAR'=Mycol[2],
                              'AR'=Mycol[7],'IV'=Mycol[4],'PA'=Mycol[5],'OLS'=Mycol[6]),
                     name="Method")+
  theme_classic()+
  theme(text = element_text(size = 16))+
  guides(color = guide_legend(override.aes = list(size = 4)))


###########################################################
# Repeated experiments in Case 1 (Case above)
###########################################################
set.seed(8217)
n_repeat = 50
mse1 = matrix(rep(0,n_repeat*9),nrow=n_repeat)
res_bs = rep(0,n_repeat)
for(k in 1:n_repeat){
  tmp = mvrnorm(n,rep(0,3),S)
  C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
  X = pnorm((V+W)/sqrt(2))
  Z = pnorm(W)-0.5
  Y = 0.75*C - 0.25*Z + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
  xx = seq(0,1,length.out=n+2)[-c(1,n+2)]
  ty = log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)
  
  s = dis()
  sx = s[1]; sz = s[2]
  
  ###########################################################
  # Kernel
  ###########################################################
  gamma = 2
  alpha1 = 1.5/sqrt(n1)
  alpha2 = 1.5/sqrt(n2)
  alpha = 1.5/sqrt(n12)
  xi = 1.5/sqrt(m)
  
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
  mse1[k,1] = mean((ty-y_KAR)^2)
  
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
  mse1[k,2] = mean((ty-y_KAR2)^2)
  
  # KIV
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  W = K1XX%*%solve(K1ZZ + n12*alpha*diag(n12))%*%K13ZZ
  y_KIV = t(ginv(W%*%t(W)+m*xi*K1XX)%*%W%*%Y[id3])%*%K14XX
  y_KIV = t(y_KIV)
  # Error
  mse1[k,3] = mean((ty-y_KIV)^2)
  
  # KPA
  aY = t(Y[id3]) + (-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
  K1ZZ = K(Z[id1],Z[id1],sz); K13ZZ = K(Z[id1],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id1],sx); K1XX = K(X[id1],X[id1],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id1],xx,sx)
  KaXX = K3XX + 
    (-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
    (-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
  y_KPA = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KPA = t(y_KPA)
  # Error
  mse1[k,4] = mean((ty-y_KPA)^2)
  
  #KernelReg
  KXX = K(X,X,sx)
  KX4 = K(X,xx,sx)
  y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
  y_KReg = t(y_KReg)
  # Error
  mse1[k,5] = mean((ty-y_KReg)^2)
  
  ###########################################################
  # Linear estimate
  ###########################################################
  # AR
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
  b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_AR = c(b_AR)*(xx-mean(X))+mean(Y)
  mse1[k,6] = mean((ty-y_AR)^2)
  
  # IV
  aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_IV = c(b_IV)*(xx-mean(X))+mean(Y)
  mse1[k,7] = mean((ty-y_IV)^2)
  
  # PA
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_PA = c(b_PA)*(xx-mean(X))+mean(Y)
  mse1[k,8] = mean((ty-y_PA)^2)
  
  # Basic
  b = solve(t(X)%*%X)%*%t(aX)%*%Y
  y_b = c(b)*(xx-mean(X))+mean(Y)
  mse1[k,9] = mean((ty-y_b)^2)
  
  # ICA-based methods
  bs = CAM(cbind(Z,X,Y))
  res_bs[k] = bs$Adj[2,3]
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lmse1 = log(mse1,base=10)
colMeans(lmse1)
mean(res_bs)

###########################################################
# Repeated experiments for Case 2
###########################################################
set.seed(8217)
n_repeat = 50
mse2 = matrix(rep(0,n_repeat*9),nrow=n_repeat)
res_bs = rep(0,n_repeat)
for(k in 1:n_repeat){
  tmp = mvrnorm(n,rep(0,3),S)
  C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
  X = pnorm((V+abs(W))/sqrt(2))
  Z = pnorm(abs(W))-0.5
  Y = 0.75*C - 0.25*Z + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
  
  xx = seq(0,1,length.out=n+2)[-c(1,n+2)]
  ty = log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)
  
  s = dis()
  sx = s[1]; sz = s[2]
  
  ###########################################################
  # Kernel
  ###########################################################
  gamma = 2
  alpha1 = 1.5/sqrt(n1)
  alpha2 = 1.5/sqrt(n2)
  alpha = 1.5/sqrt(n12)
  xi = 1.5/sqrt(m)
  
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
  mse2[k,1] = mean((ty-y_KAR)^2)
  
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
  mse2[k,2] = mean((ty-y_KAR2)^2)
  
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
  mse2[k,3] = mean((ty-y_KIV)^2)
  
  # KPA
  aY = t(Y[id3]) + (-1)*t(Y[id2])%*%solve(K(Z[id2],Z[id2],sz)+n2*alpha2*diag(n2))%*%K(Z[id2],Z[id3],sz)
  K1ZZ = K(Z[id1],Z[id1],sz); K13ZZ = K(Z[id1],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id1],sx); K1XX = K(X[id1],X[id1],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id1],xx,sx)
  KaXX = K3XX + 
    (-1)^2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K1XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX) +
    (-1)*2*K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%t(K31XX)
  Ka4 = K34XX + K31XX%*%solve(K1ZZ+n1*alpha1*diag(n1))%*%K14XX
  y_KPA = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
  y_KPA = t(y_KPA)
  # Error
  mse2[k,4] = mean((ty-y_KPA)^2)
  
  #KernelReg
  KXX = K(X,X,sx)
  KX4 = K(X,xx,sx)
  y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
  y_KReg = t(y_KReg)
  # Error
  mse2[k,5] = mean((ty-y_KReg)^2)
  
  ###########################################################
  # Linear estimate
  ###########################################################
  # AR
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
  b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_AR = c(b_AR)*(xx-0.5)
  mse2[k,6] = mean((ty-y_AR)^2)
  
  # IV
  aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_IV = c(b_IV)*(xx-mean(X))+mean(Y)
  mse2[k,7] = mean((ty-y_IV)^2)
  
  # PA
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_PA = c(b_PA)*(xx-mean(X))+mean(Y)
  mse2[k,8] = mean((ty-y_PA)^2)
  
  # Basic
  b = solve(t(X)%*%X)%*%t(aX)%*%Y
  y_b = c(b)*(xx-mean(X))+mean(Y)
  mse2[k,9] = mean((ty-y_b)^2)
  
  # CAM
  bs = CAM(cbind(Z,X,Y))
  res_bs[k] = bs$Adj[2,3]
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lmse2 = log(mse2,base=10)
colMeans(lmse2)
mean(res_bs)

###########################################################
# Visualization
###########################################################
res3 = data.frame('MSE'=c(lmse1[,1],lmse1[,2],lmse1[,3],lmse1[,4],
                         lmse1[,5],lmse1[,6],lmse1[,7],lmse1[,8],lmse1[,9]),
                 'Method'=factor(rep(c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS'),each=n_repeat),
                                 levels = c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS')),
                 'Case'=rep('Case 1',n_repeat*ncol(lmse1)))
res4 = data.frame('MSE'=c(lmse2[,1],lmse2[,2],lmse2[,3],lmse2[,4],
                          lmse2[,5],lmse2[,6],lmse2[,7],lmse2[,8],lmse2[,9]),
                  'Method'=factor(rep(c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS'),each=n_repeat),
                                  levels = c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS')),
                  'Case'=rep('Case 2',n_repeat*ncol(lmse2)))

g4 = ggplot(res3,aes(y=MSE))+
  ylab("MSE (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6],
                              'AR'=Mycol[7],'IV'=Mycol[8],'PA'=Mycol[9],'OLS'=Mycol[10]))+
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1),limits=c(-2,1.5))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
g5 = ggplot(res4,aes(y=MSE))+
  ylab("MSE (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6],
                              'AR'=Mycol[7],'IV'=Mycol[8],'PA'=Mycol[9],'OLS'=Mycol[10]))+
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1),limits=c(-2,1.5))+
  theme_classic()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

resc = rbind(res3,res4)
g6 = ggplot(resc,aes(x=Case,y=MSE))+
  ylab("MSE (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6],
                              'AR'=Mycol[7],'IV'=Mycol[8],'PA'=Mycol[9],'OLS'=Mycol[10]))+
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1),limits=c(-2,1.2))+
  theme_classic()+
  scale_x_discrete(limits=c('Case 1','Case 2'))+
  theme(axis.title.x=element_blank())
g6

###########################################################
# Shift between training and testing set
###########################################################
# Train: Z>=0 Test: Z<0
set.seed(8217)
n_repeat = 50
pre1 = matrix(rep(0,n_repeat*9),nrow=n_repeat)

for(k in 1:n_repeat){
  tmp = mvrnorm(n,rep(0,3),S)
  C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
  X = pnorm((V+W)/sqrt(2))
  Z = pnorm(W)-0.5
  for(i in 1:n){
    while(Z[i]<0){
      tmp = mvrnorm(1,rep(0,3),S)
      C[i] = tmp[1]; V[i] = tmp[2]; W[i] = tmp[3]
      X[i] = pnorm((V[i]+W[i])/sqrt(2))
      Z[i] = pnorm(W[i])-0.5
    }
  }
  Y = 0.75*C - 0.25*Z + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
  
  tmp = mvrnorm(n,rep(0,3),S)
  cc = tmp[,1]; vv = tmp[,2]; ww = tmp[,3]
  xx = pnorm((vv+ww)/sqrt(2))
  zz = pnorm(ww)-0.5
  for(i in 1:n){
    while(zz[i]>=0){
      tmp = mvrnorm(1,rep(0,3),S)
      cc[i] = tmp[1]; vv[i] = tmp[2]; ww[i] = tmp[3]
      xx[i] = pnorm((vv[i]+ww[i])/sqrt(2))
      zz[i] = pnorm(ww[i])-0.5
    }
  }
  yy = 0.75*cc - 0.25*zz + log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)
  
  s = dis()
  sx = s[1]; sz = s[2]
  
  ###########################################################
  # Kernel
  ###########################################################
  gamma = 2
  alpha1 = 1.5/sqrt(n1)
  alpha2 = 1.5/sqrt(n2)
  alpha = 1.5/sqrt(n12)
  xi = 1.5/sqrt(m)
  
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
  pre1[k,1] = mean((yy-y_KAR)^2)
  
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
  pre1[k,2] = mean((yy-y_KAR2)^2)
  
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
  pre1[k,3] = mean((yy-y_KIV)^2)
  
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
  pre1[k,4] = mean((yy-y_KPA)^2)
  
  #KernelReg
  KXX = K(X,X,sx)
  KX4 = K(X,xx,sx)
  y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
  y_KReg = t(y_KReg)
  # Error
  pre1[k,5] = mean((yy-y_KReg)^2)
  
  ###########################################################
  # Linear estimate
  ###########################################################
  # AR
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
  b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_AR = c(b_AR)*(xx-mean(X))+mean(Y)
  pre1[k,6] = mean((yy-y_AR)^2)
  
  # IV
  aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_IV = c(b_IV)*(xx-mean(X))+mean(Y)
  pre1[k,7] = mean((yy-y_IV)^2)
  
  # PA
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_PA = c(b_PA)*(xx-mean(X))+mean(Y)
  pre1[k,8] = mean((yy-y_PA)^2)
  
  # Basic
  b = solve(t(X)%*%X)%*%t(aX)%*%Y
  y_b = c(b)*(xx-mean(X))+mean(Y)
  pre1[k,9] = mean((ty-y_b)^2)
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lpre1 = log(pre1,base=10)
colMeans(lpre1)

# Train: Z<0 Test: Z>=0
set.seed(8217)
n_repeat = 50
pre2 = matrix(rep(0,n_repeat*9),nrow=n_repeat)
for(k in 1:n_repeat){
  tmp = mvrnorm(n,rep(0,3),S)
  C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
  X = pnorm((V+W)/sqrt(2))
  Z = pnorm(W)-0.5
  for(i in 1:n){
    while(Z[i]>=0){
      tmp = mvrnorm(1,rep(0,3),S)
      C[i] = tmp[1]; V[i] = tmp[2]; W[i] = tmp[3]
      X[i] = pnorm((V[i]+W[i])/sqrt(2))
      Z[i] = pnorm(W[i])-0.5
    }
  }
  Y = 0.75*C - 0.25*Z + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
  
  tmp = mvrnorm(n,rep(0,3),S)
  cc = tmp[,1]; vv = tmp[,2]; ww = tmp[,3]
  xx = pnorm((vv+ww)/sqrt(2))
  zz = pnorm(ww)-0.5
  for(i in 1:n){
    while(zz[i]<0){
      tmp = mvrnorm(1,rep(0,3),S)
      cc[i] = tmp[1]; vv[i] = tmp[2]; ww[i] = tmp[3]
      xx[i] = pnorm((vv[i]+ww[i])/sqrt(2))
      zz[i] = pnorm(ww[i])-0.5
    }
  }
  yy = 0.75*cc - 0.25*zz + log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)
  
  s = dis()
  sx = s[1]; sz = s[2]
  
  ###########################################################
  # Kernel
  ###########################################################
  gamma = 2
  alpha1 = 1.5/sqrt(n1)
  alpha2 = 1.5/sqrt(n2)
  alpha = 1.5/sqrt(n12)
  xi = 1.5/sqrt(m)
  
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
  pre2[k,1] = mean((yy-y_KAR)^2)
  
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
  pre2[k,2] = mean((yy-y_KAR2)^2)
  
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
  pre2[k,3] = mean((yy-y_KIV)^2)
  
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
  pre2[k,4] = mean((yy-y_KPA)^2)
  
  #KernelReg
  KXX = K(X,X,sx)
  KX4 = K(X,xx,sx)
  y_KReg = Y%*%solve(KXX+n*xi*diag(n))%*%KX4
  y_KReg = t(y_KReg)
  # Error
  pre2[k,5] = mean((yy-y_KReg)^2)
  
  ###########################################################
  # Linear estimate
  ###########################################################
  # AR
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X + sqrt(gamma)*Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X            
  b_AR = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_AR = c(b_AR)*(xx-mean(X))+mean(Y)
  pre2[k,6] = mean((yy-y_AR)^2)
  
  # IV
  aY = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_IV = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_IV = c(b_IV)*(xx-mean(X))+mean(Y)
  pre2[k,7] = mean((yy-y_IV)^2)
  
  # PA
  aY = Y - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%Y
  aX = X - Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X
  b_PA = solve(t(aX)%*%aX)%*%t(aX)%*%aY
  y_PA = c(b_PA)*(xx-mean(X))+mean(Y)
  pre2[k,8] = mean((yy-y_PA)^2)
  
  # Basic
  b = solve(t(X)%*%X)%*%t(aX)%*%Y
  y_b = c(b)*(xx-mean(X))+mean(Y)
  pre2[k,9] = mean((ty-y_b)^2)
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lpre2 = log(pre2,base=10)
colMeans(lpre2)

###########################################################
# Visualization
###########################################################
res5 = data.frame('MSE'=c(lpre1[,1],lpre1[,2],lpre1[,3],lpre1[,4],
                          lpre1[,5],lpre1[,6],lpre1[,7],lpre1[,8],lpre1[,9]),
                  'Method'=factor(rep(c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS'),each=n_repeat),
                                  levels = c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS')),
                  'Case'=rep('Train: Z<0 Test: Z>=0',9*n_repeat))

res6 = data.frame('MSE'=c(lpre2[,1],lpre2[,2],lpre2[,3],lpre2[,4],
                          lpre2[,5],lpre2[,6],lpre2[,7],lpre2[,8],lpre2[,9]),
                  'Method'=factor(rep(c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS'),each=n_repeat),
                                  levels = c('KAR','KAR.2','KIV','KPA','KReg','AR','IV','PA','OLS')),
                  'Case'=rep('Train: Z>=0 Test: Z<0',9*n_repeat))
resc2 = rbind(res5,res6)
gc2 = ggplot(resc2,aes(x=Case,y=MSE))+
  ylab("Prediction Error (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('KAR'=Mycol[2],'KAR.2'=Mycol[3],'KIV'=Mycol[4],'KPA'=Mycol[5],'KReg'=Mycol[6],
                              'AR'=Mycol[7],'IV'=Mycol[8],'PA'=Mycol[9],'OLS'=Mycol[10]))+
  theme_classic()+
  scale_x_discrete(limits=c('Train: Z<0 Test: Z>=0','Train: Z>=0 Test: Z<0'))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=12))


###########################################################
# Kernel IV setting
###########################################################
# Select parameters for alphas and xi
set.seed(123)
xi = 1/sqrt(m)
multig = c(0,0.5,1,2,5,10,100)
multic = c(0.01,0.05,0.1,0.5,0.8,1,2,3)
cKAR = rep(0,length(multig))
cKAR.2 = rep(0,length(multig))
cKIV = 0

tmp = mvrnorm(n,rep(0,3),S)
C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
X = pnorm((V+W)/sqrt(2))
Z = pnorm(W)
Y = C + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
xx = seq(0,1,length.out=n+2)[-c(1,n+2)]
ty = log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)

s = dis()
sx = s[1]; sz = s[2]
for(g in c(1:length(multig))){
  gamma = multig[g]
  tmp = rep(0,length(multic))
  for(c in c(1:length(multic))){
    alpha1 = multic[c]/sqrt(n1)
    alpha2 = multic[c]/sqrt(n2)
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
    tmp[c] = mean((ty-y_KAR)^2)
  }
  cKAR[g] = multic[which.min(tmp)]
  print(c(multig[g],cKAR[g],log(min(tmp),base=10)))
}
for(g in c(1:length(multig))){
  gamma = multig[g]
  tmp = rep(0,length(multic))
  for(c in c(1:length(multic))){
    alpha = multic[c]/sqrt(n12)
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
    tmp[c] = mean((ty-y_KAR2)^2)
  }
  cKAR.2[g] = multic[which.min(tmp)]
  print(c(multig[g],cKAR.2[g],log(min(tmp),base=10)))
}
tmp = rep(0,length(multic))
for(c in c(1:length(multic))){
  alpha = multic[c]/sqrt(n12)
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  W = K1XX%*%solve(K1ZZ + n12*alpha*diag(n12))%*%K13ZZ
  y_KIV = t(ginv(W%*%t(W)+m*xi*K1XX)%*%W%*%Y[id3])%*%K14XX
  y_KIV = t(y_KIV)
  # Error
  tmp[c] = mean((ty-y_KIV)^2)
  print(c(multic[c],log(tmp[c],base=10)))
} 
print(c(log(min(tmp),base=10)))
cKIV = multic[which.min(tmp)]

# Repeated experiments for KIV setting
set.seed(321)
mse3 = matrix(rep(0,n_repeat*(length(multig)*2+1)),nrow=n_repeat)
for(k in 1:n_repeat){
  tmp = mvrnorm(n,rep(0,3),S)
  C = tmp[,1]; V = tmp[,2]; W = tmp[,3]
  X = pnorm((V+W)/sqrt(2))
  Z = pnorm(W)
  Y = C + log(abs(16*(X)-8)+1)*(X-0.5)/abs(X-0.5)
  xx = seq(0,1,length.out=n+2)[-c(1,n+2)]
  ty = log(abs(16*(xx)-8)+1)*(xx-0.5)/abs(xx-0.5)
  s = dis()
  sx = s[1]; sz = s[2]
  
  # KAR
  for(g in c(1:length(multig))){
    gamma = multig[g]
    alpha1 = cKAR[g]/sqrt(n1)
    alpha2 = cKAR[g]/sqrt(n2)
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
    mse3[k,g] = mean((ty-y_KAR)^2)
  }
  
  # KAR.2
  for(g in c(1:length(multig))){
    gamma = multig[g]
    alpha = cKAR.2[g]/sqrt(n12)
    iK1ZZ = solve(K(Z[id12],Z[id12],sz)+n12*alpha*diag(n12))
    aY = t(Y[id3]) + (sqrt(gamma)-1)*t(Y[id12])%*%iK1ZZ%*%K(Z[id12],Z[id3],sz)
    K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
    K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
    K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
    KaXX = K3XX + 
      (sqrt(gamma)-1)^2*K31XX%*%iK1ZZ%*%K1XX%*%iK1ZZ%*%t(K31XX) +
      (sqrt(gamma)-1)*2*K31XX%*%iK1ZZ%*%t(K31XX)
    Ka4 = K34XX + K31XX%*%iK1ZZ%*%K14XX
    y_KAR2 = aY%*%solve(KaXX+m*xi*diag(m))%*%Ka4
    y_KAR2 = t(y_KAR2)
    # Error
    mse3[k,(g+length(multig))] = mean((ty-y_KAR2)^2)
  }
  
  # KIV
  alpha = cKIV/sqrt(n12)
  K1ZZ = K(Z[id12],Z[id12],sz); K13ZZ = K(Z[id12],Z[id3],sz)
  K3XX = K(X[id3],X[id3],sx); K31XX = K(X[id3],X[id12],sx); K1XX = K(X[id12],X[id12],sx)
  K34XX = K(X[id3],xx,sx); K14XX = K(X[id12],xx,sx)
  W = K1XX%*%solve(K1ZZ + n12*alpha*diag(n12))%*%K13ZZ
  y_KIV = t(ginv(W%*%t(W)+m*xi*K1XX)%*%W%*%Y[id3])%*%K14XX
  y_KIV = t(y_KIV)
  # Error
  mse3[k,(length(multig)*2+1)] = mean((ty-y_KIV)^2) 
  
  if(k%%10==0) print(paste(k/n_repeat*100,"%"))
}
lmse3 = log(mse3,base=10)
colMeans(lmse3)

###########################################################
# Visualization
###########################################################
res = data.frame('MSE'=c(lmse3[,1],lmse3[,2],lmse3[,3],lmse3[,4],
                         lmse3[,5],lmse3[,6],lmse3[,7],lmse3[,8],
                         lmse3[,9],lmse3[,10],lmse3[,11],lmse3[,12],
                         lmse3[,13],lmse3[,14],lmse3[,15]),
                 'Method'=factor(rep(c('gamma=0','gamma=0.5','gamma=1','gamma=2','gamma=5','gamma=10','gamma=100',
                                       'gamma=0','gamma=0.5','gamma=1','gamma=2','gamma=5','gamma=10','gamma=100','KIV'),each=n_repeat),
                                 levels = c('gamma=0','gamma=0.5','gamma=1','gamma=2','gamma=5','gamma=10','gamma=100','KIV')),
                 'Method2'=factor(c(rep("KAR",n_repeat*length(multig)), rep("KAR.2",n_repeat*length(multig)),
                                    rep("KIV",n_repeat))))
g = ggplot(res,aes(x=Method2,y=MSE))+
  ylab("MSE (log)")+
  geom_boxplot(aes(color=Method))+
  scale_color_manual(values=c('gamma=0'=Mycol[2],'gamma=0.5'=Mycol[3],'gamma=1'=Mycol[4],'gamma=2'=Mycol[5],
                              'gamma=5'=Mycol[6],'gamma=10'=Mycol[7],'gamma=100'=Mycol[8],'KIV'=Mycol[9]))+
  scale_y_continuous(breaks=c(-1.5,-1,-0.5,0,0.5,1))+
  scale_x_discrete(limits=c('KAR','KAR.2',"KIV"))+
  theme_classic()+
  theme(axis.title.x=element_blank())
g
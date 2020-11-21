## Code includes 6 methods, 
## 4 from first-order table, they are:
##  **Frank-Wolfe
##  **GD
##  **SGD
##  **AdaGrad
## and 2 from high-order table:
##  **Newton Method
##  **CG


library(tidyverse)
library(pracma)
library(ggplot2)

n=400
p=2
beta=matrix(c(1,-2),ncol = 1)
#data
y=rbinom(n,1,0.5)
x=cbind(rnorm(n),rnorm(n))

#Likelihood Function
Likeli=function(b){
  u=exp(x%*%b)/(1+exp(x%*%b))
  -sum(y*log(u/(1-u))+log(1-u))
}

#Training Error evaluator
train_err=function(b){
  sum(y*(x%*%b)<0)/n#y=1,hwr pre=0
}

#solution with optim() in R
optim(beta,Likeli)


#############***Agorithm******##########

#####----Frank Wolfe-----#####
bt=beta
g=grad(Likeli,bt)

betamesh=matrix(c(-1,-1,-1,1,1,1,1,-1),nrow  = 2)

argmin=function(g){
  fvalue=g%*%betamesh
  min_b=betamesh[,which(fvalue==min(fvalue))]
}

error=c()
for(i in 1:700){
  error[i]=train_err(bt)
  min_bt=argmin(g)
  bt=i/(i+2)*bt+2/(i+2)*(min_bt-bt)
  g=grad(Likeli,bt)
}
bt1=bt
plot(error,type = "l",col='#9400D3', xlab = "Iteration Times",ylab = "Training Error")


#########-----GD-----#########
bt=beta
g=-grad(Likeli,bt)
error=c()
for(i in 1:15){
  bt=bt+0.01*g
  g=-grad(Likeli,bt)
  error[i]=train_err(bt)
}

plot(error,type = "l",lwd=2,col="DarkTurquoise",ylim = c(0.2,0.24),lty=1,
     xlab = "Iteration Times",ylab = "Training Error")
points(error,pch=15,col="DarkTurquoise")
bt2=bt


#########----SGD-----#########
rdm_idx=sample(1:400,1)
x1=x[rdm_idx,]
y1=y[rdm_idx]

Likeli1=function(b){
  u=exp(x1%*%b)/(1+exp(x1%*%b))
  -sum(y1*log(u/(1-u))+log(1-u))
}

bt=beta
g=-grad(Likeli1,bt)

error=c()
for(i in 1:10){
  error[i]=train_err(bt)
  bt=bt+0.01*g
  g=-grad(Likeli,bt)
  
}
bt
points(error,pch=16,col="DeepPink")
lines(error,lwd=2,col= "DeepPink",lty=2)
bt3=bt

#####-------AdaGrad-------#######
bt=beta
g=-grad(Likeli,bt)
ita=0
e=1e-7
error=c()
for(i in 1:15){
  error[i]=train_err(bt)
  ita=ita+sum(g^2)
  bt=bt+1/sqrt(ita+e)*g
  g=-grad(Likeli,bt)
}
bt
points(error,pch=17,col="RosyBrown")
lines(error,lwd=2,col= "RosyBrown",lty=3)
bt4=bt

#########---Newton-----########
bt=c(1,-1)
error=c()
for (i in 1:15){
  error[i]=train_err(bt)
  g=grad(Likeli,bt)
  dir=solve(hessian(Likeli,bt))%*%g
  bt=bt-dir
}
bt

points(error,pch=18,col="SeaGreen")
lines(error,lwd=2,col= "SeaGreen",lty=4)
bt5=bt

########------CG------##########
bt=beta
g=grad(Likeli,bt)
d=(-g)%>%as.matrix()#initial d1
#g>=err
iter=1
error=c()
while(iter<10){
  error[iter]=train_err(bt)
  fvalue=c()
  sq=seq(0,0.05,0.00001)
  for (i in sq){
    fvalue[i*100000+1]=Likeli(bt+i*d)
  }
  temp=sq[which(fvalue==min(fvalue))]
  bt=bt+temp*d
  temp_g2=sum(g^2)
  g=grad(Likeli,bt)
  d=-g+(sum(g^2)/temp_g2)*d
  iter=iter+1;
}
iter
bt6=bt

points(error,pch=19,col="DarkMagenta")
lines(error,lwd=2,col= "DarkMagenta",lty=5)

legend(10,0.24,c("GD","SGD","AdaGrad","Newton","FR"),col=c("DarkTurquoise","DeepPink","RosyBrown","SeaGreen","DarkMagenta"),
       text.col=c("DarkTurquoise","DeepPink","RosyBrown","SeaGreen","DarkMagenta"),pch=c(15,16,17,18,19),lty=c(1,2,3,4,5))

bt_table=cbind(optim(beta,Likeli)[[1]],bt1,bt2,bt3,bt4,bt5,bt6)
colnames(bt_table)=c("optim in R","Frank-Wolfe","GD","SGD","AdaGrad","Newton's","FR")
bt_table

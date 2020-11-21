library(tidyverse)
n=400
p=2
beta=matrix(c(1,-2),ncol = 1)
#data
y=rbinom(n,1,0.5)
x=cbind(1,rnorm(n))


Likeli=function(b){
  u=exp(x%*%b)/(1+exp(x%*%b))
  -sum(y*log(u/(1-u))+log(1-u))
}

library(pracma)

err=1e-3#error limit
bt=beta
g=grad(Likeli,bt)

# bt1=seq(0,1,0.001)
# bt2=seq(0,1,0.001)
# cbind(rep(0,1001),bt2)
# betamesh=cbind(rep(bt1,1001),rep(bt2,each=1001))%>%t
betamesh=matrix(c(-1,-1,-1,1,1,1,1,-1),nrow  = 2)

argmin=function(g){
  fvalue=g%*%betamesh
  min_b=betamesh[,which(fvalue==min(fvalue))]
}

bt=beta
for(i in 1:500){
  min_bt=argmin(g)
  bt=bt+2/(i+2)*(min_bt-bt)
  g=grad(Likeli,bt)
}
bt
optim(beta,Likeli)

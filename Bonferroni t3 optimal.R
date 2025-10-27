####Bonferroni t3####
####Optimal #######


alpha=.05
beta=.1
v=3
mu1=1
mu2=1.5
mu3=1.5
mu4=2.5
sigma1=.75
sigma2=1
sigma3=1
sigma4=.75

h=2


delta0=.5
delta1=1.5

delta=delta1-delta0

n0=5
k=(delta0*qnorm(1-beta/3)+delta1*qnorm(1-alpha/3))/(qnorm(1-alpha/3)+qnorm(1-beta/3))

iterations = 500
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-1)
registerDoParallel(cl)

result2=foreach (i=1:iterations,.combine=rbind,.packages = c("smoothmest","NlcOptim","MASS"))%dopar%{  
  set.seed(i+2000)
  
  x1=sigma1*rt(n0,v)+mu1
  x4=sigma4*rt(n0,v)+mu4
  x2=sigma2*rt(n0,v)+mu2
  x3=sigma3*rt(n0,v)+mu3
  stopcrmet=FALSE
  while(!stopcrmet){
    obj=function(z){
      return(length(x1)*log(z[1])+((v+1)/2)*sum(log(1+((x1-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x10=c(sigma1,mu1)
    sigma1hat=solnl(x10,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu1hat=solnl(x10,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x2)*log(z[1])+((v+1)/2)*sum(log(1+((x2-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x20=c(sigma2,mu2)
    sigma2hat=solnl(x20,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu2hat=solnl(x20,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x3)*log(z[1])+((v+1)/2)*sum(log(1+((x3-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x30=c(sigma3,mu3)
    sigma3hat=solnl(x30,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu3hat=solnl(x30,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x4)*log(z[1])+((v+1)/2)*sum(log(1+((x4-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x40=c(sigma4,mu4)
    sigma4hat=solnl(x40,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu4hat=solnl(x40,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    
    
    var1=v*sigma1hat^2/(v-2)
    var2=v*sigma2hat^2/(v-2)
    var3=v*sigma3hat^2/(v-2)
    var4=v*sigma4hat^2/(v-2)
    
    psi1hat=pt((h-mu1hat)/sigma1hat,df=v)
    psi2hat=pt((h-mu2hat)/sigma2hat,df=v)
    psi3hat=pt((h-mu3hat)/sigma3hat,df=v)
    psi4hat=pt((h-mu4hat)/sigma4hat,df=v)
    
    
    pi1=sqrt(var1/psi1hat)*
      sqrt(var2*psi2hat+
             var3*psi3hat+
             var4*psi4hat)
    pi2=var2
    pi3=var3
    pi4=var4
    
    pi11=pi1/(pi1+pi2+pi3+pi4)
    pi22=pi2/(pi1+pi2+pi3+pi4)
    pi33=pi3/(pi1+pi2+pi3+pi4)
    pi44=pi4/(pi1+pi2+pi3+pi4)
    p=c(pi11,pi22,pi33,pi44)
    u=sample(c(1,2,3,4),1,prob=p)
    if(u==1){
      x1=c(x1,sigma1*rt(1,v)+mu1)
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,sigma2*rt(1,v)+mu2)
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,sigma3*rt(1,v)+mu3)
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,sigma4*rt(1,v)+mu4)
    }
    obj=function(z){
      return(length(x1)*log(z[1])+((v+1)/2)*sum(log(1+((x1-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x10=c(sigma1,mu1)
    sigma1hat=solnl(x10,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu1hat=solnl(x10,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x2)*log(z[1])+((v+1)/2)*sum(log(1+((x2-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x20=c(sigma2,mu2)
    sigma2hat=solnl(x20,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu2hat=solnl(x20,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x3)*log(z[1])+((v+1)/2)*sum(log(1+((x3-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x30=c(sigma3,mu3)
    sigma3hat=solnl(x30,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu3hat=solnl(x30,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    obj=function(z){
      return(length(x4)*log(z[1])+((v+1)/2)*sum(log(1+((x4-z[2])^2/(v*z[1]^2)))))
    }
    
    budgetconstraint<-function(z){
      f=NULL
      f=rbind(f,.5-z[1])
      return(list(ceq = NULL, c = f))
    }
    x40=c(sigma4,mu4)
    sigma4hat=solnl(x40,obj,budgetconstraint,tolX = 1e-5,
                    tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[1]
    mu4hat=solnl(x40,obj,budgetconstraint,tolX = 1e-5,
                 tolFun = 1e-5, tolCon = 1e-5, maxnFun = 1e+10,maxIter = 500)$par[2]
    
    var1=v*sigma1hat^2/(v-2)
    var2=v*sigma2hat^2/(v-2)
    var3=v*sigma3hat^2/(v-2)
    var4=v*sigma4hat^2/(v-2)
    
    if(var1/length(x1)+(var2/length(x2))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2 &&
       (var1/length(x1))+(var3/length(x3))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2  &&
       (var1/length(x1))+(var4/length(x4))<=(delta^2)/(qnorm(1-alpha/3)+qnorm(1-beta/3))^2){
      stopcrmet=TRUE
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)
  propfail=(sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h))/N
  NF=sum(x1<h)+sum(x2<h)+sum(x3<h)+sum(x4<h)
  
  if(mean(x2)-mean(x1)>k){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>k){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>k){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  
  c(N,propfail,NF,max(r2,r3),max(a4))
}

colMeans(result1)
sd(result1[,1])
sd(result1[,2])


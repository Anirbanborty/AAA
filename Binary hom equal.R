###########binary holm equal#####





stopCluster(cl)
alpha=.05
beta=.1
mu1=.3
mu2=.4
mu3=.4
mu4=.6
sigma1=sqrt(mu1*(1-mu1))
sigma2=sqrt(mu2*(1-mu2))
sigma3=sqrt(mu3*(1-mu3))
sigma4=sqrt(mu4*(1-mu4))



delta0=.1
delta1=.3

delta=delta1-delta0


n0=5


iterations = 10000
library(parallel)
library(doParallel)
library(foreach)
nc=detectCores()
cl=makeCluster(nc-2)
registerDoParallel(cl)

result=foreach (i=1:iterations,.combine=rbind,.packages = c("NlcOptim","MASS"))%dopar%{  
  
  x1=c(1,0,rbinom(5,1,mu1))
  x4=c(1,0,rbinom(5,1,mu4))
  x2=c(1,0,rbinom(5,1,mu2))
  x3=c(1,0,rbinom(5,1,mu3))
  stopcrmet=FALSE
  while(!stopcrmet){
    pi=c(1/4,1/4,1/4,1/4)
    u=sample(c(1,2,3,4),1,prob=pi)
    if(u==1){
      x1=c(x1,rbinom(1,1,mu1))
      x2=x2
      x3=x3
      x4=x4
    }else if(u==2){
      x1=x1
      x2=c(x2,rbinom(1,1,mu2))
      x3=x3
      x4=x4
    } else if (u==3) {
      x2=x2
      x1=x1
      x3=c(x3,rbinom(1,1,mu3))
      x4=x4
    } else {
      x2=x2
      x1=x1
      x3=x3
      x4=c(x4,rbinom(1,1,mu4))
      
    }
    
    d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
    if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
      d2=which(d==mean(x2))
      d3=which(d==mean(x3))
      d4=which(d==mean(x4))
    } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
      d3=1
      u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
      d4=u[1]
      d2=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)>mean(x3)){
      d4=1
      u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
      d2=1
      u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
      d4=u[1]
      d3=u[2]
    }else if (mean(x3)==mean(x2) && mean(x4)<mean(x3)){
      d4=3
      u=rank(c(mean(x3),mean(x2)),ties.method= "random")
      d2=u[1]
      d3=u[2]
    }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
      d3=3
      u=rank(c(mean(x4),mean(x2)),ties.method= "random")
      d2=u[1]
      d4=u[2]
    }else if (mean(x3)==mean(x4) && mean(x4)>mean(x2)){
      d2=3
      u=rank(c(mean(x3),mean(x4)),ties.method= "random")
      d4=u[1]
      d3=u[2]
    }else{
      u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
      d4=u[1]
      d2=u[2]
      d3=u[3]
    }
    
    if((mean(x1)*(1-mean(x1))/length(x1))+(mean(x2)*(1-mean(x2))/length(x2))<=delta^2/(qnorm(1-alpha/(3-d2+1))+
                                                                                       qnorm(1-beta/d2))^2      &&
       (mean(x1)*(1-mean(x1))/length(x1))+(mean(x3)*(1-mean(x3))/length(x3))<=delta^2/(qnorm(1-alpha/(3-d3+1))+
                                                                                       qnorm(1-beta/d3))^2   &&
       (mean(x1)*(1-mean(x1))/length(x1))+(mean(x4)*(1-mean(x4))/length(x4))<=delta^2/(qnorm(1-alpha/(3-d4+1))+
                                                                                       qnorm(1-beta/d4))^2){
      stopcrmet=TRUE
    }
    
  }
  N=length(x1)+length(x2)+length(x3)+length(x4)
  propfail=(sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0))/N
  NF= sum(x1==0)+sum(x2==0)+sum(x3==0)+sum(x4==0)
  d=sort(c(mean(x2),mean(x3),mean(x4)),decreasing = TRUE)
  if(mean(x4)!=mean(x2) && mean(x2)!=mean(x3) && mean(x4)!=mean(x3)){
    d2=which(d==mean(x2))
    d3=which(d==mean(x3))
    d4=which(d==mean(x4))
  } else if (mean(x4)==mean(x2) && mean(x3)>mean(x4)){
    d3=1
    u=1+rank(c(mean(x4),mean(x2)),ties.method= "random")
    d4=u[1]
    d2=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)>mean(x3)){
    d4=1
    u=1+rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x4) && mean(x2)>mean(x3)){
    d2=1
    u=1+rank(c(mean(x3),mean(x4)),ties.method= "random")
    d4=u[1]
    d3=u[2]
  }else if (mean(x3)==mean(x2) && mean(x4)<mean(x3)){
    d4=3
    u=rank(c(mean(x3),mean(x2)),ties.method= "random")
    d2=u[1]
    d3=u[2]
  }else if (mean(x4)==mean(x2) && mean(x4)>mean(x3)){
    d3=3
    u=rank(c(mean(x4),mean(x2)),ties.method= "random")
    d2=u[1]
    d4=u[2]
  }else if (mean(x3)==mean(x4) && mean(x4)>mean(x2)){
    d2=3
    u=rank(c(mean(x3),mean(x4)),ties.method= "random")
    d4=u[1]
    d3=u[2]
  }else{
    u=rank(c(mean(x3),mean(x2),mean(x1)),ties.method= "random")
    d4=u[1]
    d2=u[2]
    d3=u[3]
  }
  
  if(mean(x2)-mean(x1)>(delta0*qnorm(1-beta/d2)+
                        delta1*qnorm(1-alpha/(3-d2+1)))/
     (qnorm(1-alpha/(3-d2+1))+
      qnorm(1-beta/d2))){
    r2=1
    a2=0
  }else{
    r2=0
    a2=1
  }
  if(mean(x3)-mean(x1)>(delta0*qnorm(1-beta/d3)+
                        delta1*qnorm(1-alpha/(3-d3+1)))/
     (qnorm(1-alpha/(3-d3+1))+
      qnorm(1-beta/d3))){
    r3=1
    a3=0
  }else{
    r3=0
    a3=1
  }
  if(mean(x4)-mean(x1)>(delta0*qnorm(1-beta/d4)+
                        delta1*qnorm(1-alpha/(3-d4+1)))/
     (qnorm(1-alpha/(3-d4+1))+
      qnorm(1-beta/d4))){
    r4=1
    a4=0
  }else{
    r4=0
    a4=1
  }
  
  c(N,propfail,NF,max(r2,r3),max(a4))
}

colMeans(result)
colMeans(result)
sd(result[,1])
sd(result[,2])



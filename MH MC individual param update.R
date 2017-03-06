
##case2-- miss prob that decays with size...
#--p(s|gamma) = gamma[2]+(1-gamma[2])*exp(-(s-gamma[1])/gamma[3]) when s>=gamma[1]
#--o.w. p=1
# assume limit when s->inf =0...
# pmiss=exp(-(s-gamma[1])/gamma[2]) when s>=gamma[1]

##case1--miss prob which is a step function...
#--p(s|gamma) = gamma[2] when s>=gamma[1] 
#--o.w. p=1

#---------------------------------------
J=13
# pretreated data from a rectangle, assume A=1
cra <- read.table("/Users/zhang/Desktop/crater after pretreat.txt",header=T)
cra2 <- cra[cra$repeats!=1,]

N2=nrow(cra2); J=13; n=as.integer(cra2[,2]); s=cra2[,1]

s1=cra[cra[,2]==1,][,1]
dn1=pretreat[,3:4][pretreat[,3] %in% as.character(s1),]
table(dn1$observer)
levels(dn1$observer)
hist(sort(dn1$Diameter[obs==j[13]]),breaks=50) ## size in Nj for observer i

#----------------------------------------
Pmiss=function(s,gamma){
  p=NULL
  sx=s[s >= gamma[1]]
  p[s >= gamma[1]]=gamma[2]+(1-gamma[2])*exp(-(sx-gamma[1])/gamma[3]) 
  p[s < gamma[1]] =1
  return(p)
} # assume Pmiss is the same for all observers, same gamma

Pdf=function(s,beta){
  beta[1]*exp(-beta[1]*log(beta[2]))*s^(beta[1]-1)
} # 0< size <beta_2,  and 1>beta_1>0, beta_2>1


Ppha=function(s,eta){
  eta[1]*exp(-eta[1]*log(eta[2]))*s^(eta[1]-1)
} ## pdf for fake craters

#----------------------------
qI_s=function(n,s,gamma){
  (1-Pmiss(s,gamma))^n * Pmiss(s,gamma)^(J-n)
} # n= #seen ;the probability that a crater of size s is seen by n person

q2_s=function(s,gamma){
  1-qI_s(0,s,gamma)-qI_s(1,s,gamma)*(J-1)
} #  the probability that a crater of size s is seen by at least 2 analysts

qI=function(n,gamma,beta){
  f=function(s) qI_s(n,s,gamma)*Pdf(s,beta)
  int=integrate(f,lower = 0, upper = beta[2],rel.tol =10^-9)
  return(int$value)
}

q2=function(gamma,beta){
  f=function(s) q2_s(s,gamma)*Pdf(s,beta)
  int=integrate(f,lower = 0, upper = beta[2])
  return(int$value)
}

###---------------------------------------------------
#Log-likelihood
Logl=function(dat2,dat1,rho,gamma,beta,eta,n0,rhoj,Tj){
  s=dat2[,1]
  n=dat2[,2]
  l_02= -rho*q2(gamma,beta) + sum( log(qI_s(n,s,gamma)*Pdf(s,beta)*rho) ) -
    rho*qI(0,gamma,beta) + n0*log(rho*qI(0,gamma,beta))## will cancel out - log(factorial(n0))
  
  s1=dat1[,1]
  obs=dat1[,2]
  j=levels(obs)
  l_1=NULL
  for (i in 1:length(j)){
    if (sum(obs==j[i])==0) l_1[i]=0 else {
    si=s1[obs==j[i]] ## choose data for observer i
    ti=Tj[obs==j[i]]
    s_1=si[ti==1]
    s_0=si[ti==0]
    l_1[i]= -rho* qI(1,gamma,beta) - rhoj[i] + 
      sum( log(rho*qI_s(1,s_1,gamma)*Pdf(s_1,beta)) ) +
      sum( log(rhoj[i]*Ppha(s_0,eta)) )
    }
  }
  l=l_02+sum(l_1)
  return(l)
}##added rho

Logl(dat2=cra2,dat1=dn1,rho=700,gamma=c(5,0.5,5),
     beta=c(.8,3),eta=c(2,3),n0=50,rhoj=c(1:13),Tj=Tj)


##--------------------initial value for Tj
Tj=NULL
Tj[dn1$Diameter>gamma[1]]=rbinom(sum(dn1$Diameter>gamma[1]),1,0.5) 
Tj[dn1$Diameter<=gamma[1]]=0


#------------------------------------------------------------
##Gibbs sampling...

#------update Tj with bern(p_tj) conditonal on all other params.
newT=function(param,dat1=dn1){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21]
  
  p_tj=NULL
  s1=dat1[,1]
  obs=dat1[,2]
  j=levels(obs)
  for (i in 1:13){
    if (sum(obs==j[i])!=0) {
      sj=s1[obs==j[i]] ## select data for observer i
      p_tj[obs==j[i]]= rho*qI_s(1,sj,gamma)*Pdf(sj,beta)/( rho*qI_s(1,sj,gamma)*Pdf(sj,beta) + rhoj[i]*Ppha(sj,eta) )
    }
  }
  p_tj
  hist(p_tj)
  ## update Tj to Tj.1
  Tj.1=NULL
  for (i in 1:length(Tj)){
    Tj.1[i]=rbinom(1,1,p_tj[i])
  }
  return(Tj.1)
}

newT(param=xxxx[5,])
#---------------------------------------
##-----update n0 with Poisson(rho*qI(0,gamma,beta)) conditonal on all other params.
newN0=function(param){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21]
  
  n0.1=rpois(1,rho*qI(0,gamma,beta))
  return(n0.1)
}
newN0(param=xxxx[5,])

#---------------------------------------
##update other parameters using M-H algorithm. h = L(param)*g(param) 
## with dat2=cra2; dat1=dn1
## param order: rho"gamma1"gamma2"gamma3"beta1"beta2"eta1"eta2"rhoj[1:13]
#------------------------------------------
Logprior=function(param){
  rho=param[1]  ###param order!!!!!
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21]
  
  rho.p=dgamma(rho,36,6/100, log = T) ## intensity of craters
  gamma1.p=dunif(gamma[1],0,min(cra2$Diameter), log = T) ## criteria size for visible crater
  gamma2.p=dunif(gamma[2],0,1, log = T) ##limit prob for missing
  gamma3.p=dgamma(gamma[3],3, log = T)  ##decreasing rate for miss prob
  beta1.p=dunif(beta[1],0,1, log = T)  ## >0 according data <1
  beta2.p=dgamma(beta[2]-1,27*3,3/100, log = T) ## beta[2]>1, 0 < size < beta[2]  max(cra$Diameter)
  eta1.p=dunif(eta[1],0,5, log = T)
  eta2.p=dgamma(eta[2]-1,36,2/10, log = T) ## max(dn1$Diameter)
   mean.j=table(dn1$observer)+0.5
   sd.j=table(dn1$observer)/10+1
  rhoj.p=dgamma(rhoj,(mean.j^2/sd.j^2),(mean.j/sd.j^2), log = T) ## for vector rhoj[1:13]
  
  prior=rho.p+gamma1.p+gamma2.p+gamma3.p+beta1.p+beta2.p+eta1.p+eta2.p+sum(rhoj.p)
  names(prior)=c("logprior")
  return(prior)
}

init.para=c(rho=700,gamma=c(5,0.5,5),beta=c(.8,max(cra$Diameter)),
            eta=c(1.8,max(dn1$Diameter)),rhoj=table(dn1$observer)+1)
Logprior(init.para)

Loglike=function(param,n0,Tj) {  
  rho=param[1]  ###param order!!!!!
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21] 
  logl=Logl(dat2=cra2,dat1=dn1,rho,gamma=gamma,
                       beta=beta,eta=eta,n0=n0,rhoj=rhoj,Tj=Tj)
  names(logl)=c("Loglikelihood")
  return(logl)
} ##data dat2=cra2,dat1=dn1
Loglike(init.para,n0,Tj)

posterior <- function(param,n0,Tj){
  p=Loglike(param,n0,Tj) + Logprior(param)
  names(p)="Log posterior"
  return (p)
}
posterior(init.para,n0,Tj)

#-----------------------------------------------------
proposalfunction <- function(para){
  return(rnorm(1,mean = para, sd= 0.1))
}
## rho"gamma1"gamma2"gamma3"beta1"beta2"eta1"eta2"rhoj[1:13]
jacobian=function(orderpara,x){
  if (orderpara ==3||orderpara ==5){ return(exp(-x)/(1+exp(-x))^2)
  }else if (orderpara == 2) {return(min(cra2$Diameter)*exp(-x)/(1+exp(-x))^2)
  }else if(orderpara == 7) {return(5*exp(-x)/(1+exp(-x))^2)
  }else { return(exp(x))}
}
logit=function(x){log(x/(1-x))}
trans=function(orderpara,x){
  if (orderpara == 2) {return(logit(x/min(cra2$Diameter)))
  }else if(orderpara == 3||orderpara ==5){return(logit(x))
  }else if(orderpara ==7){return(logit(x/5))
  }else if(orderpara == 6||orderpara ==8){return(log(x-1))
  }else {return(log(x))}
}
backtrans=function(orderpara,x){
  if (orderpara == 2) {return(min(cra2$Diameter)/(1+exp(-x)))
  }else if(orderpara == 3||orderpara ==5){return(1/(1+exp(-x)))
  }else if(orderpara ==7){return(5/(1+exp(-x)))
  }else if(orderpara == 6||orderpara ==8){return(exp(x)+1)
  }else {return(exp(x))}
}

metropolis_MCMC <- function(startvalue, n0,Tj, iterations){
  chain = array(dim = c(iterations+1,21))
  chain[1,] = startvalue
  for (i in 1:iterations){
  for (j in 1:21){
    transpar=trans(j,chain[i,j])
    old_param= c( chain[i+1,0:(j-1)],chain[i,j], chain[i,-(0:j)] )
    repeat{
      proposal = proposalfunction(transpar)
      pro_param= c( chain[i+1,0:(j-1)],backtrans(j,proposal), chain[i,-(0:j)] )
      prob = exp(posterior(pro_param,n0,Tj) - posterior(old_param,n0,Tj)
                 + jacobian(j,proposal) - jacobian(j,transpar)) ##jump ratio
      if(!is.na(prob)){ break } else{
        cat("Error: Jumping ratio for para#",j,"is NA, redo proposal\n")}
    }
    
    if (runif(1) < prob) {
      chain[i+1,j] = backtrans(j,proposal)
    }else{
      chain[i+1,j] = chain[i,j]}
    
  }
  }
  return(chain)
}

xxxx=metropolis_MCMC(startvalue=init.para, n0=20,Tj,1)

#------update Tj with bern(p_tj) conditonal on all other params.
newT=function(param,dat1=dn1){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21]
  
  p_tj=NULL
  s1=dat1[,1]
  obs=dat1[,2]
  j=levels(obs)
  for (i in 1:13){
    if (sum(obs==j[i])!=0) {
      sj=s1[obs==j[i]] ## select data for observer i
      p_tj[obs==j[i]]= rho*qI_s(1,sj,gamma)*Pdf(sj,beta)/( rho*qI_s(1,sj,gamma)*Pdf(sj,beta) + rhoj[i]*Ppha(sj,eta) )
    }
  }
  p_tj
  ## update Tj to Tj.1
  Tj.1=NULL
  for (i in 1:length(Tj)){
    Tj.1[i]=rbinom(1,1,p_tj[i])
  }
  return(Tj.1)
}
newT(param=xxxx[2,])
#---------------------------------------
##-----update n0 with Poisson(rho*qI(0,gamma,beta)) conditonal on all other params.
newN0=function(param){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:21]
  
  n0.1=rpois(1,rho*qI(0,gamma,beta))
  return(n0.1)
}
newN0(param=xxxx[2,])
#-----------------------------------------------------
MH=function(startvalue, n0,Tj, iter){
  chain = array(dim = c(iter+1,22))
  Tjs=array(dim = c(iter+1,length(Tj)))
  chain[1,] = c(startvalue,n0=n0)
  Tjs[1,] = Tj
for (i in 1:iter){
  cat("inter=",i,"\n",chain[i,],"\n")
  repeat{
    xx= metropolis_MCMC(chain[i,1:21], chain[i,22],Tjs[i,],1)[2,]
   if(!any(is.na(xx))){ break } else {cat("redo metropolis mc for param")}
   }
  chain[i+1,1:21]=xx
  chain[i+1,22]= newN0(param=chain[i+1,1:21])
  Tjs[i+1,]=newT(param=chain[i+1,1:21])
  }
  return(list(chain=chain,Tjs=Tjs))
}
xxx=MH(init.para,20,Tj,70)


burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))










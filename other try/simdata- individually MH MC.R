load("/Users/zhang/Desktop/simdata.Rdata")

####### try simdata and plot!
tr <- rbind(simdata[1:68,1:3]*0.995,simdata[69:147,1:3],simdata[148:219,1:3]*1.002)
try <- cbind(tr,ob=simdata$V4)
plot(x,y,type = "n",xlim=c(0,10),ylim=c(0,10))
nseg=360;i=1
for (i in 1:nrow(try)){
  x.cent <- try[i,1]
  y.cent <- try[i,2]
  r <- try[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='red')
}
###clustering...it works!
hc <- hclust(dist(try))
plot(hc)
hcd <- as.dendrogram(hc)
plot(hcd,ylim = c(0,4))
abline(h=0.12,col="red")  ###plot the clustering and criteria line!


cuth <- cut(hcd, h = 0.12)   #### problem: how to set a "good" height criteria
nbranch <- length(cuth$lower)
plot(cuth$lower[[89]])      ### branches which height under 0.2
cuth$lower

try1<-as.matrix(try[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ### yeah!!! find index of rows!!!
  if (length(r)>1)  try1[r,] <- rep(colMeans(try1[r,]),each = length(r))
}

try_d <- cbind(try1,ob=try$ob)  ### pretreatment simdata


####### test ....
try_n2 <- sum(table(try_d[,3])!=1) ## repeated record N_2+
try_sj <- try_d[try_d[,1] %in% names(table(try_d[,1])[table(try_d[,1])==1]),]
## record only once
try_n2
try_yd <- table(try_sj[,"ob"])
try_yd
sum(try_n2,try_yd)

try.t =data.frame(try_d,1)
# obtain: size nrecord
(tryx=aggregate(try.t[,-4], by=list(rowSums(try.t[,-4])),FUN=sum)[,4:5]) 

tryn2 <- tryx[tryx$X1!=1,]

tryN2=nrow(tryn2); tryJ=3; tryn=as.integer(tryn2[,2]); trys=tryn2[,1]

trys1=tryx[tryx[,2]==1,][,1]
trydn1=try_d[,3:4][try_d[,3] %in% as.character(trys1),]
table(trydn1[,2])

#---------------------------------------------------------------
#---------------------------------------------------------------
#MCMC
#---------------------------------------------------------------
#----------------------------------------
J=3

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
  l_02= -rho*q2(gamma,beta) + sum( log(qI_s(n,s,gamma)*Pdf(s,beta)) ) -
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
}


##--------------------initial value for Tj
Tj=NULL
Tj[trydn1[,1] >2]=rbinom(sum(trydn1[,1] >2),1,0.5) 
Tj[trydn1[,1]<=2]=0


#------------------------------------------------------------
#---------------------------------------
##update other parameters using M-H algorithm. h = L(param)*g(param) 
## with dat2=tryn2; dat1=trydn1
## param order: rho"gamma1"gamma2"gamma3"beta1"beta2"eta1"eta2"rhoj[1:13]
#------------------------------------------
Logprior=function(param){
  rho=param[1]  ###param order!!!!!
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:length(param)]
  
  rho.p=dgamma(rho,6,6/100, log = T) ## intensity of craters
  gamma1.p=dunif(gamma[1],0, 2, log = T) ## criteria size for visible crater
  gamma2.p=dunif(gamma[2],0,1, log = T) ##limit prob for missing
  gamma3.p=dgamma(gamma[3],3, log = T)  ##decreasing rate for miss prob
  beta1.p=dunif(beta[1],0,1, log = T)  ## >0 according data <1
  beta2.p=dgamma(beta[2]-1,20,3/100, log = T) ## beta[2]>1, 0 < size < beta[2]  max(cra$Diameter)
  eta1.p=dunif(eta[1],0,5, log = T)
  eta2.p=dgamma(eta[2]-1,6,2/100, log = T) ## max(dn1$Diameter)
  
  rhoj.p=dgamma(rhoj,1,0.1, log = T) ## for vector rhoj[1:13]
  
  prior=rho.p+gamma1.p+gamma2.p+gamma3.p+beta1.p+beta2.p+eta1.p+eta2.p+sum(rhoj.p)
  names(prior)=c("logprior")
  return(prior)
}

init.para=c(rho=100,gamma=c(2,0.5,5),beta=c(.8,3),
            eta=c(1.8,2),rhoj=rep(3,3))
Logprior(init.para)

Loglike=function(param,n0,Tj) {  
  rho=param[1]  ###param order!!!!!
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:length(param)] 
  logl=Logl(dat2=tryn2,dat1=trydn1,rho,gamma=gamma,
            beta=beta,eta=eta,n0=n0,rhoj=rhoj,Tj=Tj)
  names(logl)=c("Loglikelihood")
  return(logl)
} ##data dat2=cra2,dat1=dn1
Loglike(init.para,n0=10,Tj)

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
  }else if (orderpara == 2) {return(2*exp(-x)/(1+exp(-x))^2)
  }else if(orderpara == 7) {return(5*exp(-x)/(1+exp(-x))^2)
  }else { return(exp(x))}
}
logit=function(x){log(x/(1-x))}
trans=function(orderpara,x){
  if (orderpara == 2) {return(logit(x/2))
  }else if(orderpara == 3||orderpara ==5){return(logit(x))
  }else if(orderpara ==7){return(logit(x/5))
  }else if(orderpara == 6||orderpara ==8){return(log(x-1))
  }else {return(log(x))}
}
backtrans=function(orderpara,x){
  if (orderpara == 2) {return(2/(1+exp(-x)))
  }else if(orderpara == 3||orderpara ==5){return(1/(1+exp(-x)))
  }else if(orderpara ==7){return(5/(1+exp(-x)))
  }else if(orderpara == 6||orderpara ==8){return(exp(x)+1)
  }else {return(exp(x))}
}

metropolis_MCMC <- function(startvalue, n0,Tj, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    for (j in 1:length(startvalue)){
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

x=metropolis_MCMC(startvalue=init.para, n0=10,Tj,1)

#------update Tj with bern(p_tj) conditonal on all other params.
newT=function(param,dat1=trydn1){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:length(param)]
  
  p_tj=NULL
  s1=dat1[,1]
  obs=dat1[,2]
  j=levels(as.factor(obs))
  for (i in 1:J){
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
newT(param=x[2,])
#---------------------------------------
##-----update n0 with Poisson(rho*qI(0,gamma,beta)) conditonal on all other params.
newN0=function(param){
  rho=param[1]  ###param order
  gamma=param[2:4]
  beta=param[5:6]
  eta=param[7:8]
  rhoj=param[9:length(param)]
  
  n0.1=rpois(1,rho*qI(0,gamma,beta))
  return(n0.1)
}
newN0(param=xxxx[2,])
#-----------------------------------------------------
MH=function(startvalue, n0,Tj, iter){
  n=length(startvalue)
  chain = array(dim = c(iter+1,n+1))
  Tjs=array(dim = c(iter+1,length(Tj)))
  chain[1,] = c(startvalue,n0=n0)
  Tjs[1,] = Tj
  for (i in 1:iter){
    cat("inter=",i,"\n",chain[i,],"\n")
    repeat{
      xx= metropolis_MCMC(chain[i,1:n], chain[i,n+1],Tjs[i,],1)[2,]
      if(!any(is.na(xx))){ break } else {cat("redo metropolis mc for param")}
    }
    chain[i+1,1:n]=xx
    chain[i+1,n+1]= newN0(param=chain[i+1,1:n])
    Tjs[i+1,]=newT(param=chain[i+1,1:n])
  }
  return(list(chain=chain,Tjs=Tjs))
}
xxx=MH(init.para,n0=10,Tj,30)


burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))








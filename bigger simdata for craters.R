###simulate a crater data

####set:real density rho, area A , crater number N~Poisson(rho*A)
set.seed(43); N <- rpois(1,2000) ###rho=2000
N ##crater number
set.seed(121); x <- runif(N,0,2000)
set.seed(17); y <- runif(N,0,2000) 
# x;y##location
plot(x,y,xlim=c(0,2000),ylim=c(0,2000))

##size pdf : f(s|theta0,theta1)=theta1 * exp(theta0) * s^(theta1-1) where s=diameter
### 0< s <exp(-theta0/theta1) where set 

##upper=exp(-theta0/theta1);integrate(f,lower=0,upper=upper)
#theta1=0.4; theta0=-0.2
#f <- function(s) theta1 * exp(theta0) * s^(theta1-1)
#hh <- seq(0,exp(-theta0/theta1),length.out=1000)
#plot(hh,f(hh),type="l",ylim=c(0,20),xlim=c(0,exp(-theta0/theta1)))
#set.seed(233);u <- runif(N,0,1)
## size <- (exp(-theta0)* u)^(1/theta1)   ###sampling size from s=F^-1(uniform)

#---gamma dsn for size
curve(dgamma(x,2.3,.1),0,100,col=2)  ### <--size param
set.seed(233);size=rgamma(N,2.3,.1)
round(size,4)
hist(size,breaks=50,probability = T)
head(sort(size),30)
real<- data.frame(x,y,size)

plot(x,y,type = "n",xlim=c(0,2000),ylim=c(0,2000))
nseg=360;i=1
for (i in 1:N){
  x.cent <- real[i,1]
  y.cent <- real[i,2]
  r <- real[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='red')
}


## miss prob p(s|t,r)=I(s<r)+t*I(s>=r)  set r=0.01, t=0.2 ==> Bernoulli
## phantom pdf may simply be h(s) = Unif(0,1)...

## phantom intensity : set rho1*A=50 F~Poisson(rho1*A)  

see<-real[-which(size<3),]  ##<-- observing criteria=3
nsmall <-sum(size<3)  ##miss for each

####now consider 10 observers:

set.seed(63);miss <- rbinom(10*(N-nsmall),1,0.2) ## miss 
sum(miss)/sum(10*(N-nsmall)) ### 1=miss
nsee=N-nsmall
n=NULL
for (i in 1:10){
  n[i]=sum(miss[1:nsee+nsee*(i-1)]==1)
}
n  #  miss # for observers 


set.seed(1);phantom <- rpois(10,50)  ## <--obs=10, phantom rho1*A=50
phantom ### phantom num

set.seed(9);ph.size <- rgamma(sum(phantom),3,.1)###gamma dsn for phantom size
hist(ph.size,breaks=100,probability = T)
set.seed(89);ph.x <- runif(sum(phantom),0,2000);ph.y <- runif(sum(phantom),0,2000)
plot(ph.x,ph.y,xlim=c(0,2000),ylim=c(0,2000),"p")

ph <-as.data.frame(cbind(ph.x,ph.y,ph.size,rep(1:10,phantom)))
colnames(ph) = c(colnames(real),"ob")
phdat = split(ph,ph[,4])  ## phantom data for 10 observers


### obs data
simdata=NULL
obdat=NULL
for (i in 1:10){
  x=cbind(see[which(miss[1:nsee+nsee*(i-1)]==0),],ob=i)
  obdat[[i]]=rbind(x,phdat[[i]])
  row.names(obdat[[i]]) = 1:nrow(obdat[[i]])
  simdata=rbind(simdata,obdat[[i]])
}
obdat[[7]]  ###  list for 10 observers (real-miss+phantom)

simdata  ### simdata :data frame of all observations

save(simdata,file="/Users/zhang/Desktop/simdata.Rdata",ascii = F)


nob=table(simdata$ob)

####### try simdata and plot!
set.seed(689);obsize <- simdata$size*rnorm(nrow(simdata),1,0.005)  ###add random error for size
try <- cbind(simdata[,1:2],size=obsize,ob=simdata[,4])
x=NULL;y=NULL
plot(x,y,type = "n",xlim=c(0,2000),ylim=c(0,2000))
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
dtry <- cbind(try[,1:2],try[,3]*.5)   ### dist=sqrt(dx^2+dy^2+dr^2*coef) 
fn = function(x,y) sqrt(sum((x-y)^2))/min(x[3],y[3])  ### dist/min(r1,r2) define!!
datdist <- proxy::dist(dtry, method=fn)
hc <- hclust(datdist,method="average") 
### cluster method:one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
hcd <- as.dendrogram(hc)
plot(hcd,ylim = c(0,10),xlim=c(1000,1500))
abline(h=.2,col="red")  ###plot the clustering and criteria line!


cuth <- cut(hcd, h = 0.18)   #### problem: how to set a "good" height criteria
nbranch <- length(cuth$lower)
plot(cuth$lower[[109]])      ### branches which height under 0.2
cuth$lower

try1<-as.matrix(try[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ### yeah!!! find index of rows!!!
  if (length(r)>1)  try1[r,] <- rep(colMeans(try1[r,]),each = length(r))
}

try_d <- cbind(try1,ob=try$ob)  ### pretreatment simdata!!! got it

##plot after clustering
x<-NULL;y<-NULL
plot(x,y,type = "n",xlim=c(0,2000),ylim=c(0,2000))
nseg=360;i=1
for (i in 1:nrow(try)){
  x.cent <- try[i,1]
  y.cent <- try[i,2]
  r <- try[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='green')
}

nseg=360;i=1
for (i in 1:nrow(try1)){
  x.cent <- try1[i,1]+.1
  y.cent <- try1[i,2]
  r <- try1[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='red')
}

####### test ....
try_n2 <- sum(table(try_d[,3])!=1) ## repeated record N_2+
try_sj <- try_d[try_d[,1] %in% names(table(try_d[,1])[table(try_d[,1])==1]),]
## record only once
try_n2
try_yd <- table(try_sj[,"ob"])
try_yd  ## almost the same as phantom
sum(try_n2,try_yd)  ## individual circles on the plot

try.t =data.frame(try_d,1)
# obtain: size nrecord
tryx= cbind( as.numeric(names(table(try.t[,3]))),table(try.t[,3]) ) 
row.names(tryx)=1:nrow(tryx)
tryx

tryn2 <- tryx[tryx[,2]!=1,]

#tryN2=nrow(tryn2); tryJ=3; tryn=as.integer(tryn2[,2]); trys=tryn2[,1]

trys1=tryx[tryx[,2]==1,][,1]
tryn1=try_d[,3:4][try_d[,3] %in% as.character(trys1),] ## same as try_sj

save(tryn2,file="/Users/zhang/Desktop/bsimdata-2.Rdata")
save(tryn1,file="/Users/zhang/Desktop/bsimdata-1.Rdata")

hist(size,probability = T,breaks=150)
hist(tryx[,1],probability = T,breaks=150)
curve(dgamma(x,2.3,.1),0,120,add=T,col=2)


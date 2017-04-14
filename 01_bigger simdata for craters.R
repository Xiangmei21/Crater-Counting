###simulate a crater data

plotcircle = function(dat,colour="green"){
  x<-NULL;y<-NULL
  plot(x,y,type = "n",xlim=c(min(dat[,1])-max(dat[,3]),max(dat[,1])+max(dat[,3])),
       ylim=c(min(dat[,2])-max(dat[,3]),max(dat[,2])+max(dat[,3])))
  nseg=360;i=1
  for (i in 1:nrow(dat)){
    x.cent <- dat[i,1]
    y.cent <- dat[i,2]
    r <- dat[i,3]/2
    xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
    yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
    lines(xx,yy, col=colour)
  }
}


linecircle = function(data,colour="red"){
  nseg=360;i=1
  for (i in 1:nrow(data)){
    x.cent <- data[i,1]+1
    y.cent <- data[i,2]
    r <- data[i,3]/2
    xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
    yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
    lines(xx,yy, col=colour)
  }
}



###set:real density rho, area A , crater number N~Poisson(rho*A)
set.seed(43); N <- rpois(1,2000) ###rho=2000
N ##crater number
## simulate locations
set.seed(82); locations <- matrix(runif(4*N,0,2000),byrow = T,ncol=2)
## remove overlaps
simtree <- as.dendrogram(hclust(dist(locations),method="average"))
plot(simtree,ylim = c(0,100),xlim=c(1000,1500),leaflab = "none")
cutsim <- cut(simtree, h = 19)   ## at least N non overlaping locations
nsim <- length(cutsim$lower)  
new_loc=NULL
for (i in 1:nsim){  
  leaf=cutsim$lower[[i]]
  r <- as.numeric(unlist(leaf))   
  new_loc <- rbind(new_loc,locations[r[1],]) 
}
## randomlly sample N locations out of nsim non overlapping locations
set.seed(59);real_loc <- new_loc[sample(1:nsim,N),] 
x = real_loc[,1]; y = real_loc[,2]
plot(x,y,xlim=c(0,2000),ylim=c(0,2000))
#---gamma dsn for size
curve(dgamma(x,2.3,.1),0,100,col=2)  ### <--size param
set.seed(233);size=rgamma(N,2.3,.1)
round(size,4)
hist(size,breaks=50,probability = T)
head(sort(size),30)
real<- data.frame(x,y,size)

plotcircle(real)

## miss prob p(s|t,r)=I(s<r)+t*I(s>=r)  ==> Bernoulli
## phantom pdf may simply be h(s) = Unif(0,1)...

## phantom intensity : set rho1*A=50 F~Poisson(rho1*A)  

see<-real[-which(size<3),]  ##<-- observing criteria=3
nsmall <-sum(size<3)  ## miss # for small size 

####now consider 10 observers:

set.seed(63);miss <- rbinom(10*(N-nsmall),1,0.2) ## miss 
sum(miss)/sum(10*(N-nsmall)) ### 1=miss
nsee=N-nsmall
miss_byeach=matrix(nrow=10,ncol=nsee)
miss_byall=rep(1,nsee)
for (i in 1:10){
  miss_byeach[i,]=miss[1:nsee+nsee*(i-1)]
  miss_byall = miss_byall & miss[1:nsee+nsee*(i-1)]
}
sum(miss_byall)  #  miss # for all observers 
rowSums(miss_byeach)     #  miss # for each observer
which(colSums(miss_byeach) >7)  # misspattern of each crater 


set.seed(1);phantom <- rpois(10,50)  ## <--obs=10, phantom rho1*A=50
phantom ### phantom num

set.seed(9);ph.size <- rgamma(sum(phantom),3,.1)###gamma dsn for phantom size
hist(ph.size,breaks=100,probability = T)

set.seed(59);ph_loc <- new_loc[-sample(1:nsim,N),][1:sum(phantom),]
set.seed(90);ph.x <- ph_loc[,1]; ph.y <- ph_loc[,2]
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

###add random error for sizes and locations
set.seed(689);simdata$size <- simdata$size*rnorm(nrow(simdata),1,0.001)
set.seed(687);simdata$x <- simdata$x*rnorm(nrow(simdata),1,0.0001)
set.seed(686);simdata$y <- simdata$y*rnorm(nrow(simdata),1,0.0001)

save(simdata,file="./data/simdata.Rdata",ascii = F)
# load("./data/simdata.Rdata")

nob=table(simdata$ob)

####### try simdata and plot!
try <- simdata
plotcircle(try)

###clustering...
dtry <- cbind(try[,1:2],try[,3]*0.5)   ### dist=sqrt(dx^2+dy^2+dr^2*coef) 
fn = function(x,y) sqrt(sum((x-y)^2))/(min(x[3],y[3])) ### dist/min(r1,r2) define!!
datdist <- proxy::dist(dtry, method=fn)
hc <- hclust(datdist,method="single") 
### cluster method:one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
hcd <- as.dendrogram(hc)
# save(hcd, file="./data/hcd_simdata_clustering.Rdata")
load("./data/hcd_simdata_clustering.Rdata")
load("./data/hcd_single_simdata_clustering.Rdata")
load("./data/hcd_complete_simdata_clustering.Rdata")
load("./data/hcd_ward_simdata_clustering.Rdata")

plot(hcd_single,ylim = c(0,10),xlim=c(1800,2500),leaflab = "none")
abline(h=.5,col="red")  ###plot the clustering and criteria line!

cuth <- cut(hcd_single, h = .5)   #### set a "good" height criteria
nbranch <- length(cuth$lower)
leaves=NULL
for (i in 1:nbranch){
  leaves[i]=length(unlist(cuth$lower[[i]]))
}
# ------ check inconsistency---------
max(leaves) ##need <=10
which(leaves >10 )
plot(cuth$lower[[which.max(leaves)]]) 
for (i in which(leaves >10 )){
  plot(cuth$lower[[i]],main=i) 
  plotcircle(try[unlist(cuth$lower[[i]]),])
}

#-------------------------------------

try1<-as.matrix(try[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ### yeah!!! find index of rows!!!
  if (length(r)>1)  try1[r,] <- rep(colMeans(try1[r,]),each = length(r))
}

try_d <- cbind(try1,ob=try$ob)  ### pretreatment simdata!!! got it

##plot after clustering
plotcircle(try)

linecircle(try1)

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

hist(tryn1[,1],breaks=50,probability = T)
hist(ph.size,breaks=50,probability = T)

save(tryn2,file="./data/bsimdata-2.Rdata")
save(tryn1,file="./data/bsimdata-1.Rdata")

hist(size,probability = T,breaks=150)
hist(tryx[,1],probability = T,breaks=150)
curve(dgamma(x,2.3,.1),0,120,add=T,col=2)


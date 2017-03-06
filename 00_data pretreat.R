da <- read.table("./Crater_Meas_data.txt", header = T)
s <- da[,3]
hist(s,breaks = 300,xlim = c(0,100),probability = T) ##size hist: shape like gamma dist

table(da$Observer,da$Image)
########
library(tidyverse)
library(ggforce)
da %>% filter(Image=="NAC",Observer=="Antonenko1") %>%
  ggplot(aes(x0=X,y0=Y,r=Diameter/2)) +
  geom_circle() +
  geom_circle(aes(x0=X,y0=Y,r=Diameter/2),color="red",data=da %>% filter(Image=="NAC",Observer=="Antonenko3"))
  
#########################################################################
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

##  calculate distance between rows --- not good.
#dist <- as.matrix(dist(try))
#dist[upper.tri(dist,diag = T)] <- NA
#similar <- which(dist<0.3,arr.ind=T)  ### similar pairs..do little help....
### plus how to set a "good" criteria ?

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

try_d <- cbind(try1,ob=try$ob)


####### test ....
try_n2 <- sum(table(try_d[,3])!=1) ## repeated record N_2+
try_sj <- try_d[try_d[,1] %in% names(table(try_d[,1])[table(try_d[,1])==1]),]
## record only once
try_n2
try_yd <- table(try_sj[,"ob"])
try_yd
sum(try_n2,try_yd)


############  ##--------pretreatment--------##  ################

#### choose a rectangle of crater image
head(da)
dat <- da[,1:4][ da$X>1000 & da$X<2000 & da$Y> -1500 & da$Y< -1000, ]
table(dat$Observer)
x<-NULL;y<-NULL
plot(x,y,type = "n",xlim=c(1000,2000),ylim=c(-1500,-1000))
nseg=360;i=1
for (i in 1:nrow(dat)){
  x.cent <- dat[i,1]
  y.cent <- dat[i,2]
  r <- dat[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='green')
}

##########################################start  clustering...

d_dat <- cbind(dat[,1:2],dat[,3]*0.5)   ### dist=sqrt(dx^2+dy^2+dr^2*0.25) good enough???
fn = function(x,y) sqrt(sum((x-y)^2))/min(x[3],y[3])  ### dist/min(r1,r2) define!!
datdist <- proxy::dist(d_dat, method=fn)
hc <- hclust(datdist, "average") 
### cluster method=single "min-linkage" 
plot(hc)
hcd <- as.dendrogram(hc)
plot(hcd)
abline(h=1.2,col="red")  ###plot the clustering and criteria line!


cuth <- cut(hcd, h = 1.2)   #### problem: how to set a "good" height criteria ???
nbranch <- length(cuth$lower)
plot(cuth$lower[[28]])      ### branches which height under 0.2
cuth$lower      #### str of leaf  ~~

length(unlist(cuth$lower[[28]]))

nrecord <- NULL;srecord <- NULL

dat1<-as.matrix(dat[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ### yeah!!! find index of rows!!!
  if (length(r)>1)  {dat1[r,] <- rep(colMeans(dat1[r,]),each = length(r))}
  nrecord[i] <- length(r)
  srecord[i] <- mean(dat1[r,3])
}

nrecord
srecord
plot(srecord,nrecord)
misspattern <- data.frame(srecord,nrecord)

dat_d <- data.frame(dat1,ob=dat[,4])

####### data
dat_n2 <- sum(table(rowSums(dat_d[,1:2]))!=1)  ## repeated record N_2+
jud <- rowSums(dat_d[,1:2])
dat_sj <- dat_d[jud %in% names(table(jud)[table(jud)==1]),]
## record only once
dat_n2  ### = sum(nrecord!=1)
dat_yd <- table(dat_sj[,"ob"])
dat_yd  ##
sum(dat_n2,dat_yd)
nbranch

length(table(rowMeans(dat1))) #####

##plot after clustering
x<-NULL;y<-NULL
plot(x,y,type = "n",xlim=c(1000,2000),ylim=c(-1500,-1000))
nseg=360;i=1
for (i in 1:nrow(dat)){
  x.cent <- dat[i,1]
  y.cent <- dat[i,2]
  r <- dat[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='green')
}

nseg=360;i=1
for (i in 1:nrow(dat_d)){
  x.cent <- dat_d[i,1]+1
  y.cent <- dat_d[i,2]
  r <- dat_d[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='red')
}

dat_yd  ################# Sj seems not to be a homogeneous Poisson for all observers
## rho may not be a constant in the crater image

####-----------------------------------------------------------------miss pattern
misspattern <- data.frame(srecord,nrecord)
rep <- misspattern[nrecord!=1,]
cut <- cut(rep$srecord,breaks = c(1:470)/2)
x=aggregate(rep$nrecord~cut,FUN=mean)
barplot(x[,2],xlim=c(0,90)) ## height= average # of records for size=x
x.cut=sub("]","",gsub(".*,","", as.character(x[,1])))
plot(x.cut,x[,2],xlim=c(0,50),"l")

#-----------------------------------------------------------------
dat <- da[,1:4][ da$X>1000 & da$X<2000 & da$Y> -1000 & da$Y< -500, ]
d_dat <- cbind(dat[,1:2],dat[,3]*0.5)   ### dist=sqrt(dx^2+dy^2+dr^2*0.25) 
fn = function(x,y) sqrt(sum((x-y)^2))/min(x[3],y[3])  ### dist/min(r1,r2) define!!
datdist <- proxy::dist(d_dat, method=fn)
hc <- hclust(datdist, "average") 
### cluster method=single "min-linkage" 
hcd <- as.dendrogram(hc)
nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.17, col = "blue")
plot(hcd,ylim=c(-1,20),xlim=c(0,300),ylab="Height",main="Cluster Dendrogram",
     nodePar = nodePar,leaflab = "none")
abline(h=1,col="red")  ###plot the clustering and criteria line!
text(x=-5,y=1.5,"height=1",col=2,cex=0.7)

cuth <- cut(hcd, h = 1)   #### problem: how to set a "good" height criteria 
nbranch <- length(cuth$lower)

nrecord <- NULL;srecord <- NULL
datt<-as.matrix(dat[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ### yeah!!! find index of rows!!!
  if (length(r)>1)  {datt[r,] <- rep(colMeans(datt[r,]),each = length(r))}
  nrecord[i] <- length(r)
  srecord[i] <- mean(datt[r,3])
}
##plot after clustering
x<-NULL;y<-NULL
plot(x,y,type = "n",xlim=c(1000,2000),ylim=c(-1000,-500),xlab="X",ylab="Y",main="Pretreatment Process")
nseg=360;i=1
for (i in 1:nrow(dat)){
  x.cent <- dat[i,1]
  y.cent <- dat[i,2]
  r <- dat[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='green')
}

nseg=360;i=1
for (i in 1:nrow(datt)){
  x.cent <- datt[i,1]+1
  y.cent <- datt[i,2]
  r <- datt[i,3]/2
  xx <- x.cent + r*cos( seq(0,2*pi, length.out=nseg) )
  yy <- y.cent + r*sin( seq(0,2*pi, length.out=nseg) )
  lines(xx,yy, col='red',lwd=0.5)
}


dat.t =data.frame(datt,1)
# obtain: size nrecord
(xxx=aggregate(dat.t, by=list(rowSums(dat.t)),FUN=sum)[,4:5])  


observer <- dat[,4]
pretreat <- data.frame(datt,observer)
save(pretreat,file="/Users/zhang/Desktop/crater rectangle.Rdata",ascii = F)
s1=xxx[xxx[,2]==1,][,1]

sort(s1)==sort(pretreat[,3][pretreat[,3] %in% as.character(s1)])

### Nj data only seen by 1 person
dn1=pretreat[,3:4][pretreat[,3] %in% as.character(s1),]
sort(s1)==sort(dn1[,1])

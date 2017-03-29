da <- read.table("./Crater_Meas_data.txt", header = T)
table(da$Observer,da$Image)
########
library(tidyverse)
library(ggforce)
wac=da %>% filter(Image=="WAC")
max(wac$X)
min(wac$Y)
nac=da %>% filter(Image=="NAC") ## nac = da[da$Image=="NAC",]
max(nac$X)
min(nac$Y)

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
plotcircle(nac[15100:15380,],"red")

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
linecircle(nac[15200:15380,],"blue")

#######------------start  clustering----------------------------
# deal with nac

# set tuning parameter for distance = 0.25
d.nac <- cbind(nac[,1:2],nac[,3]*0.5)   
### dist=sqrt(dx^2+dy^2+dr^2*0.25) 
fn = function(x,y) sqrt(sum((x-y)^2))/min(x[3],y[3])  ### dist/min(r1,r2) define!!
datdist <- proxy::dist(d.nac, method=fn)
hc <- hclust(datdist, "average") 
### cluster method="average" 

hcd <- as.dendrogram(hc)
save(hcd,file="nac_cluster.RData")
plot(hcd,xlim=c(1600,1700),ylim=c(-1,10))
abline(h=0.7,col="red")  ### plot the clustering and criteria line!

cuth <- cut(hcd, h = 0.7)   #### set a "good" height criteria 
nbranch <- length(cuth$lower)
cuth$lower      #### str of leaf 
leaves=NULL
for (i in 1:nbranch){
  leaves[i]=length(unlist(cuth$lower[[i]]))
}
# ------ check inconsistency---------
max(leaves) ##need <=12
plot(cuth$lower[[which.max(leaves)]])      ### eg: branch which height under cut line
maxindex=as.numeric(unlist(cuth$lower[[which.max(leaves)]]))
plotcircle(nac[maxindex,])
linecircle(nac[maxindex[-1],]) ## check inconsistency when number of leaves >12
nac[maxindex,]
plotcircle(nac[c(2390,2347),]) ## circles in maxindex by the same observer

### clustering data---------------------------------
nrecord <- NULL;srecord <- NULL
nac.clustered<-as.matrix(nac[,1:3])
for (i in 1:nbranch){
  leaf=cuth$lower[[i]]
  r <- as.numeric(unlist(leaf))      ###  find index of rows!!!
  if (length(r)>1)  {nac.clustered[r,] <- rep(colMeans(nac[r,1:3]),each = length(r))}
  nrecord[i] <- length(r)
  srecord[i] <- mean(nac[r,3])
}

nrecord  ## number of repeat records for clustered craters
srecord  ## average size for craters
plot(srecord,nrecord)
obpattern <- data.frame(srecord,nrecord)

nac_pretreat <- data.frame(nac.clustered,ob=nac[,4])
save(nac_pretreat,file="nac_pretreated.RData")
#load("./nac_pretreated.RData")

####### data extract
length(unique(nac_pretreat[,1]-nac_pretreat[,2]+nac_pretreat[,3])) ## different craters should =nbranch
uni_crater=table(nac_pretreat[,1]-nac_pretreat[,2]+nac_pretreat[,3])
nac_2 <- sum(uni_crater!=1)  ## repeated record N_2+  should = sum(nrecord!=1)
jud <- nac_pretreat[,1]-nac_pretreat[,2]+nac_pretreat[,3] ## judgment criteria of different craters
nac_sj <- nac_pretreat[jud %in% names(table(jud)[table(jud)==1]),]
## record only once
nac_2  ### = sum(nrecord!=1)
(nac_1 <- table(nac_sj[,"ob"]))
sum(nac_1)  ### = sum(nrecord==1)
#check inconsistency
nbranch == sum(nac_2,nac_1)

##plot after clustering
plotcircle(nac)

linecircle(nac_pretreat)

nac_yd  ################# Sj seems not to be a homogeneous Poisson for all observers
## rho may not be a constant in the crater image

####---------------------I:observation pattern-----------------
obpattern <- data.frame(srecord,nrecord)
rep <- nac_n2## obpattern[nrecord!=1,]
cut <- cut(rep$srecord,breaks = c(10:800)/2)
x=aggregate(rep$nrecord~cut,FUN=mean)
barplot(x[,2],xlim=c(0,100)) ## height= average # of records for size=x
x.cut=sub("]","",gsub(".*,","", as.character(x[,1])))
plot(as.numeric(x.cut),x[,2],xlim=c(0,100),"l",xlab ="crater size",ylab="number of observers who saw the cater")
lines(as.numeric(x.cut),lowess(x[,2],f=.2)$y,xlim=c(0,100),"l",add=T,col=2)
min(srecord)
hist(nac$Diameter,breaks = 500)

smoothpattern=data.frame(x=as.numeric(x.cut),y=x[,2])
ggplot(data=smoothpattern,aes(x=x,y=y))+geom_line()+geom_smooth(span=0.2)+
  theme_classic()+xlab("crater size")+ylab("number of observers who saw the cater")+
  ggtitle(expression(paste("Observation pattern")))+
  theme(plot.title = element_text(hjust = 0.5,size=16, face="bold"),
        axis.title.y=element_text(vjust = 0.5),
        axis.text=element_text(size=12),
        axis.title = element_text(color="black", size=14))

## data save
nac_n2=obpattern[nrecord!=1,]
nac_n1=nac_sj[,3:4]
save(nac_n1,file="nac_n1.RData")
save(nac_n2,file="nac_n2.RData")
write.table(nac_n1,"nac_n1.txt",row.names = F)
write.table(nac_n2,"nac_n2.txt",row.names = F)

##### plot nac_n1
table(nac_n1$ob)
plotcircle(nac_sj,"red")






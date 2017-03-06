simn2 <- sum(table(simdata[,3])!=1) ## repeated record N_2+
simsj <- simdata[simdata[,1] %in% names(table(simdata[,1])[table(simdata[,1])==1]),]
## record only once
simn2
simyd <- table(simsj[,4])
simyd  ## almost the same as phantom
sum(simn2,simyd)  ## individual circles on the plot

sim.t =data.frame(simdata,1)
# obtain: size nrecord
simx= cbind( as.numeric(names(table(sim.t[,3]))),table(sim.t[,3]) ) 
row.names(simx)=1:nrow(simx)
simx

simn2 <- simx[simx[,2]!=1,]

#simN2=nrow(simn2); simJ=3; simn=as.integer(simn2[,2]); sims=simn2[,1]

sims1=simx[simx[,2]==1,][,1]
simn1=simdata[,3:4][simdata[,3] %in% as.character(sims1),] ## same as simsj

save(simn2,file="/Users/zhang/Desktop/sim-2.Rdata")
save(simn1,file="/Users/zhang/Desktop/sim-1.Rdata")

hist(simx[,1],probability = T,breaks=150)

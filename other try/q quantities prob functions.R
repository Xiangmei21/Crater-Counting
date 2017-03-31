library(rstan)
simdata
d <- table(simdata$V4)
d ###num of obs
n2 <- sum(table(simdata$x)!=1) ## repeated record N_2+
sj <- simdata[simdata$x %in% names(table(simdata$x)[table(simdata$x)==1]),]
   ## record only once
n2
y <- table(sj$V4)
y

J=3 ## observers
theta <- c(0.4,-0.2)
gamma <- c(0.01,0.2)


## size pdf
fs <- function(s) theta[1] * exp(theta[2]) * s^(theta[1]-1) ### 0< s <exp(-theta2/theta1)
## miss prob for size=s
pmiss <- function(s) ifelse(s>=gamma[1], gamma[2],1)
## prob of a crater size=s is seen once
qjs <- function(s) pmiss(s)^(J-1)*(1-pmiss(s))
## prob of a crater is seen once 
qj <- integrate(function(s) fs(s)*qjs(s),lower=0,upper=exp(-theta[2]/theta[1]) )$value

## prob of a crater size=s is not seen
q0s <- function(s) pmiss(s)^J
## prob of a crater size=s is seen more than twice
q2s <- function(s) {1-q0s(s)-qjs(s)*J}

q2 <- integrate(function(s) fs(s)*q2s(s),lower=0,upper=exp(-theta0/theta1) )$value

qj;q2
(q0=1-qj*J-q2)




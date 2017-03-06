
load("\\\\my.files.iastate.edu\\Users\\xmzhang\\Desktop\\crater counting\\sim-1.Rdata")
load("\\\\my.files.iastate.edu\\Users\\xmzhang\\Desktop\\crater counting\\sim-2.Rdata")

(dimj=table(simn1[,2]))

simj1=simn1[simn1$V4=="ob1",1]
simj2=simn1[simn1$V4=="ob2",1]
simj3=simn1[simn1$V4=="ob3",1]
simj=simn1[,1]


#----------------------------------------
J=3


library(rstan)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)


crater.model = " 
data {
int<lower=0> J;
int<lower=0> N2;              //number of repeated craters
int<lower=0> N1;  
int<lower=0> k[J]; 
vector[N2] dat2s;  
int<lower=1> dat2n[N2];     //pairs of size and n(#repeats)
vector[k[1]] datj1;  //size for observerid =1 
vector[k[2]] datj2;
vector[k[3]] datj3;
}
parameters {
real<lower=0> rho;
vector<lower=0>[J] rhoj;
vector<lower=0>[2] gamma;  //para in pmiss
vector<lower=0>[2] beta;  // para in gamma size density
vector<lower=0>[2] eta;   // para in phantom size pdf
real<lower=0> N0;
}
transformed parameters { 
real<lower=0,upper=1> p;
real<lower=0> qj;
real<lower=0> q2;
real<lower=0> q0;

p <- gamma[2];
qj <- (1-p)*p^(J-1)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q2 <- (1-J*p^(J-1)+(J-1)*p^J)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q0 <- gamma_cdf(gamma[1],beta[1],beta[2])+p^J*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
}
model {

gamma[1] ~ exponential(10);           // ??? gamma1 :small crater r 
gamma[2] ~ uniform(0,1);            //prior ?? gamma2 :prob miss when s>r
beta[1] ~ gamma(3,.1);              //+
beta[2] ~ gamma(3,.1);             //
eta[1] ~ gamma(3,.1);              //+
eta[2] ~ gamma(3,.1); 

N2 ~ poisson(q2*rho);
for (i in 1:N2){
dat2s[i] ~ gamma(beta[1],beta[2]);   
dat2n[i] ~ binomial(J,(1-p)) T[2,];
}

k[1] ~ poisson(qj*rho+rhoj[1]);
k[2] ~ poisson(qj*rho+rhoj[2]);
k[3] ~ poisson(qj*rho+rhoj[3]);
for (i in 1:k[1]){
if(datj1[i]<gamma[1]){
increment_log_prob( bernoulli_log(0 ,qj*rho/(qj*rho+rhoj[1])) + gamma_log(datj1[i] ,eta[1],eta[2]));
} else { 
increment_log_prob(log_mix(qj*rho/(qj*rho+rhoj[1]),
                  gamma_log(datj1[i] ,beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                  gamma_log(datj1[i] ,eta[1],eta[2])));
}}

for (j in 1:k[2]){
if(datj2[j]<gamma[1]){
increment_log_prob( bernoulli_log(0 ,qj*rho/(qj*rho+rhoj[2])) + gamma_log(datj2[j] ,eta[1],eta[2]));
} else { 
increment_log_prob( log_mix(qj*rho/(qj*rho+rhoj[2]),
                  gamma_log(datj2[j] ,beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                  gamma_log(datj2[j] ,eta[1],eta[2])));
}}
for (n in 1:k[3]){
if(datj3[n]<gamma[1]){
increment_log_prob( bernoulli_log(0 ,qj*rho/(qj*rho+rhoj[3])) + gamma_log(datj3[n] ,eta[1],eta[2]));
} else { 
increment_log_prob( log_mix(qj*rho/(qj*rho+rhoj[3]),
                 gamma_log(datj3[n] ,beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                 gamma_log(datj3[n] ,eta[1],eta[2])));
}}

increment_log_prob( -q0*rho+N0*log(q0*rho)-lgamma(N0+1));
}
"  
 
sdat = list(J=3,N2=nrow(simn2),N1=nrow(simn1),dat2s=simn2[,1],dat2n=simn2[,2],datj= simj)
fit = stan(model_code = crater.model, data = sdat, iter = 1000, chains=4)
fit
##para values
#rho = 200 , rhoj=10
#theta <- c(0.4,-0.2)  #gamma <- c(0.01,0.2)
#miss 21 (<gamma1) + (36 34 32)randomly miss
plot(fit)
e <- extract(fit)
plot(x=1:length(e$gamma[,2]),e$gamma[,2],"l")
plot(x=e$gamma[,2],y=e$rho)

beta1=mean(e$beta[,1])
beta2=mean(e$beta[,2])
curve(dgamma(x,beta1,beta2),0,2,ylim=c(0,3))
theta1=0.4; theta0=-0.2
f <- function(s) theta1 * exp(theta0) * s^(theta1-1)
curve(f,0,2,col=2,add=T)
hist(simdata$size,freq = F,breaks=50,add=T)
fit


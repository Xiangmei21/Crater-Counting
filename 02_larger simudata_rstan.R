
load("\\\\my.files.iastate.edu\\Users\\xmzhang\\Desktop\\crater counting\\bsimdata-1.Rdata")
load("\\\\my.files.iastate.edu\\Users\\xmzhang\\Desktop\\crater counting\\bsimdata-2.Rdata")

(dj=table(tryn1[,2]))


tryj=tryn1[,1]


#----------------------------------------


library(rstan)

crater.model = " 
data {
int<lower=0> J;            // number of observers
int<lower=0> N2;              //number of repeated craters
int<lower=0> N1;  

vector[N2] dat2s;  
int<lower=1> dat2n[N2];     //pairs of size and n(#repeats)
vector[N1] datj;  //size  

}
parameters {
real<lower=0> rho;
real<lower=0> rhoj;
vector<lower=0>[2] gamma;  //para in pmiss
vector<lower=0>[2] beta;  // para in gamma size density
vector<lower=0>[2] eta;   // para in phantom size pdf
real<lower=0> N0;
}
transformed parameters { 
real<lower=0> qj;
real<lower=0> q2;
real<lower=0> q0;
real<lower=0,upper=1> p;
p=gamma[2];
qj = (1-p)*p^(J-1)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q2 = (1-J*p^(J-1)+(J-1)*p^J)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q0 = gamma_cdf(gamma[1],beta[1],beta[2])+p^J*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
}
model {

gamma[1] ~ exponential(.5);           // ??? gamma1 :small crater r 
gamma[2] ~ beta(2,2);                //prior ?? gamma2 :prob miss when s>r
beta[1] ~ gamma(1,.5);              //+
beta[2] ~ gamma(1,.5);              //
eta[1] ~ gamma(1,.5);              //+
eta[2] ~ gamma(1,.5); 

N2 ~ poisson(q2*rho);
for (i in 1:N2){
dat2s[i] ~ gamma(beta[1],beta[2]);   
dat2n[i] ~ binomial(J,1-p) T[2,];
}

N1 ~ poisson(qj*rho+rhoj);
for (i in 1:N1){
if(datj[i]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj)) + gamma_lpdf(datj[i] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj),
                  gamma_lpdf(datj[i] |beta[1],beta[2]) + log((1-p)*p^(J-1)) - log(qj),
                  gamma_lpdf(datj[i] |eta[1],eta[2]));
}}

target += -q0*rho+N0*log(q0*rho)-lgamma(N0+1);
}
"  
 
dat = list(J=10,N2=nrow(tryn2),N1=sum(dj),dat2s=tryn2[,1],dat2n=tryn2[,2],datj= tryj)
bfit = stan(model_code = crater.model, data = dat, iter = 1000, chains=4)
bfit
## para values for large simulation data J=10
#  rho = 2000 , rhoj=50
#  # ## missing prob =1 when size<3 / =0.2 when size>3
#  size pdf gamma(2.3,0.1)
## phantom size gamma(3,0.1)
###in simulation DATA: N=rpois(rho)=1998 #size<3=29=N0
plot(bfit)
e <- extract(bfit)
plot(e$gamma[,2],e$rho)
plot(1:length(e$gamma[,2]),e$gamma[,2],"l")
plot(1:length(e$rhoj),e$gamma[,1],"l")
plot(1:length(e$gamma[,2]),e$rho,"l")
plot(1:length(e$gamma[,2]),e$N0,"l")
plot(1:length(e$rho),e$beta[,2],"l")
plot(1:length(e$rho),e$beta[,1],"l")
plot(e$q2,e$rho)

quantile(e$q2*e$rho,c(0.025,0.975))


beta1=mean(e$beta[,1])
beta2=mean(e$beta[,2])
hist(c(tryn1[,1],tryn2[,1]),breaks=50,probability = T)
curve(dgamma(x,beta1,beta2),0,200,ylim=c(0,3),add=T,col=3)
curve(dgamma(x,2.3,.1),0,200,ylim=c(0,3),add=T,col=2)## red true para

bfit
J=10;p=0.2;gamma <- c(3,0.2);beta=c(2.3,0.1)
(1-p)*p^(J-1)*(1-pgamma(gamma[1],beta[1],beta[2])) ##qj
(1-J*p^(J-1)+(J-1)*p^J)*(1-pgamma(gamma[1],beta[1],beta[2])) ##q2 
pgamma(gamma[1],beta[1],beta[2])+p^J*(1-pgamma(gamma[1],beta[1],beta[2])) ##q0

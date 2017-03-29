#load("./nac_n1.Rdata")
#load("./nac_n1.Rdata")
nac_n1=read.table("./nac_n1.txt", header = T)
nac_n2=read.table("./nac_n2.txt", header = T)
nac_n2[nac_n2$nrecord>12,2]<-12  ## fix value when #repeats >12 #obs

(dj=table(nac_n1[,2]))

nacj=nac_n1[,1]

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

gamma[1] ~ gamma(3,1);              //  gamma1 :small crater r 
gamma[2] ~ beta(2,2);                // gamma2 :prob miss when s>r
beta[1] ~ gamma(5,.3);              //+
beta[2] ~ gamma(1,.5);              //
eta[1] ~ gamma(5,.3);              //+
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

dat = list(J=12,N2=nrow(nac_n2),N1=sum(dj),dat2s=nac_n2[,1],dat2n=nac_n2[,2],datj= nacj)
bfit = stan(model_code = crater.model, data = dat, iter = 1000, chains=4)
bfit

plot(bfit)
e <- extract(bfit)
load("./bfit6.RData")
load("./e6.RData")
bfit

library(ggplot2)
n0 = data.frame(N0=e$N0,gamma1=e$gamma[,1],rho=e$rho)
ggplot(data = n0,aes(gamma1,N0))+geom_point(alpha=0.3)+
  theme_classic()+labs(title=expression(paste("Posterior draws for ",N[0]," and ",gamma[1])),
                       x=expression(gamma[1]),y=expression(N[0])) +
  theme(plot.title = element_text(hjust = 0.5,size=18, face="bold"),
        axis.text=element_text(size=12),
        axis.title = element_text(color="black", face="bold", size=16))

ggplot(data = n0,aes(gamma1,rho))+geom_point(alpha=0.3)+
  theme_classic()+labs(title=expression(paste("Posterior draws for ",rho," and ",gamma[1])),
                       x=expression(gamma[1]),y=expression(rho)) +
  theme(plot.title = element_text(hjust = 0.5,size=18, face="bold"),
        axis.title.y=element_text(angle=0,vjust = 0.5),
        axis.text=element_text(size=12),
        axis.title = element_text(color="black", face="bold", size=16))
length=length(e$rho)
plot(1:length,e$gamma[,2],"l")
plot(1:length,e$gamma[,1],"l")
plot(1:length,e$rho,"l")
plot(1:length,e$rhoj,"l")
plot(1:length,e$N0,"l")
plot(1:length,e$beta[,2],"l")
plot(1:length,e$beta[,1],"l")

quantile(e$q2*e$rho,c(0.025,0.975))


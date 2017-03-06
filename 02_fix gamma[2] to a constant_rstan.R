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
real<lower=0> p;
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
real<lower=0> qj;
real<lower=0> q2;
real<lower=0> q0;

qj = (1-p)*p^(J-1)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q2 = (1-J*p^(J-1)+(J-1)*p^J)*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
q0 = gamma_cdf(gamma[1],beta[1],beta[2])+p^J*(1-gamma_cdf(gamma[1],beta[1],beta[2]));
}
model {

gamma[1] ~ exponential(10);           // ??? gamma1 :small crater r 
gamma[2] ~ beta(5,5);             //prior ?? gamma2 :prob miss when s>r
beta[1] ~ gamma(10,10);              //+
beta[2] ~ gamma(10,10);             //
eta[1] ~ gamma(10,10);              //+
eta[2] ~ gamma(10,10); 

N2 ~ poisson(q2*rho);
for (i in 1:N2){
dat2s[i] ~ gamma(beta[1],beta[2]);   
dat2n[i] ~ binomial(J,1-p) T[2,];
}

N1 ~ poisson(qj*rho+rhoj);
for (i in 1:k[1]){
if(datj1[i]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[1])) + gamma_lpdf(datj1[i] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[1]),
gamma_lpdf(datj1[i] |beta[1],beta[2]) + log((1-p)*p^(J-1)) - log(qj),
gamma_lpdf(datj1[i] |eta[1],eta[2]));
}}

for (j in 1:k[2]){
if(datj2[j]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[2])) + gamma_lpdf(datj2[j] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[2]),
gamma_lpdf(datj2[j] |beta[1],beta[2]) + log((1-p)*p^(J-1)) - log(qj),
gamma_lpdf(datj2[j] |eta[1],eta[2]));
}}
for (n in 1:k[3]){
if(datj3[n]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[3])) + gamma_lpdf(datj3[n] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[3]),
gamma_lpdf(datj3[n] |beta[1],beta[2]) + log((1-p)*p^(J-1)) - log(qj),
gamma_lpdf(datj3[n] |eta[1],eta[2]));
}}

target += -q0*rho+N0*log(q0*rho)-lgamma(N0+1);
}
"  

dat = list(p=0.2,J=3,N2=nrow(tryn2),N1=sum(dj),k=dj,dat2s=tryn2[,1],dat2n=tryn2[,2],datj1= tryj1,datj2= tryj2,datj3= tryj3)
fit = stan(model_code = crater.model, data = dat, iter = 1000, chains=6) #p=0.2
fit
dat1 = list(p=0.25,J=3,N2=nrow(tryn2),N1=sum(dj),k=dj,dat2s=tryn2[,1],dat2n=tryn2[,2],datj1= tryj1,datj2= tryj2,datj3= tryj3)
fit1 = stan(model_code = crater.model, data = dat1, iter = 1000, chains=4) #p=0.25
fit1
dat2 = list(p=0.3,J=3,N2=nrow(tryn2),N1=sum(dj),k=dj,dat2s=tryn2[,1],dat2n=tryn2[,2],datj1= tryj1,datj2= tryj2,datj3= tryj3)
fit2 = stan(model_code = crater.model, data = dat2, iter = 1000, chains=4) #p=0.3
fit2
dat3 = list(p=0.18,J=3,N2=nrow(tryn2),N1=sum(dj),k=dj,dat2s=tryn2[,1],dat2n=tryn2[,2],datj1= tryj1,datj2= tryj2,datj3= tryj3)
fit3 = stan(model_code = crater.model, data = dat3, iter = 1000, chains=4) #p=0.18
fit3
## para values
#  rho = 200 , rhoj=10
#  #gamma <- c(0.05,0.2)
#  size pdf gamma(1,1)
plot(fit)
e <- extract(fit)
plot(e$gamma[,2],e$rho)
plot(1:length(e$rho),e$gamma[,1],"l")
plot(1:length(e$rho),e$gamma[,2],"l")
plot(1:length(e$rho),e$rho,"l")
plot(1:length(e$rho),e$N0,"l")
plot(1:length(e$rho),e$beta[,2],"l")
plot(1:length(e$rho),e$beta[,1],"l")
plot(e$q2,e$rho)

quantile(e$q2*e$rho,c(0.025,0.975))

hist(e$gamma[,2],breaks=50,probability = T,xlim=c(0,1))
curve(dbeta(x,5,5),0,1,add=T)

beta1=mean(e$beta[,1])
beta2=mean(e$beta[,2])
hist(c(tryn1[,1],tryn2[,1]),breaks=50,probability = T)
curve(dgamma(x,beta1,beta2),0,15,ylim=c(0,3),add=T,col=3)
curve(dgamma(x,1,1),0,15,ylim=c(0,3),add=T,col=2)## red true para
curve(dgamma(x,1.2,1),add=T)   ## if change beta[2]=1
p=0.2
p=0.65
(1-J*p^(J-1)+(J-1)*p^J)*(1-pgamma(0.05,1,1)) ##highly depend on gamma[2]
fit

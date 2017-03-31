## for small dataset
load("/Users/zhang/Desktop/crater analysis/simdata-1.Rdata")
load("/Users/zhang/Desktop/crater analysis/simdata-2.Rdata")

#tryn2 <- tryx[tryx$X1!=1,]
#trys1=tryx[tryx[,2]==1,][,1]
#tryn1=try_d[,3:4][try_d[,3] %in% as.character(trys1),]
d1=sum(tryn1[,2]==1)
d2=sum(tryn1[,2]==2)
d3=sum(tryn1[,2]==3)
tryj1=tryn1[1:d1,1]
tryj2=tryn1[1:d2+d1,1]
tryj3=tryn1[1:d3+d2+d1,1]
dj=c(d1,d2,d3)

(dimj=table(tryn1[,2]))
#----------------------------------------
J=3


library(rstan)

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

gamma[1] ~ exponential(10);           // ??? gamma1 :small crater r 
gamma[2] ~ uniform(0,1);            //prior ?? gamma2 :prob miss when s>r
beta[1] ~ gamma(1,.1);              //+
beta[2] ~ gamma(1,.1);             //
eta[1] ~ gamma(1,.1);              //+
eta[2] ~ gamma(1,.1); 

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
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[1])) + gamma_lpdf(datj1[i] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[1]),
                  gamma_lpdf(datj1[i] |beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                  gamma_lpdf(datj1[i] |eta[1],eta[2]));
}}

for (j in 1:k[2]){
if(datj2[j]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[2])) + gamma_lpdf(datj2[j] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[2]),
                  gamma_lpdf(datj2[j] |beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                  gamma_lpdf(datj2[j] |eta[1],eta[2]));
}}
for (n in 1:k[3]){
if(datj3[n]<gamma[1]){
target += bernoulli_lpmf(0 |qj*rho/(qj*rho+rhoj[3])) + gamma_lpdf(datj3[n] |eta[1],eta[2]);
} else { 
target += log_mix(qj*rho/(qj*rho+rhoj[3]),
                 gamma_lpdf(datj3[n] |beta[1],beta[2]) + log(1-p)+(J-1)*log(p) - log(qj),
                 gamma_lpdf(datj3[n] |eta[1],eta[2]));
}}

target += -q0*rho+N0*log(q0*rho)-lgamma(N0+1);
}
"  
 
dat = list(J=3,N2=nrow(tryn2),N1=sum(dj),k=dj,dat2s=tryn2[,1],dat2n=tryn2[,2],datj1= tryj1,datj2= tryj2,datj3= tryj3)
fit = stan(model_code = crater.model, data = dat, iter = 1000, chains=3)
fit
##para values
#rho = 200 , rhoj=10
#size pdf gamma(1,1) 
#gamma <- c(0.05,0.2)
#miss 10 (<gamma1) + (38 34 35)randomly miss
## phantom num 8 10 7
plot(fit)
e <- extract(fit)
plot(x=1:length(e$rho),e$gamma[,2],"l")
plot(x=1:length(e$rho),e$rho,"l")
plot(x=e$gamma[,2],y=e$rho)
plot(x=e$q2,y=e$rho)
plot(x=e$N0,y=e$rho)
plot(x=e$N0,y=e$q2)
plot(x=log(e$gamma[,1]),y=e$q2)
### 1-integrate(function(x) dgamma(x,1,1),0.05,Inf)$value  1-F(x)

p=0.2
p=0.55
(1-J*p^(J-1)+(J-1)*p^J)*(1-pgamma(0.05,1,1)) #q2 !!!check q2
J*(1-p)*p^(J-1)*(1-pgamma(0.05,1,1)) #qj  ###add J*  !!!!
pgamma(0.05,1,1)+p^J*(1-pgamma(0.05,1,1)) #q0


quantile(e$q2*e$rho,c(0.025,0.975))

beta1=mean(e$beta[,1])
beta2=mean(e$beta[,2])
curve(dgamma(x,beta1,beta2),0,5,ylim=c(0,1),n=400)
curve(dgamma(x,1,1),0,5,col=2,add=T)
hist(simdata$size,freq = F,breaks=50,add=T)  ##!!!check beta2
curve(dgamma(x,1.2,1),0,5,col=3,add=T,n=200)

# ------------------------------
# Lídia André - January 2020
# ------------------------------

install.packages("rjags")
install.packages("runjags")
install.packages("coda")
library(rjags)
library(runjags)
library(coda)

# =====
# Data
# =====

castelo<-read.table("Cod_1200570 CASTELO_BRANCO.txt", header=T)

# ===================
# Data Preprocessing
# ===================

s<-function(dados){
  subset(dados,dados$year>=2000 & dados$year<=2012)
}

castelo2<-s(castelo)

autumn<-function(dados){
  a<-dados[(dados$month==9 & dados$day >= 22) | dados$month==10 | dados$month==11 | (dados$month==12 & dados$day<22),]
  n1<-dim(a)[[1]]
  a1<-seq(1,n1,5)
  a2<-a[a1,]
  return(a2)
}
winter<-function(dados){
  a<-dados[(dados$month==12 & dados$day >= 22) | dados$month==1 | dados$month==2 | (dados$month==3 & dados$day<21),]
  n1<-dim(a)[[1]]
  a1<-seq(1,n1,5)
  a2<-a[a1,]
  return(a2)
}
spring<-function(dados){
  a<-dados[(dados$month==3 & dados$day >= 21) | dados$month==4 | dados$month==5 | (dados$month==6 & dados$day<22),] 
  n1<-dim(a)[[1]]
  a1<-seq(1,n1,5)
  a2<-a[a1,]
  return(a2)
}
summer<-function(dados){
  a<-dados[(dados$month==6 & dados$day >= 22) | dados$month==7 | dados$month==8 | (dados$month==9 & dados$day<22),]
  n1<-dim(a)[[1]]
  a1<-seq(1,n1,5)
  a2<-a[a1,]
  return(a2)
}

castelo_a<-autumn(castelo2)
castelo_w<-winter(castelo2)
castelo_sp<-spring(castelo2)
castelo_su<-summer(castelo2)

# =================
# Bayesian Apprach
# =================

dados<-function(data){
  n<-dim(data)[[1]]
  x<-data$obs
  y<-data$sim
  zeros<-rep(0,n)
  return(list("zeros"=zeros,"n"=n,"x"=x,"y"=y))
}

# ---------------
# Castelo Branco
# ---------------

data<-dados(castelo_a)
inits1<-list(alpha=2.0,beta=0.5,mu=1.5,sigma=0.3)
inits2<-list(alpha=5.0,beta=1.5,mu=2.0,sigma=0.7)
inits<-list(inits1,inits2)

inits1_c<-list(rho=0.5)
inits2_c<-list(rho=0.2)
inits_c<-list(inits1_c,inits2_c)

marginals<-"model{
          for (i  in 1:n){
            x[i]~dgamma(alpha,beta)
            y[i]~dlnorm(mu,prec)
          }
          prec<-1/sigma^2
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            u[i]<-pgamma(x[i],alpha,beta)
            v[i]<-plnorm(y[i],mu,prec)
            l1[i]<-alpha*log(beta)+(alpha-1)*log(x[i])-beta*x[i]-loggam(alpha)
            l2[i]<--log(y[i]*sigma*sqrt(2*3.141593))-0.5*((log(y[i])-mu)^2/sigma^2)
            logL[i]<-l1[i]+l2[i]
            phi[i]<--logL[i]+C
          }
          alpha~dgamma(0.001,0.001)
          beta~dgamma(0.001,0.001)
          mu~dnorm(0.0,0.0001)
          sigma~dgamma(0.001,0.001)
        }"
marginals_a<-autorun.jags(marginals,monitor=c("alpha","beta","mu","sigma"), data=data,n.chains=2, inits=inits,startsample=4000,thin=30)
summary(marginals_a)
plot(marginals_a)

gaussian<-"model{
          for (i  in 1:n){
            x[i]~dgamma(alpha,beta)
            y[i]~dlnorm(mu,prec)
          }
          alpha<-*         # *substitute with the values of the command summary(marginals_a)[x,4] because BUGS language do not recognise R language
          beta<-*          #  x -> row which corresponds to the parameter in question
          mu<-*
          prec<-*^2
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            s[i]<-qnorm(u[i],0,1)
            t[i]<-(log(y[i])-mu)/sigma
            l3[i]<--0.5*log(1-rho^2)-0.5*(rho^2*s[i]^2+rho^2*t[i]^2-2*rho*s[i]*t[i])/(1-rho^2)
            logL[i]<-l3[i]
            phi[i]<--logL[i]+C
          }
          rho~dunif(-1.0,1.0)
        }"
gaussian_a<-autorun.jags(gaussian,monitor=c("rho"), data=data,n.chains=2, inits=inits_c,startsample=4000,thin=30)
summary(gaussian_a)
plot(gaussian_a)

# -----

data<-dados(castelo_w)
inits1<-list(mu1=1.5,sigma1=1.5,mu2=1.5,sigma2=0.3)
inits2<-list(mu1=2.0,sigma1=0.7,mu2=2.0,sigma2=0.7)
inits<-list(inits1,inits2)

inits1_c<-list(theta=0.5)
inits2_c<-list(theta=0.2)
inits_c<-list(inits1_c,inits2_c)

marginals<-"model{
          for (i  in 1:n){
            x[i]~dlnorm(mu1,prec1)
            y[i]~dlnorm(mu2,prec2)
          }
          prec<-1/sigma^2
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            u[i]<-plnorm(x[i],mu1,prec1)
            v[i]<-plnorm(y[i],mu2,prec2)
            l1[i]<--log(x[i]*sigma1*sqrt(2*3.141593))-0.5*((log(x[i])-mu1)^2/sigma1^2)
            l2[i]<--log(y[i]*sigma2*sqrt(2*3.141593))-0.5*((log(y[i])-mu2)^2/sigma2^2)
            logL[i]<-l1[i]+l2[i]
            phi[i]<--logL[i]+C
          }
          mu1~dnorm(0.0,0.0001)
          sigma1~dgamma(0.001,0.001)
          mu2~dnorm(0.0,0.0001)
          sigma2~dgamma(0.001,0.001)
        }"
marginals_w<-autorun.jags(marginals,monitor=c("mu1","sigma1","mu2","sigma2"), data=data,n.chains=2, inits=inits,startsample=4000,thin=30)
summary(marginals_w)
plot(marginals_w)

gumbel<-"model{
          for (i  in 1:n){
            x[i]~dlnorm(mu1,prec1)
            y[i]~dlnorm(mu2,prec2)
          }
          mu1<-*          # *substitute with the values of the command summary(marginals_w)[x,4] because BUGS language do not recognise R language
          mu2<-*          #  x -> row which corresponds to the parameter in question
          prec1<-*^2
          prec2<-*^2
          alphac<-1/theta
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            a[i]<--log(u[i])
            b[i]<--log(v[i])
            w[i]<-a[i]^(alphac)+b[i]^(alphac)
            l3[i]<--log(u[i]*v[i])+(alphac-1)*log(a[i]*b[i])+log(w[i]^(2*(1/alphac)-2)+(alphac-1)*w[i]^((1/alphac)-2))-w[i]^(1/alphac)
            logL[i]<-l3[i]
            phi[i]<--logL[i]+C
          }
          theta~dbeta(0.5,0.5)
        }"
gumbel_w<-autorun.jags(gumbel,monitor=c("alphac"), data=data,n.chains=2, inits=inits_c,startsample=4000)
summary(gumbel_w)
plot(gumbel_w)

# -----

data<-dados(castelo_sp)
inits1<-list(k=1.0,c=5.0,lambda=0.1,alpha=1.5,beta=0.3)
inits2<-list(k=6.0,c=1.0,lambda=0.5,alpha=2.5,beta=0.7)
inits<-list(inits1,inits2)

inits1_c<-list(rho=0.5)
inits2_c<-list(rho=0.2)
inits_c<-list(inits1_c,inits2_c)

marginals<-"model{
          for (i  in 1:n){
            x[i]~dpar4(k,w,0,gamma)
            y[i]~dgamma(alpha,beta)
          }
          w<-1/lambda
          gamma<-1/c
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            u[i]<-ppar4(x[i],k,w,0,gamma)
            v[i]<-pgamma(y[i],alpha,beta)
            l1[i]<-log(c)+log(k)+log(lambda)+(c-1)*(log(x[i])+log(lambda))-(k+1)*log(1+(x[i]*lambda)^c)
            l2[i]<-alpha*log(beta)+(alpha-1)*log(y[i])-beta*y[i]-loggam(alpha)
            logL[i]<-l1[i]+l2[i]
            phi[i]<--logL[i]+C
          }
          c~dgamma(0.001,0.001)
          k~dgamma(0.001,0.001)
          lambda~dgamma(0.001,0.001)
          alpha~dnorm(0.0,0.0001)
          beta~dgamma(0.001,0.001)
        }"
marginals_sp<-autorun.jags(marginals,monitor=c("k","c","lambda","alpha","beta"), data=data,n.chains=2, inits=inits,startsample=4000,thin=30)
summary(marginals_sp)
plot(marginals_sp)

gaussian<-"model{
            for (i in 1:n){
              x[i]~dpar4(k,w,0,gamma)
              y[i]~dgamma(alpha,beta)
            }
            k<-*          # *substitute with the values of the command summary(marginals_sp)[x,4] because BUGS language do not recognise R language
            w<-1/*        #  x -> row which corresponds to the parameter in question
            gamma<-1/*
            alpha<-*
            beta<-*
            C<-10000000
            for(i in 1:n){
              zeros[i]~dpois(phi[i])
              s[i]<-qnorm(u[i],0,1)
              t[i]<-qnorm(v[i],0,1)
              l3[i]<--0.5*log(1-rho^2)-0.5*(rho^2*s[i]^2+rho^2*t[i]^2-2*rho*s[i]*t[i])/(1-rho^2)
              logL[i]<-l3[i]
              phi[i]<--logL[i]+C
            }
            rho~dunif(-1.0,1.0)
          }"
gaussian_sp<-autorun.jags(gaussian,monitor=c("rho"), data=data,n.chains=2, inits=inits_c,startsample=4000,thin=30)
summary(gaussian_sp)
plot(gaussian_sp)

# -----

data<-dados(castelo_su)
inits1<-list(alpha=1.5,beta=0.3,omega=1.5,b=0.001)
inits2<-list(alpha=2.5,beta=0.7,omega=1.0,b=0.005)
inits<-list(inits1,inits2)

inits1_c<-list(rho=0.5)
inits2_c<-list(rho=0.2)
inits_c<-list(inits1_c,inits2_c)

marginals<-"model{
          for (i  in 1:n){
            x[i]~dgamma(alpha,beta)
            y[i]~dweib(omega,b)
          }
          delta<-b^(-1/omega)
          C<-10000000
          for(i in 1:n){
            zeros[i]~dpois(phi[i])
            u[i]<-pgamma(x[i],alpha,beta)
            v[i]<-pweib(y[i],omega,b)
            l1[i]<-alpha*log(beta)+(alpha-1)*log(x[i])-beta*x[i]-loggam(alpha)
            l2[i]<-log(omega)+log(b)+(omega-1)*(log(y[i]))-b*(y[i])^omega
            logL[i]<-l1[i]+l2[i]
            phi[i]<--logL[i]+C
          }
          alpha~dnorm(0.0,0.0001)
          beta~dgamma(0.001,0.001)
          omega~dgamma(0.001,0.001)
          b~dgamma(0.001,0.001)
        }"
marginals_su<-autorun.jags(marginals,monitor=c("alpha","beta","omega","delta"), data=data,n.chains=2, inits=inits,startsample=4000,thin=30)
summary(marginals_su)
plot(marginals_su)

gaussian<-"model{
            for (i in 1:n){
              x[i]~dgamma(alpha,beta)
              y[i]~dweib(omega,b)
            }
            alpha<-*          # *substitute with the values of the command summary(marginals_sp)[x,4] because BUGS language do not recognise R language
            beta<-*           #  x -> row which corresponds to the parameter in question
            omega<-*
            b<-omega^*        # in this case substitute both omega and delta (the * in this row) with the corresponded values obtained with marginals_su
            C<-10000000
            for(i in 1:n){
              zeros[i]~dpois(phi[i])
              s[i]<-qnorm(u[i],0,1)
              t[i]<-qnorm(v[i],0,1)
              l3[i]<--0.5*log(1-rho^2)-0.5*(rho^2*s[i]^2+rho^2*t[i]^2-2*rho*s[i]*t[i])/(1-rho^2)
              logL[i]<-l3[i]
              phi[i]<--logL[i]+C
            }
            rho~dunif(-1.0,1.0)
          }"
gaussian_su<-autorun.jags(gaussian,monitor=c("rho"), data=data,n.chains=2, inits=inits_c,startsample=4000,thin=30)
summary(gaussian_su)
plot(gaussian_su)


# ------------------------------
# Lídia André - September 2019
# ------------------------------

install.packages("copula")
install.packages("VineCopula")
install.packages("MASS")
install.packages("fitdistrplus")
install.packages("scatterplot3d")
install.packages("gofCopula")
install.packages("dplyr")
install.packages("actuar")
install.packages("extraDistr")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("car")
install.packages("EnvStats")
install.packages("gridExtra")
install.packages("KScorrect")
install.packages("goft")
install.packages("latticeExtra")
install.packages("rjags")
install.packages("runjags")
install.packages("coda")
library(copula)
library(VineCopula)
library(MASS)
library(fitdistrplus)
library(scatterplot3d)
library(gofCopula)
library(dplyr)
library(actuar)
library(extraDistr)
library(ggplot2)
library(RColorBrewer)
library(car)
library(EnvStats)
library(gridExtra)
library(KScorrect)
library(goft)
library(latticeExtra)
library(rjags)
library(runjags)
library(coda)

# ===================
# Data Preprocessing
# ===================

s<-function(dados){
  subset(dados,dados$year>=2000 & dados$year<=2012)
}

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

# ========
# Boxplot
# ========

boxplot<-function(a,w,sp,su,titulo1,titulo2){
  n_a<-dim(a)[[1]]
  n_w<-dim(w)[[1]]
  n_sp<-dim(sp)[[1]]
  n_su<-dim(su)[[1]]
  dat1<-data.frame(Season<-factor(c(rep("Autumn",n_a),rep("Winter",n_w),rep("Spring",n_sp),rep("Summer",n_su))),Wind=c(a$obs,w$obs,sp$obs,su$obs))
  dat2<-data.frame(Season<-factor(c(rep("Autumn",n_a),rep("Winter",n_w),rep("Spring",n_sp),rep("Summer",n_su))),Wind=c(a$sim,w$sim,sp$sim,su$sim))
  o<-ggplot(dat1, aes(x=Season, y=Wind, color=Season,fill=Season)) + geom_boxplot(alpha=0.5)+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+coord_cartesian(ylim=c(0,20))+ggtitle(titulo1)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  s<-ggplot(dat2, aes(x=Season, y=Wind, color=Season,fill=Season)) + geom_boxplot(alpha=0.5)+scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+coord_cartesian(ylim=c(0,20))+ggtitle(titulo2)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")
  grid.arrange(o, s, ncol=2)
}

# ====
# ACF
# ====

autocorr<-function(a,w,sp,su,titulo1,titulo2,titulo3,titulo4,titulo5,titulo6,titulo7,titulo8){
  par(mfrow=c(2,4))
  acf(a$obs, main=titulo1)
  acf(w$obs, main=titulo2)
  acf(sp$obs, main=titulo3)
  acf(su$obs, main=titulo4)
  acf(a$sim, main=titulo5)
  acf(w$sim, main=titulo6)
  acf(sp$sim, main=titulo7)
  acf(su$sim, main=titulo8)
}

# ============================
# Marginal Distributions: GoF
# ============================

fit<-function(data,dist){
  x<-data$obs
  y<-data$sim
  if(dist=="burr"){
    f_x<-fitdist(x,dist,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    f_y<-fitdist(y,dist,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  }
  else{
    f_x<-fitdist(x,dist)
    f_y<-fitdist(y,dist)
  }
  estimate_x<-round(f_x$estimate,4)
  estimate_y<-round(f_y$estimate,4)
  aic_x<-round(f_x$aic,4)
  aic_y<-round(f_y$aic,4)
  return(c(estimate_x,aic_x,estimate_y,aic_y))
}

gof<-function(data,fun){
  x<-data$obs
  y<-data$sim
  if(fun=="lnorm"){
    f_x<-fitdist(x,fun)
    f_y<-fitdist(y,fun)
    crit<-round(c(.8, .85, .9, .95, .99) * 4999, 0)
    Do_x<-LcKS(x,"plnorm")
    Ds_x<-sort(Do_x$D.sim)
    critic_x<-round(Ds_x[crit], 3)
    Do_y<-LcKS(y,"plnorm")
    Ds_y<-sort(Do_y$D.sim)
    critic_y<-round(Ds_y[crit], 3)
    ks<-data.frame(Do_x$D.obs<critic_x[5],Do_y$D.obs<critic_y[5])
    gof_x<-data.frame(round(gofstat(f_x)$chisqpvalue,4))
    gof_y<-data.frame(round(gofstat(f_y)$chisqpvalue,4))
    a_x<-lnorm_test(x)
    a_y<-lnorm_test(y)
    shapiro_x<-a_x$p.value>0.05
    shapiro_y<-a_y$p.value>0.05
    return(c(gof_x,gof_y,ks,shapiro_x,shapiro_y))
  }
  else if(fun=="gamma"){
    f_x<-fitdist(x,fun)
    f_y<-fitdist(y,fun)
    crit<-round(c(.8, .85, .9, .95, .99) * 4999, 0)
    Do_x<-LcKS(x,"pgamma")
    Ds_x<-sort(Do_x$D.sim)
    critic_x<-round(Ds_x[crit], 3)
    Do_y<-LcKS(y,"pgamma")
    Ds_y<-sort(Do_y$D.sim)
    critic_y<-round(Ds_y[crit], 3)
    ks<-data.frame(Do_x$D.obs<critic_x[5],Do_y$D.obs<critic_y[5])
    gof_x<-data.frame(round(gofstat(f_x)$chisqpvalue,4),gofstat(f_x)$cvmtest,gofstat(f_x)$adtest)
    gof_y<-data.frame(round(gofstat(f_y)$chisqpvalue,4),gofstat(f_y)$cvmtest,gofstat(f_y)$adtest)
    return(c(gof_x,gof_y,ks))
  }
  else if(fun=="weibull"){
    f_x<-fitdist(x,fun)
    f_y<-fitdist(y,fun)
    crit<-round(c(.8, .85, .9, .95, .99) * 4999, 0)
    Do_x<-LcKS(x,"pweibull")
    Ds_x<-sort(Do_x$D.sim)
    critic_x<-round(Ds_x[crit], 3)
    Do_y<-LcKS(y,"pweibull")
    Ds_y<-sort(Do_y$D.sim)
    critic_y<-round(Ds_y[crit], 3)
    ks<-data.frame(Do_x$D.obs<critic_x[5],Do_y$D.obs<critic_y[5])
    gof_x<-data.frame(round(gofstat(f_x)$chisqpvalue,4),gofstat(f_x)$cvmtest,gofstat(f_x)$adtest)
    gof_y<-data.frame(round(gofstat(f_y)$chisqpvalue,4),gofstat(f_y)$cvmtest,gofstat(f_y)$adtest)
    return(c(gof_x,gof_y,ks))
  }
  else{
    f_x<-fitdist(x,fun,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    f_y<-fitdist(y,fun,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    gof_x<-data.frame(round(gofstat(f_x)$chisqpvalue,4))
    gof_y<-data.frame(round(gofstat(f_y)$chisqpvalue,4))
    return(c(gof_x,gof_y)) 
  }
}

Aic<-function(a,b,c,d){
  x<-data.frame(a[3],b[3],c[3],d[4])
  y<-data.frame(a[6],b[6],c[6],d[8])
  min_x<-which.min(x)
  min_y<-which.min(y)
  return(c(min_x,min_y))
}

rsquared<-function(data,fun){
  n<-dim(data)[[1]]
  x<-data$obs
  y<-data$sim
  p<-(1:n)/(n+1)
  if(fun=="lnorm"){
    qq<-qnorm(p)
    sort_x<-sort(log(x))
    sort_y<-sort(log(y))
    r_x<-round(cor(qq,sort_x)^2,4)
    r_y<-round(cor(qq,sort_y)^2,4)
    return(c(r_x,r_y))
  }
  else if(fun=="weibull"){
    qq<-log(-log(1-p))
    sort_x<-sort(log(x))
    sort_y<-sort(log(y))
    r_x<-round(cor(qq,sort_x)^2,4)
    r_y<-round(cor(qq,sort_y)^2,4)
    return(c(r_x,r_y))
  }
  else if(fun=="gamma"){
    f1<-fitdist(x,fun)
    f2<-fitdist(y,fun)
    p_i<--log(1-p)
    sort_x<-sort(x)
    sort_y<-sort(y)
    f_x<-pgamma(sort_x,shape=f1$estimate[1],rate=f1$estimate[2])
    f_y<-pgamma(sort_y,shape=f2$estimate[1],rate=f2$estimate[2])
    e_x<--log(1-f_x)
    e_y<--log(1-f_y)
    r_x<-round(cor(p_i,e_x)^2,4)
    r_y<-round(cor(p_i,e_y)^2,4)
    return(c(r_x,r_y))
  }
  else{
    f1<-fitdist(x,fun,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    f2<-fitdist(y,fun,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    p_i<--log(1-p)
    sort_x<-sort(x)
    sort_y<-sort(y)
    f_x<-pburr(sort_x,shape1=f1$estimate[1],shape2=f1$estimate[2],rate=f1$estimate[3])
    f_y<-pburr(sort_y,shape1=f2$estimate[1],shape2=f2$estimate[2],rate=f2$estimate[3])
    e_x<--log(1-f_x)
    e_y<--log(1-f_y)
    r_x<-round(cor(p_i,e_x)^2,4)
    r_y<-round(cor(p_i,e_y)^2,4)
    return(c(r_x,r_y)) 
  }
}

comp<-function(data,fun1,fun2,fun3,fun4,titulo){
  x<-data$obs
  y<-data$sim
  f1_x<-fitdist(x,fun1)
  f1_y<-fitdist(y,fun1)
  f2_x<-fitdist(x,fun2)
  f2_y<-fitdist(y,fun2)
  f3_x<-fitdist(x,fun3)
  f3_y<-fitdist(y,fun3)
  f4_x<-fitdist(x,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  f4_y<-fitdist(y,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  par(mfrow=c(2,4))
  plot.legend<-c("Lognormal", "Gamma", "Weibull","Burr")
  denscomp(list(f1_x,f2_x,f3_x,f4_x), legendtext=plot.legend,main="Observed Wind")
  qqcomp(list(f1_x,f2_x,f3_x,f4_x), legendtext=plot.legend)
  cdfcomp(list(f1_x,f2_x,f3_x,f4_x), legendtext=plot.legend)
  ppcomp(list(f1_x,f2_x,f3_x,f4_x), legendtext=plot.legend)
  mtext(titulo, side = 3, line = -1.5, outer = TRUE)
  plot.legend<-c("Lognormal", "Gamma", "Weibull","Burr")
  denscomp(list(f1_y,f2_y,f3_y,f4_y), legendtext=plot.legend,main="Simulated Wind")
  qqcomp(list(f1_y,f2_y,f3_y,f4_y), legendtext=plot.legend)
  cdfcomp(list(f1_y,f2_y,f3_y,f4_y), legendtext=plot.legend)
  ppcomp(list(f1_y,f2_y,f3_y,f4_y), legendtext=plot.legend)
  par(mfrow=c(1,1))
}

qq<-function(data,fun1,fun2,fun3,fun4,titulo){
  x<-data$obs
  y<-data$sim
  f1_x<-fitdist(x,fun1)
  f1_y<-fitdist(y,fun1)
  f2_x<-fitdist(x,fun2)
  f2_y<-fitdist(y,fun2)
  f3_x<-fitdist(x,fun3)
  f3_y<-fitdist(y,fun3)
  f4_x<-fitdist(x,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  f4_y<-fitdist(y,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  par(mfrow=c(2,4))
  car::qqPlot(x,distribution=fun1,meanlog=f1_x$estimate[1],sdlog=f1_x$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(x,distribution=fun2,shape=f2_x$estimate[1],rate=f2_x$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(x,distribution=fun3,shape=f3_x$estimate[1],scale=f3_x$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(x,distribution=fun4,shape1=f4_x$estimate[1],shape2=f4_x$estimate[2],rate=f4_x$estimate[3],id=F,pch=21,cex=0.7,lwd=1)
  mtext(titulo, side = 3, line = -1.5, outer = TRUE)
  car::qqPlot(y,distribution=fun1,meanlog=f1_y$estimate[1],sdlog=f1_y$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(y,distribution=fun2,shape=f2_y$estimate[1],rate=f2_y$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(y,distribution=fun3,shape=f3_y$estimate[1],scale=f3_y$estimate[2],id=F,pch=21,cex=0.7,lwd=1)
  car::qqPlot(y,distribution=fun4,shape1=f4_y$estimate[1],shape2=f4_y$estimate[2],rate=f4_y$estimate[3],id=F,pch=21,cex=0.7,lwd=1)
  par(mfrow=c(1,1))  
}

kernel<-function(data,fun1,fun2,fun3,fun4,titulo){
  x<-data$obs
  y<-data$sim
  f1_x<-fitdist(x,fun1)
  f1_y<-fitdist(y,fun1)
  f2_x<-fitdist(x,fun2)
  f2_y<-fitdist(y,fun2)
  f3_x<-fitdist(x,fun3)
  f3_y<-fitdist(y,fun3)
  f4_x<-fitdist(x,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  f4_y<-fitdist(y,fun4,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  par(mfrow=c(2,4))
  plot(density(x),main="Kernel vs Lognormal",xlab="Observed Wind",lwd=1,ylim=c(0,0.4))
  curve(dlnorm(x, meanlog=f1_x$estimate[1],sdlog=f1_x$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Lognormal'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(x),main="Kernel vs Gamma",xlab="Observed Wind",lwd=1,ylim=c(0,0.4))
  curve(dgamma(x, shape=f2_x$estimate[1],rate=f2_x$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Gamma'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(x),main="Kernel vs Weibull",xlab="Observed Wind",lwd=1,ylim=c(0,0.4))
  curve(dweibull(x, shape=f3_x$estimate[1],scale=f3_x$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Weibull'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(x),main="Kernel vs Burr",xlab="Observed Wind",lwd=1,ylim=c(0,0.4))
  curve(dburr(x, shape1=f4_x$estimate[1],shape2=f4_x$estimate[2],rate=f4_x$estimate[3]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Burr'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  mtext(titulo, side = 3, line = -1.5, outer = TRUE)
  plot(density(y),main="Kernel vs Lognormal",xlab="Simulated Wind",lwd=1,ylim=c(0,0.4))
  curve(dlnorm(x, meanlog=f1_y$estimate[1],sdlog=f1_y$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Lognormal'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(y),main="Kernel vs Gamma",xlab="Simulated Wind",lwd=1,ylim=c(0,0.4))
  curve(dgamma(x, shape=f2_y$estimate[1],rate=f2_y$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Gamma'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(y),main="Kernel vs Weibull",xlab="Simulated Wind",lwd=1,ylim=c(0,0.4))
  curve(dweibull(x, shape=f3_y$estimate[1],scale=f3_y$estimate[2]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Weibull'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  plot(density(y),main="Kernel vs Burr",xlab="Simulated Wind",lwd=1,ylim=c(0,0.4))
  curve(dburr(x, shape1=f4_y$estimate[1],shape2=f4_y$estimate[2],rate=f4_y$estimate[3]), col="red3",lty=3,lwd=2,add=T)
  legend('topright',c('Kernel','Burr'),col=c('black','red3'),lwd=c(1,2),lty=c(1,3),cex=0.7)
  par(mfrow=c(1,1))  
}

graficos<-function(data,fun1,fun2,fun3,fun4,titulo){
  comp(data,fun1,fun2,fun3,fun4,titulo)
  qq(data,fun1,fun2,fun3,fun4,titulo)
  kernel(data,fun1,fun2,fun3,fun4,titulo)
}

# ============================
# Marginal Distributions: Fit
# ============================

confidence<-function(data,dist){
  x<-data$obs
  y<-data$sim
  if(dist=="burr"){
    f_x<-fitdist(x,dist,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
    f_y<-fitdist(y,dist,start=list(shape1=0.1,shape2=0.1,rate=0.01)) # the start values have to be different in some of the cases
  }
  else{
    f_x<-fitdist(x,dist)
    f_y<-fitdist(y,dist)
  }
  a<-bootdist(f_x,niter=4e3)
  b<-bootdist(f_y,niter=4e3)
  ci_x<-summary(a)
  ci_y<-summary(b)
  return(c(ci_x,ci_y))
}

marginals<-function(data,fun1,fun2,titulo){
  n<-dim(data)[[1]]
  x<-data$obs
  y<-data$sim
  if(fun1=="burr" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
  }
  else if(fun1=="burr" || fun2=="burr"){
    if(fun1=="burr"){
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      c_x<-fit(data,fun1)[3]
      a_y<-fit(data,fun2)[4]
      b_y<-fit(data,fun2)[5]
    }
    else{
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      a_y<-fit(data,fun2)[5]
      b_y<-fit(data,fun2)[6]
      c_y<-fit(data,fun2)[7]
    }
  }
  else{
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
  }
  df<-data.frame(Wind=factor(rep(c("Observed","Simulated"),each=n)),Data=c(x,y))
  if(fun1=="lnorm" && fun2=="lnorm"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dlnorm(x, meanlog=a_x,sdlog=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dlnorm(x, meanlog=a_y,sdlog=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="lnorm" && fun2=="gamma"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dlnorm(x, meanlog=a_x,sdlog=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dgamma(x, shape=a_y,rate=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="lnorm" && fun2=="weibull"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dlnorm(x, meanlog=a_x,sdlog=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dweibull(x, shape=a_y,scale=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="lnorm" && fun2=="burr"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dlnorm(x, meanlog=a_x,sdlog=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dburr(x, shape1=a_y,shape2=b_y,rate=c_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="gamma" && fun2=="lnorm"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dgamma(x, shape=a_x,rate=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dlnorm(x, meanlog=a_y,sdlog=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="gamma" && fun2=="gamma"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dgamma(x, shape=a_x,rate=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dgamma(x, shape=a_y,rate=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="gamma" && fun2=="weibull"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dgamma(x, shape=a_x,rate=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dweibull(x, shape=a_y,scale=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="gamma" && fun2=="burr"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dgamma(x, shape=a_x,rate=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dburr(x, shape1=a_y,shape2=b_y,rate=c_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="weibull" && fun2=="lnorm"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dweibull(x, shape=a_x,scale=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dlnorm(x, meanlog=a_y,sdlog=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="weibull" && fun2=="gamma"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dweibull(x, shape=a_x,scale=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dgamma(x, shape=a_y,rate=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="weibull" && fun2=="weibull"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dweibull(x, shape=a_x,scale=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dweibull(x, shape=a_y,scale=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="weibull" && fun2=="burr"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dweibull(x, shape=a_x,scale=b_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dburr(x, shape1=a_y,shape2=b_y,rate=c_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="burr" && fun2=="lnorm"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dburr(x, shape1=a_x,shape2=b_x,rate=c_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dlnorm(x, meanlog=a_y,sdlog=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="burr" && fun2=="gamma"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dburr(x, shape1=a_x,shape2=b_x,rate=c_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dgamma(x, shape=a_y,rate=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else if(fun1=="burr" && fun2=="weibull"){
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dburr(x, shape1=a_x,shape2=b_x,rate=c_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dweibull(x, shape=a_y,scale=b_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
  else{
    ggplot(df,aes(x=Data,fill=Wind,color=Wind))+geom_histogram(alpha=0.5,bins=12,position="identity",aes(y=..density..))+scale_color_brewer(palette="Dark2")+
      scale_fill_brewer(palette="Dark2")+stat_function(fun = function(x) dburr(x, shape1=a_x,shape2=b_x,rate=c_x),color = "darkgreen", size = 0.7)+stat_function(fun = function(x) dburr(x, shape1=a_y,shape2=b_y,rate=c_y),color = "darkorange", size = 0.7)+coord_cartesian(ylim=c(0,0.4))+ggtitle(titulo)+theme(plot.title = element_text(hjust = 0.5))
  }
}

# ==================
# Copula: Selection
# ==================

copula<-function(data){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  u<-pobs(mat)[,1]
  v<-pobs(mat)[,2]
  fam<-BiCopSelect(u,v,familyset=NA,se=T)$family
  if(fam==2){
    par1<-BiCopSelect(u,v,familyset=fam,se=T)$par
    par2<-BiCopSelect(u,v,familyset=fam,se=T)$par2
    se1<-BiCopSelect(u,v,familyset=fam,se=T)$se
    se2<-BiCopSelect(u,v,familyset=fam,se=T)$se2
    return(c(fam,par1,se1,par2,se2))
  }
  else if(fam==7 || fam==8 || fam==9 || fam==10 || fam==17 || fam==18 || fam==19 || fam==20 || fam==27 || fam==28 || fam==29 || fam==30 || fam==37 || fam==38 || fam==39 || fam==40 || fam==104 || fam==114 || fam==124 || fam==134 || fam==204 || fam==214 || fam==224 || fam==234){
    aic1<-BiCopSelect(u,v,familyset=1,se=T)$AIC
    aic2<-BiCopSelect(u,v,familyset=2,se=T)$AIC
    aic3<-BiCopSelect(u,v,familyset=3,se=T)$AIC
    aic4<-BiCopSelect(u,v,familyset=4,se=T)$AIC
    aic5<-BiCopSelect(u,v,familyset=5,se=T)$AIC
    aic6<-BiCopSelect(u,v,familyset=6,se=T)$AIC
    aic<-data.frame(aic1,aic2,aic3,aic4,aic5,aic6)
    min<-which.min(aic)
    if(min==1){
      family<-BiCopSelect(u,v,familyset=1,se=T)$family
      par<-BiCopSelect(u,v,familyset=1,se=T)$par
      se<-BiCopSelect(u,v,familyset=1,se=T)$se
      return(c(family,par,se))
    }
    else if(min==2){
      family<-BiCopSelect(u,v,familyset=2,se=T)$family
      par1<-BiCopSelect(u,v,familyset=2,se=T)$par
      par2<-BiCopSelect(u,v,familyset=2,se=T)$par2
      se1<-BiCopSelect(u,v,familyset=2,se=T)$se
      se2<-BiCopSelect(u,v,familyset=2,se=T)$se2
      return(c(family,par1,se1,par2,se2))
    }
    else if(min==3){
      family<-BiCopSelect(u,v,familyset=3,se=T)$family
      par<-BiCopSelect(u,v,familyset=3,se=T)$par
      se<-BiCopSelect(u,v,familyset=3,se=T)$se
      return(c(family,par,se))
    }
    else if(min==4){
      family<-BiCopSelect(u,v,familyset=4,se=T)$family
      par<-BiCopSelect(u,v,familyset=4,se=T)$par
      se<-BiCopSelect(u,v,familyset=4,se=T)$se
      return(c(family,par,se))
    }
    else if(min==5){
      family<-BiCopSelect(u,v,familyset=5,se=T)$family
      par<-BiCopSelect(u,v,familyset=5,se=T)$par
      se<-BiCopSelect(u,v,familyset=5,se=T)$se
      return(c(family,par,se))
    }
    else{
      family<-BiCopSelect(u,v,familyset=6,se=T)$family
      par<-BiCopSelect(u,v,familyset=6,se=T)$par
      se<-BiCopSelect(u,v,familyset=6,se=T)$se
      return(c(family,par,se))
    }
  }
  else{
    par<-BiCopSelect(u,v,familyset=fam,se=T)$par
    se<-BiCopSelect(u,v,familyset=fam,se=T)$se
    return(c(fam,par,se))
  }
}

# ============
# Copula: Fit
# ============

fitcop<-function(data){
  cop<-copula(data)
  if(cop[1]==1){
    rho<-cop[2]
    gaus<-normalCopula(param=rho,dim=2)
  }
  else if(cop[1]==2){
    rho<-cop[2]
    eta<-cop[4]
    t<-tCopula(param=rho,df=eta,dim=2)
  }
  else if(cop[1]==3){
    alpha<-cop[2]
    c<-claytonCopula(param=alpha,dim=2)
  }
  else if(cop[1]==4){
    alpha<-cop[2]
    g<-gumbelCopula(param=alpha,dim=2)
  }
  else if(cop[1]==5){
    alpha<-cop[2]
    f<-frankCopula(param=alpha,dim=2)
  }
  else if(cop[1]==6){
    alpha<-cop[2]
    j<-joeCopula(param=alpha,dim=2)
  }
  else if(cop[1]==13){
    alpha<-cop[2]
    sc<-rotCopula(claytonCopula(param=alpha,dim=2))
  }
  else{
  	alpha<-cop[2]
  	sg<-rotCopula(gumbelCopula(param=alpha,dim=2))
  }
}

mple<-function(data){
	x<-data$obs
	y<-data$sim
	mat<-cbind(x,y)
	#cop<-fitcop(data)
	c<-copula(data)
	if(c[1]==1){
		est<-fitCopula(normalCopula(),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==2){
		est<-fitCopula(tCopula(),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==3){
		est<-fitCopula(claytonCopula(),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==13){
		est<-fitCopula(rotCopula(claytonCopula()),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==4){
		est<-fitCopula(gumbelCopula(),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==14){
		est<-fitCopula(rotCopula(gumbelCopula()),pobs(mat),method="mpl")
		return(est)
	}
	else if(c[1]==5){
		est<-fitCopula(frankCopula(),pobs(mat),method="mpl")
		return(est)
	}
	else{
		est<-fitCopula(joeCopula(),pobs(mat),method="mpl")
		return(est)
	}
}

int_mple<-function(data){
	a<-mple(data)
	return(confint(a))
}

pdf_plot<-function(data,titulo){
  cop<-copula(data)
  c<-fitcop(data)
  if(cop[1]==1){
    col="chartreuse2"
  }
  else if(cop[1]==2){
    col="chocolate1"
  }
  else if(cop[1]==3){
    col="orangered1"
  }
  else if(cop[1]==13){
    col="deepskyblue"
  }
  else if(cop[1]==4){
    col="skyblue1"
  }
  else if(cop[1]==14){
    col="red1"
  }
  else if(cop[1]==5){
    col="darkseagreen1"
  }
  else{
    col="lightskyblue1"
  }
  persp(c,dCopula,main=titulo,theta=45, phi = 20, col = col,ltheta = 120, shade = 0.8, ticktype = "detailed",xlab = "u", ylab = "v", zlab = "c(u,v)",zlim=c(0,14))
}

# ========================
# Copula: Tail Dependence
# ========================

tail<-function(data){
  lambda(fitcop(data))
}

# ============
# Copula: GOF
# ============

gofcop<-function(data,name){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  u<-pobs(mat)[,1]
  v<-pobs(mat)[,2]
  c<-copula(data)
  if(c[1]==2){
    rho<-c[2]
    eta<-c[4]
    t<-tCopula(param=rho,df=round(eta,0),df.fixed=T)
    kendallcvm<-gofKendallCvM(name,mat,param=rho,df=eta)
    rosemblatt<-gofRosenblattSnB(name,mat,param=rho,df=eta)
    sn<-gofCopula(t,pobs(mat))
    return(c(round(rosemblatt$erg.tests[1],4),round(sn$p.value,4),round(kendallcvm$erg.tests[1],4)))
  }
  else if(c[1]==6){
    alpha<-c[2]
    j<-joeCopula(param=alpha,dim=2)
    sn<-gofCopula(j,pobs(mat))
    kendallcvm<-BiCopGofTest(u,v,family=6,method="kendall")
    return(c(round(sn$p.value,4),round(kendallcvm$p.value.KS,4)))
  }
  else{
    rosemblatt<-gofRosenblattSnB(name,mat,param=c[2])
    sn<-gofSn(name,mat,param=c[2])
    kendallcvm<-gofKendallCvM(name,mat,param=c[2])
    return(c(round(rosemblatt$erg.tests[1],4),round(sn$erg.tests[1],4),round(kendallcvm$erg.tests[1],4)))
  }
}

# ===========================
# Copula: Estimation Methods
# ===========================

mle<-function(data,fun1,fun2){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  cop<-fitcop(data)
  c<-copula(data)
  if(fun1=="burr" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    mv<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y,c_y)))
    if(c[1]==1){
      st<-c(0.1,0.1,0.01,0.1,0.1,0.01,0.1)
    }
    else if(c[1]==2){
        st<-c(0.1,0.1,0.01,0.1,0.1,0.01,0.1,2)
    }
    else if(c[1]==3){
        st<-c(0.1,0.1,0.01,0.1,0.1,0.01,0.1)
    }
    else if(c[1]==4){
        st<-c(0.1,0.1,0.01,0.1,0.1,0.01,1.1)
    }
    else if(c[1]==5){
        st<-c(0.1,0.1,0.01,0.1,0.1,0.01,1)
    }
    else{
        st<-c(0.1,0.1,0.01,0.1,0.1,0.01,1.1)
    }
  }
  else if(fun1=="burr" || fun2=="burr"){
    if(fun1=="burr"){
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      c_x<-fit(data,fun1)[3]
      a_y<-fit(data,fun2)[5]
      b_y<-fit(data,fun2)[6]
      mv<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y)))
      if(c[1]==1){
        st<-c(0.01,0.1,0.01,1,1,0.1)
      }
      else if(c[1]==2){
        st<-c(0.01,0.1,0.01,1,1,0.1,2)
      }
      else if(c[1]==3){
        st<-c(0.01,0.1,0.01,1,1,0.1)
      }
      else if(c[1]==4){
        st<-c(0.01,0.1,0.01,1,1,1.1)
      }
      else if(c[1]==5){
        st<-c(0.01,0.1,0.01,1,1,1)
      }
      else{
        st<-c(0.01,0.1,0.01,1,1,1.1)
      }
    }
    else{
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      a_y<-fit(data,fun2)[4]
      b_y<-fit(data,fun2)[5]
      c_y<-fit(data,fun2)[6]
      mv<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y,c_y)))
      if(c[1]==1){
        st<-c(1,1,0.01,0.1,0.01,0.1)
      }
      else if(c[1]==2){
        st<-c(1,1,0.01,0.1,0.01,0.1,2)
      }
      else if(c[1]==3){
        st<-c(1,1,0.01,0.1,0.01,0.1)
      }
      else if(c[1]==4){
        st<-c(1,1,0.01,0.1,0.01,1.1)
      }
      else if(c[1]==5){
        st<-c(1,1,0.01,0.1,0.01,1)
      }
      else{
        st<-c(1,1,0.01,0.1,0.01,1.1)
      }
    }
  }
  else{
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    mv<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y)))
    if(c[1]==1){
      st<-c(1,1,1,1,0.1)
    }
    else if(c[1]==2){
      st<-c(1,1,1,1,0.1,2)
    }
    else if(c[1]==3){
      st<-c(1,1,1,1,0.1)
    }
    else if(c[1]==4){
      st<-c(1,1,1,1,1.1)
    }
    else if(c[1]==5){
      st<-c(1,1,1,1,1)
    }
    else{
      st<-c(1,1,1,1,1.1)
    }
  }
  est<-fitMvdc(mat,mv,st)
  return(est)
}

int_mle<-function(data,fun1,fun2){
  a<-mle(data,fun1,fun2)
  return(confint(a))
}

ifme<-function(data,fun1,fun2){
  x<-data$obs
  y<-data$sim
  c<-copula(data)
  if(fun1=="burr" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    d<-cbind(pburr(x,shape1=a_x,shape2=b_x,rate=c_x),pburr(y,shape1=a_y,shape2=b_y,rate=c_y))
  }
  else if(fun1=="burr" && fun2=="lnorm"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pburr(x,shape1=a_x,shape2=b_x,rate=c_x),plnorm(y,meanlog=a_y,sdlog=b_y))
  }
  else if(fun1=="burr" && fun2=="gamma"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pburr(x,shape1=a_x,shape2=b_x,rate=c_x),pgamma(y,shape=a_y,rate=b_y))
  }
  else if(fun1=="burr" && fun2=="weibull"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pburr(x,shape1=a_x,shape2=b_x,rate=c_x),pweibull(y,shape=a_y,scale=b_y))
  }
  else if(fun1=="lnorm" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    d<-cbind(plnorm(x,meanlog=a_x,sdlog=b_x),pburr(y,shape1=a_y,shape2=b_y,rate=c_y))
  }
  else if(fun1=="gamma" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    d<-cbind(pgamma(x,shape=a_x,rate=b_x),pburr(y,shape1=a_y,shape2=b_y,rate=c_y))
  }
  else if(fun1=="weibull" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    d<-cbind(pweibull(x,shape=a_x,scale=b_x),pburr(y,shape1=a_y,shape2=b_y,rate=c_y))
  }
  else if(fun1=="lnorm" && fun2=="lnorm"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(plnorm(x,meanlog=a_x,sdlog=b_x),plnorm(y,meanlog=a_y,sdlog=b_y))
  }
  else if(fun1=="lnorm" && fun2=="gamma"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(plnorm(x,meanlog=a_x,sdlog=b_x),pgamma(y,shape=a_y,rate=b_y))
  }
  else if(fun1=="lnorm" && fun2=="weibull"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(plnorm(x,meanlog=a_x,sdlog=b_x),pweibull(y,shape=a_y,scale=b_y))
  }
  else if(fun1=="gamma" && fun2=="lnorm"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pgamma(x,shape=a_x,rate=b_x),plnorm(y,meanlog=a_y,sdlog=b_y))
  }
  else if(fun1=="weibull" && fun2=="lnorm"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pweibull(x,shape=a_x,scale=b_x),plnorm(y,meanlog=a_y,sdlog=b_y))
  }
  else if(fun1=="gamma" && fun2=="gamma"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pgamma(x,shape=a_x,rate=b_x),pgamma(y,shape=a_y,rate=b_y))
  }
  else if(fun1=="gamma" && fun2=="weibull"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pgamma(x,shape=a_x,rate=b_x),pweibull(y,shape=a_y,scale=b_y))
  }
  else if(fun1=="weibull" && fun2=="gamma"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pweibull(x,shape=a_x,scale=b_x),pgamma(y,shape=a_y,rate=b_y))
  }
  else{
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    d<-cbind(pweibull(x,shape=a_x,scale=b_x),pweibull(y,shape=a_y,scale=b_y))
  }
  if(c[1]==1){
    est<-fitCopula(normalCopula(),d,method="ml")
    return(est)
  }
  else if(c[1]==2){
    est<-fitCopula(tCopula(),d,method="ml")
    return(est)
  }
  else if(c[1]==3){
    est<-fitCopula(claytonCopula(),d,method="ml")
    return(est)
  }
  else if(c[1]==4){
    est<-fitCopula(gumbelCopula(),d,method="ml")
    return(est)
  }
  else if(c[1]==5){
    est<-fitCopula(frankCopula(),d,method="ml")
    return(est)
  }
  else{
    est<-fitCopula(joeCopula(),d,method="ml")
    return(est)
  }
}

int_ifme<-function(data,fun1,fun2){
  a<-ifme(data,fun1,fun2)
  return(confint(a))
}

mm<-function(data){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  c<-copula(data)
  tau<-cor(x,y,method="kendall")
  rho<-cor(x,y,method="spearman")
  if(c[1]==1){
    itau<-fitCopula(normalCopula(),pobs(mat),method="itau")
    irho<-fitCopula(normalCopula(),pobs(mat),method="irho")
  }
  else if(c[1]==2){
    itau<-fitCopula(tCopula(),pobs(mat),method="itau")
    irho<-fitCopula(tCopula(),pobs(mat),method="irho")
  }
  else if(c[1]==3){
    itau<-fitCopula(claytonCopula(),pobs(mat),method="itau")
    irho<-fitCopula(claytonCopula(),pobs(mat),method="irho") 
  }
  else if(c[1]==4){
    itau<-fitCopula(gumbelCopula(),pobs(mat),method="itau")
    irho<-fitCopula(gumbelCopula(),pobs(mat),method="irho")
  }
  else if(c[1]==5){
    itau<-fitCopula(frankCopula(),pobs(mat),method="itau")
    irho<-fitCopula(frankCopula(),pobs(mat),method="irho")
  }
  else{
    itau<-fitCopula(joeCopula(),pobs(mat),method="itau")
    irho<-fitCopula(joeCopula(),pobs(mat),method="irho")
  }
return(c(itau,irho))
}

int_mm<-function(data){
  a<-mm(data)
  return(c(confint(a[[1]]),confint(a[[2]])))
}

empirical<-function(data,titulo){
  n<-dim(data)[[1]]
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  cop<-fitcop(data)
  r<-rCopula(n,cop)
  v<-pobs(mat)
  ec<-C.n(v,r)
  true<-pCopula(v,cop)
  mce<-round(mean(abs(true-ec)/true)*100,2)
  return(mce)
}

# ==================
# Copula: Bivariate
# ==================

biv<-function(data,fun1,fun2){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  cop<-fitcop(data)
  if(fun1=="burr" && fun2=="burr"){
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    c_x<-fit(data,fun1)[3]
    a_y<-fit(data,fun2)[5]
    b_y<-fit(data,fun2)[6]
    c_y<-fit(data,fun2)[7]
    dist<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y,c_y)))
  }
  else if(fun1=="burr" || fun2=="burr"){
    if(fun1=="burr"){
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      c_x<-fit(data,fun1)[3]
      a_y<-fit(data,fun2)[4]
      b_y<-fit(data,fun2)[5]
      dist<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y)))
    }
    else{
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      a_y<-fit(data,fun2)[5]
      b_y<-fit(data,fun2)[6]
      c_y<-fit(data,fun2)[7]
      dist<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y,c_y)))
    }
  }
  else{
    a_x<-fit(data,fun1)[1]
    b_x<-fit(data,fun1)[2]
    a_y<-fit(data,fun2)[4]
    b_y<-fit(data,fun2)[5]
    dist<-mvdc(cop,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y)))}
}

# ===================================
# Copula: Random Sample w/ Bivariate
# ===================================

rsb<-function(data,fun1,fun2,titulo1,titulo2){
  n<-dim(data)[[1]]
  x<-data$obs
  y<-data$sim
  dist<-biv(data,fun1,fun2)
  am<-rMvdc(n,dist)
  pdf<-dMvdc(am,dist)
  par(mfrow=c(1,3))
  plot(x,y,xlab="Observed Wind",ylab="Simulated Wind",main=titulo1,col='royalblue',pch=0,cex=0.5,xlim=c(0,max(am,x,y)),ylim=c(0,max(am,x,y)))
  points(am[,1],am[,2],col="red3",pch=8,cex=0.5)
  legend('topleft',c('Original data','Simulated data\nby Copula'),col=c('royalblue','red3'),pch=c(0,8),cex=0.8)
  scatterplot3d(am[,1],am[,2], pdf, color="red3", main="Density", xlab = "Observed Wind", ylab="Simulated Wind", zlab="PDF",pch=8)
  #par(mfrow=c(1,1))
  contour(dist, dMvdc, xlim = c(0, max(x)), ylim=c(0, max(y)), xlab="Observed Wind", ylab="Simulated Wind", main = titulo2)
}

prob<-function(data,fun1,fun2,p){
  x<-data$obs
  y<-data$sim
  mat<-cbind(x,y)
  c<-copula(data)
  if(c[1]==2){
    rho<-c[2]
    eta<-c[4]
    t<-tCopula(param=rho,df=round(eta,0),df.fixed=T)
    if(fun1=="burr" && fun2=="burr"){
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      c_x<-fit(data,fun1)[3]
      a_y<-fit(data,fun2)[5]
      b_y<-fit(data,fun2)[6]
      c_y<-fit(data,fun2)[7]
      dist<-mvdc(t,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y,c_y)))
    }
    else if(fun1=="burr" || fun2=="burr"){
      if(fun1=="burr"){
        a_x<-fit(data,fun1)[1]
        b_x<-fit(data,fun1)[2]
        c_x<-fit(data,fun1)[3]
        a_y<-fit(data,fun2)[4]
        b_y<-fit(data,fun2)[5]
        dist<-mvdc(t,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x,c_x),list(a_y,b_y)))
      }
      else{
        a_x<-fit(data,fun1)[1]
        b_x<-fit(data,fun1)[2]
        a_y<-fit(data,fun2)[5]
        b_y<-fit(data,fun2)[6]
        c_y<-fit(data,fun2)[7]
        dist<-mvdc(t,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y,c_y)))
      }
    }
    else{
      a_x<-fit(data,fun1)[1]
      b_x<-fit(data,fun1)[2]
      a_y<-fit(data,fun2)[4]
      b_y<-fit(data,fun2)[5]
      dist<-mvdc(t,margins=c(fun1,fun2),paramMargins=list(list(a_x,b_x),list(a_y,b_y)))}
  }
  else{
    dist<-biv(data,fun1,fun2)
  }
  H<-pMvdc(c(p,p),dist)
  HC<-1-H
  return(c(H,HC))
}

# ================================
# Copula: Random Sample w/ Copula
# ================================

rsc<-function(data,titulo1,titulo2){
  n<-dim(data)[[1]]
  cop<-fitcop(data)
  copr<-rCopula(n,cop)
  pdf<-dCopula(copr,cop)
  par(mfrow=c(1,2))
  plot(copr[,1], copr[,2], col="red3", main=titulo1, xlab = "u", ylab = "v",pch=8,cex=0.5)
  scatterplot3d(copr[,1], copr[,2], pdf, color="red3", main="Density", xlab ="u1", ylab="u2", zlab="PDF", pch="*")
  par(mfrow=c(1,1))
  contour(cop, dCopula, xlim = c(0, 1), ylim=c(0, 1), xlab="u1", ylab="u2", main = titulo2)
}

# ============
# Correlation
# ============

corr<-function(data){
  x<-data$obs
  y<-data$sim
  pearson<-cor(x,y,method="pearson")
  kendall<-cor(x,y,method="kendall")
  spearman<-cor(x,y,method="spearman")
  return(c(round(pearson,4),round(kendall,4),round(spearman,4)))
}













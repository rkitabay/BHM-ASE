analytic_Tdist<-function(n,pi1,pi2,Aq01.star){
  #output the quantiles of Tq0 and Tq1 for each q0.star and q1.star in Step 1 of Section 2.1
  #n: the number of patients for each cancer type
  #pi1: the null response rate
  #pi1: the alternative response rate
  #Aq01.star: q0.star and q1.star common array 
  
  pi1pi2<-dbinom(x=0:n,size=n,prob=pi1)%*%t(dbinom(x=0:n,size=n,prob=pi2))
  pi1pi1<-dbinom(x=0:n,size=n,prob=pi1)%*%t(dbinom(x=0:n,size=n,prob=pi1))
  pi2pi2<-dbinom(x=0:n,size=n,prob=pi2)%*%t(dbinom(x=0:n,size=n,prob=pi2))
  Ttable<-matrix(nrow=n+1,ncol=n+1)
  for(x1 in 0:n){
    for(x2 in 0:n){
      table<-matrix(c(x1,n-x1,x2,n-x2),ncol=2,byrow=T)
      Ttable[x1+1,x2+1]<-scaled_X2(table)
      if(x1==x2&&(x1==0||x1==n))Ttable[x1+1,x2+1]<-0
    }
  }
  Tp<-data.frame(Tkk=c(Ttable),pi1pi2=c(pi1pi2),pi1pi1=c(pi1pi1),pi2pi2=c(pi2pi2))
  sort.Tp<-arrange(Tp,Tkk)
  cums.Tp12<-data.frame(Tkk=sort.Tp$Tkk,pi1pi2=sort.Tp$pi1pi2,cums12=cumsum(sort.Tp$pi1pi2))
  ATq1<-c()
  for(i in 1:length(Aq01.star)){
    up.Tp12<-cums.Tp12[cums.Tp12$cums12>Aq01.star[i],]
    ATq1[i]<-up.Tp12$Tkk[1]
  }
  ATq1[ATq1<0.01]<-0.01
  
  Tp1122<-data.frame(Tkk=rep(sort.Tp$Tkk,2),pp=c(sort.Tp$pi1pi1/2,sort.Tp$pi2pi2/2))
  sort.Tp1122<-arrange(Tp1122,Tkk)
  cums.Tp1122<-data.frame(Tkk=sort.Tp1122$Tkk,pp=sort.Tp1122$pp,cums=cumsum(sort.Tp1122$pp))
  ATq0<-c()
  for(i in 1:length(Aq01.star)){
    up.Tp1122<-cums.Tp1122[cums.Tp1122$cums>Aq01.star[i],]
    ATq0[i]<-up.Tp1122$Tkk[1]
  }
  ATq0[ATq0<0.01]<-0.01
  
  list(Tq1=ATq1,Tq0=ATq0)
}

scaled_X2<-function(m){
  #output the scaled Chi-squared test statistic T in the equation (1)
  #m: 2 x 2 division table
  xa<-m[1,1]
  ya<-m[1,2]
  na<-xa+ya
  xb<-m[2,1]
  yb<-m[2,2]
  nb<-xb+yb
  ea0<-na*(ya+yb)/(na+nb)
  eb0<-nb*(ya+yb)/(na+nb)
  ea1<-na*(xa+xb)/(na+nb)
  eb1<-nb*(xa+xb)/(na+nb)
  t<-((ya-ea0)^2/ea0+(yb-eb0)^2/eb0+(xa-ea1)^2/ea1+(xb-eb1)^2/eb1)/(max(c(na,nb)))^(1/4)
  t
}


omega_m<-function(data,Nbask,pi.null,pi.alt){
  #output omega_m and m in the equation (7)
  #The arguments are the same as those defined respectively in implementation.R
  
  model.all<-dec2bin(1:2^Nbask-1,digit=Nbask)
  Neff<-rowSums(model.all)
  
  pri.scen<-c()
  for(i in 0:Nbask){
    pri.scen[i+1]<-1/(Nbask+1)/sum(Neff==i)
  }
  Alike<-c()
  pri.all<-c()
  for(i in 1:nrow(model.all)){
    Api<-rep(pi.null,Nbask)
    Api[model.all[i,]==1]<-pi.alt
    llike<-log(choose(data$n,data$x))+data$x*log(Api)+(data$n-data$x)*log(1-Api)
    Alike[i]<-exp(sum(llike))
    pri.all[i]<-pri.scen[Neff[i]+1]
  }
  sum.model.p<-sum(Alike*pri.all)
  Ap<-c()
  for(i in 1:nrow(model.all)){
    Ap[i]<-Alike[i]*pri.all[i]/sum.model.p
  }
  pes.tmp<-c()
  for(i in 0:Nbask){
    pes.tmp[i+1]<-sum(Ap[Neff==i])
  }
  data.frame(Neff=0:Nbask,p=pes.tmp)
}

dec2bin <- function(num, digit=10){
  #Convert a binary number to a decimal number
  #num: binary number
  
  bin.tab<-matrix(0,nrow=length(num),ncol=digit)
  for(j in 1:length(num)){
    bin<-c()
    tmp<-num[j]
    for(i in digit:1){
      bin[i]<-tmp%/%2^(i-1)
      tmp<-tmp-bin[i]*2^(i-1)
    }
    bin.tab[j,]<-bin
  }
  bin.tab
}

determine_Tq<-function(gs.result,wm,Nbask,pi.null,pi.alt,loss.par){
  #output Tq0 and Tq1 in Step 5 of Section 2.1
  #The arguments are the same as those defined respectively in implementation.R
  
  alpha<-loss.par$ga
  beta<-loss.par$gb
  
  m<-1
  for(i in 1:nrow(gs.result)){
    
    if(gs.result$Neff[i]==1)loss.retain<-0
    loss.retain<-loss.retain+loss_function(res=gs.result[i,-(1:3)],Nbask,wm,
                                           alpha=alpha,beta=beta,Neff=gs.result$Neff[i],
                                           loss.par$wa,loss.par$wb)
    
    if(gs.result$Neff[i]==Nbask){
      res<-data.frame(q0=gs.result$q0[i],q1=gs.result$q1[i],lf=loss.retain)
      if(m==1) gs.lf<-res
      else gs.lf<-rbind(gs.lf,res)
      m<-m+1 
    }
  }
  gs.lf[which.min(gs.lf$lf),c("q0","q1")]
}


loss_function<-function(res,Nbask,wm,alpha,beta,Neff,wa,wb){
  #output the values of the loss function in the equation (6)
  
  l<-0
  for(j in 1:Nbask){
    if(j<=(Nbask-Neff)){
      l<-l+wm[wm$Neff==Neff,"p"]*wa*max(res[j]-alpha,0)
    }else{
      l<-l-wm[wm$Neff==Neff,"p"]*wb*min(res[j]-beta,0)
    }
  }
  l
}


elastic_function<-function(t,Tq0,Tq1,c0,c1){
  #convert the statistic T to the similarity
  #t: the statistic T in the equation (1)
  #The arguments except for t are the same as those defined respectively in implementation.R
  
  b <- log((1-c0)*c1/((1-c1)*c0))/((log(Tq0))-(log(Tq1)))
  a <- log((1-c0)/c0)-b*(log(Tq0))
  
  l<-1/(1+exp(a+b*log(t)))
  l
}

#install.packages(tolerance)
library(tolerance)
arl=function(l,n,L,r,d1,d2){
  rl=function(l,n,L){
    w=0.1
    mu0=0;sigma0=1
    x=matrix(r2exp(250000,d1,d2),ncol=5);
    mu_hat=apply(x,1,min);sigma_hat=apply(x,1,mean)-mu_hat
    y1=2*n*(mu_hat-mu0)/sigma0;y2=2*n*sigma_hat/sigma0
    z1=qnorm(pchisq(y1,2));z2=qnorm(pchisq(y2,2*(n-1)))
    u1=0;v1=0;W=0;u0=0;v0=0;vw=0;UCL=0;LCL=0
    for(i in 1:length(z1)){
      u1[i]=(1-l)*u0+l*z1[i]
      u0=u1[i]
      v1[i]=(1-l)*v0+l*z2[i]
      v0=v1[i]
      W[i]=w*u1[i]+(1-w)*v1[i]
      vw[i]=(1-2*w+2*w^2)*l*(1-(1-l)^(2*i))/(2-l)
      UCL[i]=L*sqrt(vw[i])
      LCL[i]=-L*sqrt(vw[i])
    }
    
    count=0
    for(i in 1:length(W)){
      if(W[i]>LCL[i] & W[i]<UCL[i])
      {
        count=count+1
      }else{
        print(count)
        break
      }
    }
    RL=count+1
    RL
  }
  ana_new1=replicate(r,rl(l,n,L))
  aa=mean(ana_new1);aa1=sd(ana_new1);aa2=aa1/sqrt(r)
  metric=c(aa,aa1,aa2);metric
}

lim_k2=function(arl0,l,n,w,r){
  L1=2.7;L2=2.8
  arl1=arl(l,n,L1,w,r)-arl0;arl2=arl(l,n,L2,w,r)-arl0
  print(c(arl1,arl2))
  
  while(abs(L1-L2)>0.01){
    L=(L1+L2)/2;arl3=arl(l,n,L,w,r)-arl0
    if(arl3<0) L1=L else L2=L
  }
  L
}


##to make table for different m and delta for AMRL and SDMRL
d2=c(0.1,0.25,0.5,1.0);d1=c(0.4,0.6,0.8,0.9,1,1.1,1.2,1.5,2);r=100000;n=5;l=0.1;L=2.82
#col_names=c("ARL","SDRL","SERL","ARL","SDRL","SERL","ARL","SDRL","SERL","ARL","SDRL","SERL","ARL","SDRL","SERL","ARL","SDRL","SERL");row_names=c("0.4","0.6","0.8","0.9","1","1.1","1.2","1.5","2")
ARL=matrix(0,length(d1),length(d2));SDRL=matrix(0,length(d1),length(d2));SERL=matrix(0,length(d1),length(d2))

for(j in 1:length(d1))
{
  for (i in 1:length(d2))
  {
    metric=arl(l,n,L,r,d1[j],d2[i])
    ARL[j,i]=metric[1]
    SDRL[j,i]=metric[2]
    SERL[j,i]=metric[3]
  }
  
}

write.table(cbind(ARL,SDRL,SERL),"tab_0.1l.csv",sep=",")

#bined=matrix(0,length(d1),18)
#bined[,seq(1,18,3)]=ARL
#bined[,seq(2,18,3)]=SDRL
#bined[,seq(3,18,3)]=SERL
#write.table(bined,"bined.csv",sep=",")

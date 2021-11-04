library(tidyverse)
library(GGally)
library(ggplot2)
library(rayshader)
library(sfsmisc)
library(xtable)
library(dr)
####################
n.b=13
r=1.3
B.star=3
B=B.star
for (i in 1:(n.b-1)) {
  B.star=B.star/r
  B=c(B,B.star)
}

B1=B
B2=B

Fold=5
#####exponetial kernels#######
Ker_exp=function(x){
  p=length(x)
  return(1/sqrt(2*pi)^p*exp(-1/2*t(x)%*%(x)))
}
#############Kernels###################


KR=function(E,Y,X,Ker=Ker_exp,B)
{ 
  
  
  Y=as.matrix(Y)
  X=as.matrix(X)
  p=dim(X)[2]
  n=dim(X)[1]
  E=as.matrix(E)
  nE=dim(E)[1]
  hat_Y=NULL
  
  for (i in 1:nE) {
    temp_N=0
    temp_D=0
    for (j in 1:n) {
      temp_N=temp_N+Ker((E[i,]-X[j,])/B)*Y[j]
      temp_D=temp_D+Ker((E[i,]-X[j,])/B)
      
    }
    hat_Y=c(hat_Y,temp_N/temp_D)
    
    
    
  }
  
  
  return(hat_Y)
}



CKR=function(E,Y,X,gh,Constraints=TRUE,C,Y_KR,Ker=Ker_exp,B_2)
{
  Y=as.matrix(Y)
  X=as.matrix(X)
  E=as.matrix(E)
  n=dim(X)[1]
  C=as.matrix(C)
  size=dim(E)[1]
  if(Constraints==FALSE) {
    P=diag(n)-diag(n)
  }else{
    P=C%*%solve(t(C)%*%C)%*%t(C)
  }
  
  
  
  Y_hat=(diag(n)-P)%*%Y_KR+P%*%gh
  
  return( KR(E,Y_hat,X,Ker ,B_2))
}


###################
D.E=read_csv("data3.csv")
D.E=D.E%>%filter(sysbp_min>50)
D.E=D.E%>%dplyr::select(sysbp_mean,hr_mean,resprate_mean,spo2_mean)
D.E=D.E%>%drop_na()
#D.E=D.E%>%mutate_at(2:4,scale)
#####################

direct=dr(sysbp_mean~.,data=D.E,method="save")
dd=direct$evectors[,1]
D.E=D.E%>%mutate( V.X=dd[1]*hr_mean+dd[2]*resprate_mean+dd[3]*spo2_mean )

var(dr.direction(direct,x=D.E[,-c(1,4)])[,1])
var(dr.direction(direct,x=D.E[,-c(1,4)])[,2])
var(dr.direction(direct,x=D.E[,-c(1,4)])[,3])

aa=apply(dr.direction(direct,x=D.E[,-c(1,4)]),2,var)
aa/sum(aa)

fit=lm(sysbp_mean~V.X,data=D.E)



#####################

#################


V.Z=c("etO2","inO2")
V.Y="NBP (Sys)"
B.record=NULL

  D.I=read_csv("Data_32PP.csv")
  D.I=D.I%>%rename(inCO2=imCO2)
  
  D.I2=D.I%>%dplyr::select(V.Y,V.Z,HR,SpO2,awRR)%>%drop_na()%>%mutate( V.X=dd[1]*HR+dd[2]*awRR+dd[3]*SpO2      )
  
  ###########
  DD=D.I2%>%dplyr::select( V.Y,V.Z,V.X)
  
  Z=DD%>%dplyr::select(V.Z)%>%scale()
  Y=DD%>%dplyr::select(V.Y)
  X=as.matrix(DD%>%dplyr::select(V.X)%>%scale)
  D=tibble(Y,X,Z)
  
  D.internal=D[complete.cases(D),]
  n=dim(D.internal)[1]
  #########CV############
  
  
  
  
  
  Out=NULL
  RMSE.KR=array(0,dim =c(length(B),Fold)  )
  RMSE.DKR=array(0,dim = c(length(B1),length(B2),Fold))
  RMSE.CKR=array(0,dim = c(length(B1),length(B2),Fold))
  # RMSE.CKR2=array(0,dim = c(length(B1),length(B2),length(B3),Time))
  
  
  SS=sample(1:n)
  
  for (t in 1:Fold) {
    ####seperate test and train data
    #Data=cbind(Y=y,x[,c(sample(1:100,1),  sample(101:200,1),sample(201:300,1))] )
    
    t1=Sys.time()
    
    
    S=SS[((t-1)*floor(n/Fold)+1):(t*floor(n/Fold))]
    Data.test=as.matrix(D.internal[S,])
    Data.train=as.matrix(D.internal[-S,])
    n.train=dim(Data.train)[1]
    ###KR###
    
    
    for (b in 1:length(B)) {
      
      
      P.KR=KR(E=Data.test[,-1],  X=Data.train[,-1]  ,Y=Data.train[,1],B=B[b])
      
      e.KR=(P.KR-Data.test[,1])^2
      
      RMSE.KR[b,t]=sqrt(mean(e.KR,na.rm = TRUE))
      
    }
    
    #CKR-Level1##
    
    for (b1 in 1:length(B1)) {
      for (b2 in 1:length(B2)) {
        
        
        Y_KR=KR(E=Data.train[,-1],  
                X=Data.train[,-1],Y=Data.train[,1],B=B1[b1])
        
        a=fit$coefficients
        mean.X=attr(X,"scaled:center")
        SD.X=attr(X,"scaled:scale")
        a1=a[1]+a[2]*mean.X
        a2=a[2]*SD.X
        
        
        
        
        hhat.linear=cbind(rep(1,n.train),Data.train[,2])%*%c(a1,a2)
        Constraints0=cbind(rep(1, n.train),Data.train[,2]) 
        
        P.CKR=CKR(E=Data.test[,-1],Y=Data.train[,1],
                  X=Data.train[,-1],gh=hhat.linear,Constraints=TRUE,C=Constraints0,Y_KR=Y_KR,Ker=Ker_exp,B_2=B2[b2])
        
        
        
        e.CKR=(P.CKR-Data.test[,1])^2
        
        RMSE.CKR[b1,b2,t]=sqrt(mean(e.CKR,na.rm = TRUE))
        
        
        P.DKR=CKR(E=Data.test[,-1],Y=Data.train[,1],
                  X=Data.train[,-1],gh=hhat.linear,Constraints=FALSE,C=Constraints0,Y_KR=Y_KR,Ker=Ker_exp,B_2=B2[b2])
        
        
        e.DKR=(P.DKR-Data.test[,1])^2
        
        RMSE.DKR[b1,b2,t]=sqrt(mean(e.DKR,na.rm = TRUE))
        
        
        
        
      }  
    }
    
    
    
    
    t4=Sys.time()
    
    print(paste( "The time period of"   ,"t=", t,  " simulation is ",t4-t1))
    
    
  }
  
  
  RMSE=c(min(apply(RMSE.KR, 1, mean)),
         min(apply(RMSE.DKR, c(1,2), mean)),
         min(apply(RMSE.CKR, c(1,2), mean))
  )
  
  
  
  
  index.KR=which.min(apply(RMSE.KR, 1, mean))
  index.DKR=as.vector(which(apply(RMSE.DKR, c(1,2), mean)==min(apply(RMSE.DKR, c(1,2), mean)),arr.ind = TRUE))
  index.CKR=as.vector( which(apply(RMSE.CKR, c(1,2), mean)==min(apply(RMSE.CKR, c(1,2), mean)),arr.ind = TRUE))
  
  
  ##################
  Data.train=as.matrix(D.internal)
  
  n=dim(Data.train)[1]
  hhat.linear=cbind(rep(1, n ),Data.train[,2])%*%fit$coefficients
  Constraints0=cbind(rep(1, n),Data.train[,2]) 
  #####################
  
  B.KR=B[index.KR]
  B.CK=c(B1[index.CKR[1]],B2[index.CKR[2]])
  B.DK=c(B1[index.DKR[1]],B2[index.DKR[2]])
  B.temp=NULL
  B.temp=cbind( Z=c(NA,NA),B.KR=c(B.KR,NA), B.CK=B.CK,B.DK=B.DK )
  
  B.record=rbind(B.record,B.temp)
  
  
  ##################
  
  
  Gx=range(Data.train[,2])[1]+((0:100)/100)*(range(Data.train[,2])[2]-range(Data.train[,2])[1])
  Gz1=range(Data.train[,3])[1]+((0:100)/100)*(range(Data.train[,3])[2]-range(Data.train[,3])[1])
  Gz2=range(Data.train[,4])[1]+((0:100)/100)*(range(Data.train[,4])[2]-range(Data.train[,4])[1])

   E=expand.grid( Gx  ,Gz1,Gz2)
  
  
  P.KR=KR(E=E,  X=Data.train[,-1]  ,Y=Data.train[,1],B=B.KR)
  
  Y_KR.CK=KR(E=Data.train[,-1],  
             X=Data.train[,-1],Y=Data.train[,1],B=B.CK[1])
  
  P.CKR=CKR(E=E,Y=Data.train[,1],
            X=Data.train[,-1],gh=hhat.linear,Constraints=TRUE,C=Constraints0,Y_KR=Y_KR.CK,Ker=Ker_exp,B_2=B.CK[2])
  
  Y_KR=KR(E=Data.train[,-1],  
          X=Data.train[,-1],Y=Data.train[,1],B=B.DK[1])
  
  P.DKR=CKR(E=E,Y=Data.train[,1],
            X=Data.train[,-1],gh=hhat.linear,Constraints=FALSE,C=Constraints0,Y_KR=Y_KR,Ker=Ker_exp,B_2=B.DK[2])
  ##############
  
  
  X=as.matrix(DD%>%dplyr::select(V.X))
  
  Z=as.matrix(DD%>%dplyr::select(V.Z))
  
  Out.1=tibble(CK=P.CKR, KR=P.KR, DK=P.DKR,X=sd(X)*(E[,1])+mean(X),Z1=sd(Z[,1],na.rm =TRUE)*E[,2]+mean(Z[,1],na.rm =TRUE)
               ,Z2=sd(Z[,2],na.rm =TRUE)*E[,3]+mean(Z[,2],na.rm =TRUE))
  
  Out.2=Out.1%>%dplyr::select(-Z2)%>%group_by(X,Z1)%>%summarise(CK= mean(CK),KR=mean(KR),DK=mean(DK))
  Out.2=Out.2%>%gather(key="Method",value="value",-X,-Z1)
  
  
  P=Out.2%>%ggplot(aes(x=X,y=Z1))+
    geom_raster(aes(fill=value))+
    facet_wrap(~Method)+ylab(V.Z[1])+
    labs(fill="Sys" )+ggtitle(" ")+
    scale_fill_gradientn(  colours=c("white",  "black"))
  P
 ggsave(paste0(V.Y,"&",V.Z[1],"2.png"))
  

 Out.2=Out.1%>%dplyr::select(-Z1)%>%group_by(X,Z2)%>%summarise(CK= mean(CK),KR=mean(KR),DK=mean(DK))
 Out.2=Out.2%>%gather(key="Method",value="value",-X,-Z2)
 
 
 
 P=Out.2%>%ggplot(aes(x=X,y=Z2))+
   geom_raster(aes(fill=value))+
   facet_wrap(~Method)+ylab(V.Z[2])+
   labs(fill="Sys" )+ggtitle(" ")+
   scale_fill_gradientn(  colours=c("white",  "black"))
 P
 ggsave(paste0(V.Y,"&",V.Z[2],"2.png"))
 
  

B.record
########################

rm(list=ls())
library('akima') #for linear interpolation
library('data.table')

#get covariates
setwd('U:\\GIT_models\\git_ssf\\fake data')
source('aux functions.R')

tmp=fread('fake data xmat.csv')
xdim=1000
ydim=1000
ncov=3
xmat=array(unlist(tmp),c(xdim,ydim,ncov))

#get data
coord=data.matrix(read.csv('fake data coord.csv',as.is=T))

#parameters
alpha=c(1,-1,0,0) #parameters for time model

#calculate mean time
mean.time=matrix(NA,xdim,ydim)
for (i in 1:xdim){
  for (j in 1:ydim){
    mean.time[i,j]=exp(c(1,xmat[i,j,])%*%alpha)
  }
}

#important variables
nsim=nrow(coord)
window1=20
ndados=((window1*2)+1)^2
store.calc=matrix(NA,(nsim-1)*ndados,5+ncov)
oooo=1
for (i in 1:(nsim-1)){
  print(i)
  
  #calculate cumulative time across 8 directions
  res=get.cumtime.8dir(coord1=coord[i,],mean.time=mean.time,window1=window1)
  
  #calculate prob of staying in the sample place
  tmp=matrix(coord[i,c('y','x')],1,2)
  tempo.not.move=mean.time[tmp]
  res1=rbind(res,c(tempo.not.move,coord[i,c('x','y')]))
  
  #interpolate cumulative time across a window1 x window1 area
  minx=coord[i,'x']-window1; maxx=coord[i,'x']+window1
  miny=coord[i,'y']-window1; maxy=coord[i,'y']+window1
  res2 = interp(x=res1[,'x'],y=res1[,'y'],z=res1[,'cum.time'],
                xo = minx:maxx, yo = miny:maxy, 
                linear=T, extrap = F)
  res3=interp2xyz(res2)
  colnames(res3)[3]='cum.time'
  
  #get covariates
  aux=numeric()
  for (j in 1:ncov){
    tmp=cbind(res3[,'y'],res3[,'x'],j)
    tmp1=xmat[tmp]
    aux=cbind(aux,tmp1)
  }
  colnames(aux)=paste0('cov',1:ncov)
  res4=cbind(res3,aux)

  #store results  
  ind=which(res4[,'x']==coord[i,'x'] & res4[,'y']==coord[i,'y'])
  zeros=rep(0,nrow(res4))
  zeros[ind]=1
  
  seq1=oooo:(oooo+ndados-1)
  store.calc[seq1,]=cbind(res4,zeros,i)  
  oooo=oooo+ndados
}
colnames(store.calc)=c(colnames(res4),'selected','mov.id')

setwd('U:\\GIT_models\\git_ssf\\fake data')
write.csv(store.calc,'cumul_time and covs.csv',row.names=F)

#for this code to work, we need the x and y coordinates to function as indexes which can be used to subset the data


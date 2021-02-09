
library(akima) #for linear interpolation
library(here)

set.seed(1)

source(here('fake data', 'aux functions.R'))


#basic settings
ncov=3
xdim=1000
ydim=1000
nloc=xdim*ydim
xmat=array(runif(nloc*ncov,min=-1,max=1),dim=c(xdim,ydim,ncov))
ind.loc=matrix(1:nloc,xdim,ydim)  #indexes from top to bottom, left to right

#parameters
alpha=c(1,-1,0,0) #parameters for time model
betas=c(1,-1,0) #parameters for resource selection function
gamma.b=2

#calculate mean time and preference
mean.time=mean.pref=matrix(NA,xdim,ydim)
for (i in 1:xdim){
  for (j in 1:ydim){
    mean.time[i,j]=exp(c(1,xmat[i,j,])%*%alpha)
    mean.pref[i,j]=exp(xmat[i,j,]%*%betas)
  }
}

#simulate movement
current.pos=ind.loc[500,500]
time.int=4
nsim=1000
coord=matrix(NA,nsim,2)
coord[1,]=c(500,500)
colnames(coord)=c('x','y')
window1=20

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

  #get preference
  tmp=cbind(res3[,'y'],res3[,'x'])
  pref=mean.pref[tmp]
  res4=cbind(res3,pref)

  #calculate probabilities
  # tempo.rel=res4[,'cum.time']/max(res4[,'cum.time'])
  # plot(res4[,'x'],res4[,'y'],col=grey(tempo.rel))
  
  prob=dgamma(time.int,gamma.b*res4[,'cum.time'],gamma.b)*res4[,'pref']
  soma=sum(prob)
  prob=prob/soma
  #plot(prob,type='h')
  # prob1=prob-min(prob)
  # prob1=prob1/max(prob1)
  #plot(res4[,'x'],res4[,'y'],col=grey(1-prob1))
  
  #decide where to move
  tmp=rmultinom(1,size=1,prob=prob)
  ind=which(tmp==1)
  res5=res4[ind,]
  
  coord[i+1,c('x','y')]=res5[c('x','y')]
}

#show trajectory
plot(coord[,'x'],coord[,'y'],type='l')

#save results
write.csv(coord, here('fake data', 'fake data coord.csv'), row.names=F)
xmat1=matrix(xmat,1,xdim*ydim*ncov)
write.csv(xmat1, here('fake data', 'fake data xmat.csv'), row.names=F)

rm(list=ls())
set.seed(1)

#basic settings
ncov=3
xdim=1000
ydim=1000
nloc=xdim*ydim
xmat=array(runif(nloc*ncov,min=-1,max=1),dim=c(xdim,ydim,ncov))
ind.loc=matrix(1:nloc,xdim,ydim)

#parameters
alpha=c(1,-1,0,0) #parameters for time model
betas=c(-1,1,2) #parameters for resource selection function
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
time.int=2
nsim=1000
coord=matrix(NA,nsim,2)
coord[1,]=c(500,500)
colnames(coord)=c('x','y')

for (i in 1:(nsim-1)){
  print(i)
  
  #calculate prob 4 directions
  #up
  final.coord=coord[i,'y']+20
  seq1=(coord[i,'y']+1):final.coord
  tempo.up=cumsum(mean.time[seq1,coord[i,'x']])
  pref.up=mean.pref[seq1,coord[i,'x']]  
  tmp.up=cbind(tempo.up,pref.up,coord[i,'x'],seq1)
  
  #down
  final.coord=coord[i,'y']-20
  seq1=(coord[i,'y']-1):final.coord
  tempo.do=cumsum(mean.time[seq1,coord[i,'x']])
  pref.do=mean.pref[seq1,coord[i,'x']]  
  tmp.do=cbind(tempo.do,pref.do,coord[i,'x'],seq1)
  
  #right
  final.coord=coord[i,'x']+20
  seq1=(coord[i,'x']+1):final.coord
  tempo.ri=cumsum(mean.time[seq1,coord[i,'y']])
  pref.ri=mean.pref[coord[i,'y'],seq1]  
  tmp.ri=cbind(tempo.ri,pref.ri,seq1,coord[i,'y'])
  
  #left
  final.coord=coord[i,'x']-20
  seq1=(coord[i,'x']-1):final.coord
  tempo.le=cumsum(mean.time[seq1,coord[i,'y']])
  pref.le=mean.pref[coord[i,'y'],seq1]  
  tmp.le=cbind(tempo.le,pref.le,seq1,coord[i,'y'])

  #calculate prob 4 diagonals

  #upper right
  final.coord.y=coord[i,'y']+20
  final.coord.x=coord[i,'x']+20
  seq1.y=(coord[i,'y']+1):final.coord.y
  seq1.x=(coord[i,'x']+1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  tempo.ur=cumsum(mean.time[tmp])
  pref.ur=mean.pref[tmp]  
  tmp.ur=cbind(tempo.ur,pref.ur,seq1.x,seq1.y)
  
  #upper left
  final.coord.y=coord[i,'y']+20
  final.coord.x=coord[i,'x']-20
  seq1.y=(coord[i,'y']+1):final.coord.y
  seq1.x=(coord[i,'x']-1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  tempo.ul=cumsum(mean.time[tmp])
  pref.ul=mean.pref[tmp]  
  tmp.ul=cbind(tempo.ul,pref.ul,seq1.x,seq1.y)
  
  #lower right
  final.coord.y=coord[i,'y']-20
  final.coord.x=coord[i,'x']+20
  seq1.y=(coord[i,'y']-1):final.coord.y
  seq1.x=(coord[i,'x']+1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  tempo.lr=cumsum(mean.time[tmp])
  pref.lr=mean.pref[tmp]  
  tmp.lr=cbind(tempo.lr,pref.lr,seq1.x,seq1.y)
  
  #lower left
  final.coord.y=coord[i,'y']-20
  final.coord.x=coord[i,'x']-20
  seq1.y=(coord[i,'y']-1):final.coord.y
  seq1.x=(coord[i,'x']-1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  tempo.ll=cumsum(mean.time[tmp])
  pref.ll=mean.pref[tmp]  
  tmp.ll=cbind(tempo.ll,pref.ll,seq1.x,seq1.y)
  
  #calculate prob of staying in the sample place
  tmp=matrix(coord[i,c('y','x')],1,2)
  tempo.not.move=mean.time[tmp]
  pref.not.move=mean.pref[tmp]
  tmp.not.move=c(tempo.not.move,pref.not.move,coord[i,'x'],coord[i,'y'])
  
  #combine all probabilities
  final=rbind(tmp.up,tmp.do,tmp.ri,tmp.le,
              tmp.ur,tmp.ul,tmp.lr,tmp.ll,
              tmp.not.move)
  colnames(final)=c('cum.time','pref','x','y')

  # tempo.rel=final[,'cum.time']/max(final[,'cum.time'])
  # plot(final[,'x'],final[,'y'],col=grey(tempo.rel))
  
  prob=dgamma(time.int,final[,'cum.time'],gamma.b)*final[,'pref']
  soma=sum(prob)
  prob=prob/soma
  #plot(prob,type='h')
  
  #decide where to move
  tmp=rmultinom(1,size=1,prob=prob)
  ind=which(tmp==1)
  final1=final[ind,]
  
  coord[i+1,c('x','y')]=final1[c('x','y')]
}

#show trajectory
plot(coord[,'x'],coord[,'y'],type='l')


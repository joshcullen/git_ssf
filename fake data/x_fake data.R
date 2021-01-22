rm(list=ls())

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
  soma.tempo=cumsum(mean.time[seq1,coord[i,'x']])
  p.up=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[seq1,coord[i,'x']]  
  plot(p.up,type='h')
  tmp.up=cbind(p.up,coord[i,'x'],seq1)
  
  #down
  final.coord=coord[i,'y']-20
  seq1=(coord[i,'y']-1):final.coord
  soma.tempo=cumsum(mean.time[seq1,coord[i,'x']])
  p.do=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[seq1,coord[i,'x']]  
  plot(p.do,type='h')
  tmp.do=cbind(p.do,coord[i,'x'],seq1)
  
  #right
  final.coord=coord[i,'x']+20
  seq1=(coord[i,'x']+1):final.coord
  soma.tempo=cumsum(mean.time[seq1,coord[i,'y']])
  p.ri=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[coord[i,'y'],seq1]  
  plot(p.ri,type='h')
  tmp.ri=cbind(p.ri,seq1,coord[i,'y'])
  
  #left
  final.coord=coord[i,'x']-20
  seq1=(coord[i,'x']-1):final.coord
  soma.tempo=cumsum(mean.time[seq1,coord[i,'y']])
  p.le=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[coord[i,'y'],seq1]  
  plot(p.le,type='h')
  tmp.le=cbind(p.le,seq1,coord[i,'y'])

  #calculate prob 4 diagonals

  #upper right
  final.coord.y=coord[i,'y']+20
  final.coord.x=coord[i,'x']+20
  seq1.y=(coord[i,'y']+1):final.coord.y
  seq1.x=(coord[i,'x']+1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  soma.tempo=cumsum(mean.time[tmp])
  p.ur=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[tmp]  
  plot(p.ur,type='h')
  tmp.ur=cbind(p.ur,seq1.x,seq1.y)
  
  #upper left
  final.coord.y=coord[i,'y']+20
  final.coord.x=coord[i,'x']-20
  seq1.y=(coord[i,'y']+1):final.coord.y
  seq1.x=(coord[i,'x']-1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  soma.tempo=cumsum(mean.time[tmp])
  p.ul=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[tmp]  
  plot(p.ul,type='h')
  tmp.ul=cbind(p.ul,seq1.x,seq1.y)
  
  #lower right
  final.coord.y=coord[i,'y']-20
  final.coord.x=coord[i,'x']+20
  seq1.y=(coord[i,'y']-1):final.coord.y
  seq1.x=(coord[i,'x']+1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  soma.tempo=cumsum(mean.time[tmp])
  p.lr=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[tmp]  
  plot(p.lr,type='h')
  tmp.lr=cbind(p.lr,seq1.x,seq1.y)
  
  #lower left
  final.coord.y=coord[i,'y']-20
  final.coord.x=coord[i,'x']-20
  seq1.y=(coord[i,'y']-1):final.coord.y
  seq1.x=(coord[i,'x']-1):final.coord.x
  tmp=cbind(seq1.y,seq1.x); colnames(tmp)=c('y','x')
  soma.tempo=cumsum(mean.time[tmp])
  p.ll=dgamma(time.int,soma.tempo,gamma.b)*mean.pref[tmp]  
  plot(p.ll,type='h')
  tmp.ll=cbind(p.ll,seq1.x,seq1.y)
  
  #calculate prob of staying in the sample place
  tmp=matrix(coord[i,c('y','x')],1,2)
  p.not.move=dgamma(time.int,mean.time[tmp],gamma.b)*mean.pref[tmp]
  tmp.not.move=c(p.not.move,coord[i,'x'],coord[i,'y'])
  
  #combine all probabilities
  final=rbind(tmp.up,tmp.do,tmp.ri,tmp.le,
              tmp.ur,tmp.ul,tmp.lr,tmp.ll,
              tmp.not.move)
  colnames(final)=c('prob','x','y')
  soma=sum(final[,'prob'])
  final[,'prob']=final[,'prob']/soma
  
  #decide where to move
  tmp=rmultinom(1,size=1,prob=final[,'prob'])
  ind=which(tmp==1)
  final1=final[ind,]
  
  coord[i+1,c('x','y')]=final1[c('x','y')]
}

#show trajectory
plot(coord[,'x'],coord[,'y'],type='l')


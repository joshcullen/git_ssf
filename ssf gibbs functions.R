acceptMH <- function(p0,p1,x0,x1){   #accept for M, M-H
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(1)
  accept=0
  if (z<a) {x0 = x1; accept=1}
  list(x = x0, accept = accept)
}
#----------------------------------------------------------------------------------------------
get.llk=function(betas,xmat,time.prob,mov.id,nmov.id,nobs,log.time.prob.cond.sel,cond.sel){
  #get rsf
  log.rsf=xmat%*%betas
  rsf=exp(log.rsf)
  
  #get denominator for each mov.id
  mult=rsf*time.prob
  tot=GetTotal(id=mov.id-1, vec=mult, nid=nmov.id,nobs=nobs)
  
  #calculate llk
  log.num1=log.time.prob.cond.sel+log.rsf[cond.sel]
  log.den1=log(tot)
  sum(log.num1-log.den1)
}
#--------------------------------------------------------
sample.betas=function(betas,ncov,jump1,xmat,time.prob,mov.id,nmov.id,nobs,log.time.prob.cond.sel,cond.sel){
  betas.old=betas
  prior.old=dnorm(betas.old,mean=0,sd=10,log=T)
  accept1=rep(0,ncov)
  for (i in 1:ncov){
    betas.new=betas.old
    betas.new[i]=rnorm(1,mean=betas.old[i],sd=jump1[i])
    llk.old=get.llk(betas=betas.old,xmat=xmat,time.prob=time.prob,mov.id=mov.id,
                    nmov.id=nmov.id,nobs=nobs,log.time.prob.cond.sel=log.time.prob.cond.sel,
                    cond.sel=cond.sel)  
    llk.new=get.llk(betas=betas.new,xmat=xmat,time.prob=time.prob,mov.id=mov.id,
                    nmov.id=nmov.id,nobs=nobs,log.time.prob.cond.sel=log.time.prob.cond.sel,
                    cond.sel=cond.sel)  
    tmp=acceptMH(p0=llk.old+prior.old[i],
                 p1=llk.new+dnorm(betas.new[i],mean=0,sd=10,log=T),
                 x0=betas.old[i],x1=betas.new[i])
    if (tmp$accept==1) {betas.old[i]=tmp$x; accept1[i]=1}
  }
  llk=ifelse(accept1[ncov]==1,llk.new,llk.old)
  list(betas=betas.old,accept1=accept1,llk=llk)
}
#--------------------------------------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  cond=(accept1/accept.output)>0.4 & jump1<10000
  jump1[cond] = jump1[cond]*2       
  cond=(accept1/accept.output)<0.2 & jump1>0.001
  jump1[cond] = jump1[cond]*0.5
  
  accept1[]=0
  return(list(jump1=jump1,accept1=accept1))
}

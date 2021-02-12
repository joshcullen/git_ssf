ssf_gibbs=function(dat,ngibbs,nburn,xmat){
  #get probability time
  time.prob=dat$time.prob
  log.time.prob=log(dat$time.prob)
  
  #things useful to calculate llk
  cond.sel=dat$selected==1
  mov.id=dat$mov.id
  nmov.id=max(mov.id)
  log.time.prob.cond.sel=log.time.prob[cond.sel]
  nobs=nrow(dat)
  
  #initial values
  ncov=ncol(xmat)
  betas=rep(0,ncov)
  
  #things for gibbs
  jump1=rep(1,ncov)
  accept1=rep(0,ncov)
  accept.output=50
  store.betas=matrix(NA,ngibbs,ncov)
  store.llk=matrix(NA,ngibbs,1)
  
  for (i in 1:ngibbs){
    print(i)
    
    #sample betas
    tmp=sample.betas(betas=betas,ncov=ncov,jump1=jump1,xmat=xmat,time.prob=time.prob,
                     mov.id=mov.id,nmov.id=nmov.id,nobs=nobs,log.time.prob.cond.sel=log.time.prob.cond.sel,
                     cond.sel=cond.sel)
    betas=tmp$betas
    accept1=accept1+tmp$accept1
    
    #adapt MH algorithm
    if (i%%accept.output==0 & i<nburn){
      k=print.adapt(accept1z=accept1,jump1z=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #store results
    store.betas[i,]=betas
    store.llk[i]=tmp$llk
  }
  seq1=nburn:ngibbs
  list(betas=store.betas[seq1,],
       llk=store.llk)
}
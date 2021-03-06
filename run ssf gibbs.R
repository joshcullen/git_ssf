
library(data.table)
library(Rcpp)
library(here)
library(tictoc)

sourceCpp('aux1.cpp')
source('ssf gibbs main function.R')
source('ssf gibbs functions.R')

set.seed(123)


dat=as.data.frame(fread(here('fake data', 'cumul_time and covs.csv')))

nomes=paste0('cov',1:3)
xmat=data.matrix(dat[,nomes])

#time model details
time.int=4
gamma.b=2

#parameters for gibbs sampler
ngibbs=1000
nburn=ngibbs/2

#get probability time
time.prob=dgamma(time.int,gamma.b*dat$cum.time,gamma.b)
log.time.prob=dgamma(time.int,gamma.b*dat$cum.time,gamma.b,log=T)


for (i in 1:length(unique(dat$mov.id))) {
  print(i)
  
  ind<- which(dat$mov.id == i)
  ind.u<- which(dat$mov.id == i & dat$selected == 1)
  
  cond<- which(time.prob[ind] < (time.prob[ind.u] * 0.01))
  time.prob[ind[cond]]<- NA
  log.time.prob[ind[cond]]<- NA
}

ind<- which(is.na(time.prob))
time.prob<- time.prob[-ind]
log.time.prob<- log.time.prob[-ind]
dat<- dat[-ind,]
xmat<- xmat[-ind,]

tic()
mod1=ssf_gibbs(time.int=time.int,gamma.b=gamma.b,dat=dat,ngibbs=ngibbs,nburn=nburn,xmat=xmat)
toc()
# takes ~ 13 sec for 1000 iter

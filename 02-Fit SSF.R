
library(data.table)
library(Rcpp)
library(here)
library(tictoc)
library(raster)

sourceCpp('aux1.cpp')
source('ssf gibbs main function.R')
source('ssf gibbs functions.R')

set.seed(123)


dat=as.data.frame(fread('Giant Armadillo Time and Covs.csv'))

nomes=paste0('cov',1:2)
xmat=data.matrix(dat[,nomes])

#time model details
time.int=4
gamma.b=0.312  #mean from hierarchical JAGS model

#parameters for gibbs sampler
ngibbs=1000
nburn=ngibbs/2

#get probability time
time.prob=dgamma(time.int,gamma.b*dat$cum.time,gamma.b)
log.time.prob=dgamma(time.int,gamma.b*dat$cum.time,gamma.b,log=T)

tic()
for (i in 1:length(unique(dat$mov.id))) {
  print(i)
  
  ind<- which(dat$mov.id == i)
  ind.u<- which(dat$mov.id == i & dat$selected == 1)
  
  cond<- which(time.prob[ind] < (time.prob[ind.u] * 0.01))
  time.prob[ind[cond]]<- NA
  log.time.prob[ind[cond]]<- NA
}
toc()  #takes 27 min to filter out low prob cells

res.grid <- raster(ncol=41,nrow=41, crs="+proj=utm +units=m")
values(res.grid) <- time.prob[ind]
plot(res.grid)

ind<- which(is.na(time.prob))
time.prob<- time.prob[-ind]
log.time.prob<- log.time.prob[-ind]
dat<- dat[-ind,]
xmat<- xmat[-ind,]

tic()
mod1=ssf_gibbs(time.int=time.int,gamma.b=gamma.b,dat=dat,ngibbs=ngibbs,nburn=nburn,xmat=xmat)
toc()


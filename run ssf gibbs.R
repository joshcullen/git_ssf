rm(list=ls())
library('data.table')
library('Rcpp')

setwd('U:\\GIT_models\\git_ssf')
sourceCpp('aux1.cpp')
source('ssf gibbs main function.R')
source('ssf gibbs functions.R')

setwd('U:\\GIT_models\\git_ssf\\fake data')
dat=as.data.frame(fread('cumul_time and covs.csv'))

nomes=paste0('cov',1:3)
xmat=data.matrix(dat[,nomes])

#time model details
time.int=4
gamma.b=2

#parameters for gibbs sampler
ngibbs=1000
nburn=ngibbs/2

mod1=ssf_gibbs(time.int=time.int,gamma.b=gamma.b,dat=dat,ngibbs=ngibbs,nburn=nburn,xmat=xmat)

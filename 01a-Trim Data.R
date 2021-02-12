library(data.table)
library(tictoc)
library(raster)

dat=as.data.frame(fread('Giant Armadillo Time and Covs.csv'))

nomes=paste0('cov',1:2)
xmat=data.matrix(dat[,nomes])

#time model details
time.int=4
gamma.b=0.312  #mean from hierarchical JAGS model

#get probability time
dat$time.prob=dgamma(time.int,gamma.b*dat$cum.time,gamma.b)

#get time.prob for selected pixels
cond=dat$selected==1
aux=dat[cond,c('mov.id','time.prob')]
colnames(aux)[2]='time.prob.sel'

#merge and subset
dat1=merge(dat,aux,all=T); dim(dat); dim(dat1)
cond=dat1$time.prob>(dat1$time.prob.sel*0.01)
dat2=dat1[cond,]; dim(dat1); dim(dat2); mean(cond)

#standardize covariates
dat2$cov1=(dat2$cov1-mean(dat2$cov1))/sd(dat2$cov1)
dat2$cov2=(dat2$cov2-mean(dat2$cov2))/sd(dat2$cov2)

#export
ind=which(colnames(dat2)%in%c('x','y','cum.time','time.prob.sel'))
write.csv(dat2[,-ind],'Giant Armadillo Time and Covs trimmed.csv',row.names=F)
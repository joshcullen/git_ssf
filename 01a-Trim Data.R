library(data.table)
library(tidyverse)
library(tictoc)
library(raster)

dat=as.data.frame(fread('Giant Armadillo Time and Covs.csv'))

#Get value for gamma.b from speed model
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist")

store.betas<- read.csv("Giant Armadillo Resistance Results.csv", as.is = T)

nomes=paste0('cov',1)
xmat=data.matrix(dat[,nomes])

#time model details
time.int=4
gamma.b=mean(store.betas$gamma.b)

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
# dat2$cov1=(dat2$cov1-mean(dat2$cov1))/sd(dat2$cov1)
# dat2$cov2=(dat2$cov2-mean(dat2$cov2))/sd(dat2$cov2)


### Viz probability surface per step
#1031 total steps
sel<- dat2 %>%  #plot points for true steps
  filter(mov.id %in% 120:132) %>% 
  filter(selected == 1)
sel$mov.id<- sel$mov.id - 1

ggplot() +
  geom_tile(data = dat2 %>% filter(mov.id %in% 120:131), aes(x, y, fill = time.prob)) +
  scale_fill_viridis_c() +
  geom_path(data = sel[,-1], aes(x, y), col = "grey50") +
  geom_point(data = sel[,-1], aes(x, y), col = "grey50", shape = 1) +
  geom_point(data = sel[-1,], aes(x, y), col = "red", size = 2) +
  theme_bw() +
  facet_wrap(~mov.id)




#export
setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_ssf")

ind=which(colnames(dat2)%in%c('x','y','cum.time','time.prob.sel'))
# write.csv(dat2[,-ind],'Giant Armadillo Time and Covs trimmed.csv',row.names=F)


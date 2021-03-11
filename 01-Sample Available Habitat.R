
library(tidyverse)
library(lubridate)
library(raster)
library(akima) #for linear interpolation
library(here)

source(here('fake data', 'aux functions.R'))


# Load armadillo location data and environ covars
setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")

dat<- read.csv('Giant Armadillo state estimates.csv', as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing)
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
# dat$season<- ifelse(dat$month %in% month.abb[1:7], "Flood", "Dry")
dat.em<- dat %>% 
  filter(id == "emanuel")


# Load environ covars
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist_avg")

evi<- brick('GiantArm_evi_monthly.grd')
evi<- crop(evi, extent(dat.em %>% 
                         summarize(xmin = min(x) - 3000,
                                   xmax = max(x) + 3000,
                                   ymin = min(y) - 3000,
                                   ymax = max(y) + 3000) %>% 
                         unlist()))
evi[getValues(evi) > 1 | getValues(evi) < -1]<- NA  #mask pixels where values are outside of accepted range
evi<- evi[[which(names(evi) %in% unique(dat.em$month))]]
evi.s<- scale(evi)

# Extract cell numbers from coordinates
cell.num<- cellFromXY(evi.s[[1]], dat.em[,c("x","y")])

# Store seasonal values for each covar in an array
xmat<- array(NA, dim = c(dim(evi.s)[1:2], 1, nlayers(evi.s)))
for (i in 1:nlayers(evi.s)) {
  covs<- stack(evi.s[[i]])
  xmat[,,,i]<- as.array(covs)
}

xdim=dim(xmat)[1]
ydim=dim(xmat)[2]
ncov=dim(xmat)[3]
# nseasons=dim(xmat)[4]



# Load resistance model results
setwd("~/Documents/Snail Kite Project/Data/R Scripts/ValleLabUF/resist")

resist<- read.csv('Giant Armadillo Resistance Surfaces.csv', as.is = T)
resist1<- pivot_wider(resist, id_cols = c(x,y), names_from = month, values_from = time)
mean.time.all<- rasterFromXYZ(resist1)
mean.time<- as.array(mean.time.all)  #reduce for single season and convert to matrix


#important variables
nsim=length(cell.num)
window1=15
ndados=((window1*2)+1)^2
store.calc=matrix(NA,(nsim-1)*ndados,5+ncov)
oooo=1
months1<- names(evi.s)


for (i in 1:(nsim-1)){  #extracts covar values per season
  print(i)
  
  #convert cell number to row and col coords
  coord<- as.vector(rowColFromCell(evi.s[[1]], cell.num[i]))
  names(coord)<- c("y","x")  #y=row,  x=column
  cond<- which(months1 %in% dat.em$month[i])
  
  #calculate cumulative time across 8 directions
  res=get.cumtime.8dir(coord1=coord,mean.time=mean.time[,,cond],window1=window1)
  
  #calculate prob of staying in the sample place
  tmp=matrix(coord,1,2)
  tempo.not.move=mean.time[tmp[1], tmp[2], cond]
  res1=rbind(res,c(tempo.not.move,coord[c('x','y')]))
  
  #interpolate cumulative time across a window1 x window1 area
  minx=coord['x']-window1; maxx=coord['x']+window1
  miny=coord['y']-window1; maxy=coord['y']+window1
  res2 = interp(x=res1[,'x'],y=res1[,'y'],z=res1[,'cum.time'],
                xo = minx:maxx, yo = miny:maxy, 
                linear=T, extrap = F)
  res3=interp2xyz(res2)
  colnames(res3)[3]='cum.time'
  
  #get covariates
  aux=numeric()
  for (j in 1:ncov){
    tmp=cbind(res3[,c('y','x')],j,cond)
    tmp1=xmat[tmp]
    aux=cbind(aux,tmp1)
  }
  colnames(aux)=paste0('cov',1:ncov)
  res4=cbind(res3,aux)
  
  #store results  
  ind=which(res4[,'x']==coord['x'] & res4[,'y']==coord['y'])
  zeros=rep(0,nrow(res4))
  zeros[ind]=1
  
  seq1=oooo:(oooo+ndados-1)
  store.calc[seq1,]=cbind(res4,zeros,i)  
  oooo=oooo+ndados
}


colnames(store.calc)=c(colnames(res4),'selected','mov.id')


# Export data
setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_ssf")

# write.csv(store.calc, 'Giant Armadillo Time and Covs.csv', row.names=F)

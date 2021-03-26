
library(lubridate)
library(tidyverse)
library(raster)


# Load armadillo location data and environ covars
setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")

dat<- read.csv('Giant Armadillo state estimates.csv', as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing)
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])

dat.em<- dat %>% 
  filter(id == "emanuel")


#load EVI data
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


# Load beta coefficients from model
setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_ssf")

betas<- read.csv("SSF coefficients.csv", as.is = T)




### Calculate SSF Surface ###

ssfSurf<- vector("list", 5)
names(ssfSurf)<- names(evi.s)

for (j in 1:nlayers(evi.s)) {
  print(names(evi.s)[j])
  
  #create matrix of covars
  cov.mat<- cbind(evi = raster::values(evi.s[[j]]))
  
  #calc habitat preference
  ssf.res<- evi.s[[1]]
  w.hat<- exp(cov.mat %*% betas$mean)
  raster::values(ssf.res)<- w.hat/(1 + w.hat)
  
  #create as data frame
  ssf.res.df<- as.data.frame(ssf.res, xy=T)
  names(ssf.res.df)[3]<- "sel"
  
  ssfSurf[[j]]<- ssf.res.df
   
}


ssfSurf.df<- bind_rows(ssfSurf, .id = "month")
ssfSurf.df$month<- factor(ssfSurf.df$month, levels = unique(ssfSurf.df$month))


# Plot seasonal predictive surfaces
ggplot() +
  geom_raster(data = ssfSurf.df, aes(x, y, fill = sel)) +
  geom_path(data = dat.em, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
  scale_fill_viridis_c("Selection", option = "inferno",
                       na.value = "transparent", limits = c(0,1)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="Easting", y="Northing", title = "Habitat Selection") +
  theme_bw() +
  coord_equal() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(fill = guide_colourbar(barwidth = 30, barheight = 1)) +
  facet_wrap( ~ month)
  
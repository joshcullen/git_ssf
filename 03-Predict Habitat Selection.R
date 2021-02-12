
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
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))


green<- brick('GiantArm_tcgreen_season.grd')
wet<- brick('GiantArm_tcwet_season.grd')

# scaled covars
green_s<- scale(green)
wet_s<- scale(wet)


# Load beta coefficients from model
setwd("~/Documents/Snail Kite Project/Data/R Scripts/git_ssf")

betas<- read.csv("SSF coefficients.csv", as.is = T)




### Calculate SSF Surface ###

ssfSurf<- vector("list", 4)
names(ssfSurf)<- names(green_s)

for (j in 1:nlayers(green_s)) {
  print(names(green_s)[j])
  
  #create matrix of covars
  cov.mat<- cbind(green = raster::values(green_s[[j]]), wet = raster::values(wet_s[[j]]))
  
  #calc habitat preference
  ssf.res<- green[[1]]
  w.hat<- exp(cov.mat %*% betas$mean)
  raster::values(ssf.res)<- w.hat/(1 + w.hat)
  
  #create as data frame
  ssf.res.df<- as.data.frame(ssf.res, xy=T)
  names(ssf.res.df)[3]<- "sel"
  
  ssfSurf[[j]]<- ssf.res.df
   
}


ssfSurf.df<- bind_rows(ssfSurf, .id = "season")
ssfSurf.df$season<- factor(ssfSurf.df$season, levels = unique(ssfSurf.df$season))


# Plot seasonal predictive surfaces
ggplot() +
  geom_raster(data = ssfSurf.df, aes(x, y, fill = sel)) +
  geom_path(data = dat, aes(x, y, group = id), alpha = 0.75, color = "chartreuse") +
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
  facet_wrap( ~ season, nrow = 2)
  
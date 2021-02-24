
library(data.table)
library(Rcpp)
library(tictoc)
library(tidyverse)
library(fishualize)
library(coda)

sourceCpp('aux1.cpp')
source('ssf gibbs main function.R')
source('ssf gibbs functions.R')

set.seed(123)


#################
### Load Data ###
#################

dat=as.data.frame(fread('Giant Armadillo Time and Covs trimmed.csv'))
nomes=paste0('cov',1:2)
xmat=data.matrix(dat[,nomes])


#################
### Run Model ###
#################

#parameters for gibbs sampler
ngibbs=1000
nburn=ngibbs/2

#get probability time
time.prob=dat$time.prob
log.time.prob=log(time.prob)

tic()
mod1=ssf_gibbs(dat=dat,ngibbs=ngibbs,nburn=nburn,xmat=xmat)
toc()  #takes 19 min to run 1000 iterations


#######################
### Inspect Results ###
#######################

plot(mod1$llk,type='l')
seq1=nburn:length(mod1$llk)
plot(mod1$llk[seq1],type='l')

par(mfrow=c(2,1))
for (i in 1:2) plot(mod1$betas[,i],type='l')

apply(mod1$betas,2,mean)


### Make (pretty) caterpillar plot
betas.post<- data.frame(mod1$betas)
names(betas.post)<- c("ndvi", "ndwi")
betas<- betas.post %>% 
  pivot_longer(., cols = c(ndvi, ndwi), names_to = "coeff", values_to = "value") %>% 
  group_by(coeff) %>% 
  summarize(mean = mean(value)) %>% 
  ungroup()
hpdi<- as.mcmc(betas.post) %>% 
  HPDinterval() %>% 
  data.frame()
betas<- cbind(betas, hpdi)


ggplot(data=betas, aes(x=coeff, y=mean, ymin=lower, ymax=upper, color=coeff)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(position = position_dodge(0.55), width = 0, size = 0.75) +
  geom_point(position = position_dodge(0.55), size=2) +
  scale_x_discrete(labels = c("NDVI", "NDWI")) +
  scale_color_fish_d("", option = "Scarus_tricolor") +
  coord_flip() +
  theme_bw() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())


##################################
### Viz partial response plots ###
##################################

## NDVI

#Generate sequence along green
rango1<- dat %>% 
  filter(selected == 1) %>% 
  dplyr::select(cov1) %>% 
  range()
seq.ndvi<- seq(rango1[1], rango1[2], length.out = 100)


#Create design matrix where 0s added for all other vars besides green
design.mat<- cbind(seq.ndvi, 0)

# Take cross-product of design matrix with betas and exponentiate to calc response
y.mu<- exp(design.mat %*% betas$mean)  
y.low<- exp(design.mat %*% betas$lower)
y.up<- exp(design.mat %*% betas$upper)


# Add results to data frame
ndvi.mu.df<- data.frame(x = seq.ndvi,
                     y = y.mu,
                     ymin = y.low,
                     ymax = y.up)


# Plot relationship
ggplot(data = ndvi.mu.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), fill = "forestgreen", alpha =  0.3) +
  geom_line(aes(x, y), color = "forestgreen", size = 1) +
  labs(x = "\nStandardized NDVI", y = "Habitat Preference\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))




## NDWI

#Generate sequence along wet
rango1<- dat %>% 
  filter(selected == 1) %>% 
  dplyr::select(cov2) %>% 
  range()
seq.ndwi<- seq(rango1[1], rango1[2], length.out = 100)


#Create design matrix where 0s added for all other vars besides ndwi
design.mat<- cbind(0, seq.ndwi)

# Take cross-product of design matrix with betas and exponentiate to calc response
y.mu<- exp(design.mat %*% betas$mean)  #inverse logit
y.low<- exp(design.mat %*% betas$lower)
y.up<- exp(design.mat %*% betas$upper)


# Add results to data frame
ndwi.mu.df<- data.frame(x = seq.ndwi,
                         y = y.mu,
                         ymin = y.low,
                         ymax = y.up)


# Plot relationship
ggplot(data = ndwi.mu.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), fill = "cadetblue", alpha =  0.3) +
  geom_line(aes(x, y), color = "cadetblue", size = 1) +
  labs(x = "\nStandardized Wetness", y = "Habitat Preference\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))


######################
### Export Results ###
######################

# write.csv(betas, "SSF coefficients.csv", row.names = F)

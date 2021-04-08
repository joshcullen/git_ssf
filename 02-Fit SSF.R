
library(data.table)
library(Rcpp)
library(tictoc)
library(tidyverse)
library(fishualize)
library(coda)
library(splines)

sourceCpp('aux1.cpp')
source('ssf gibbs main function.R')
source('ssf gibbs functions.R')

set.seed(123)


#################
### Load Data ###
#################

dat=as.data.frame(fread('Giant Armadillo Time and Covs trimmed.csv'))
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = unique(dat$month))

# Add B-spline (w/ 2 internal knots) for 'EVI'
rango<- range(dat$cov1)
knot.locs<- seq(rango[1], rango[2], length.out = 4)[2:3]
spline.evi<- as.data.frame(bs(dat$cov1, degree=2, intercept = FALSE,
                              knots = knot.locs))
names(spline.evi)<- paste("spline", 1:ncol(spline.evi), sep = ".")
dat<- cbind(dat, spline.evi)

#what proportion of rows have sums < 1?
foo<- rowSums(spline.evi)
length(which(foo < 1))/length(foo)  #9.97%
length(which(foo < 0.95))/length(foo)  #0.39%

#create dummy variable for month
month.dumm<- model.matrix(~dat$month + 0)
 

ind<- c(paste("spline", 1:ncol(spline.evi), sep = "."))
# ind<- "cov1"
# xmat=data.matrix(dat[,ind])
xmat<- data.matrix(cbind(month.dumm[,-1], dat[,ind]))  #treat May as ref
colnames(xmat)[1:4]<- c("Jun","Sep","Oct","Nov")



#################
### Run Model ###
#################

#parameters for gibbs sampler
ngibbs=2000
nburn=ngibbs/2

#get probability time
time.prob=dat$time.prob
log.time.prob=log(time.prob)

tic()
mod1=ssf_gibbs(dat=dat,ngibbs=ngibbs,nburn=nburn,xmat=xmat)
toc()  #takes 50 min to run 5000 iterations


#######################
### Inspect Results ###
#######################

plot(mod1$llk,type='l')
seq1=nburn:length(mod1$llk)
plot(mod1$llk[seq1],type='l')

par(mfrow=c(2,2))
for (i in 1:ncol(mod1$betas)) plot(mod1$betas[,i],type='l')
# plot(mod1$betas, type = "l")

# apply(mod1$betas,2,mean)
mean(mod1$betas)

### Make (pretty) caterpillar plot
betas.post<- data.frame(mod1$betas)
names(betas.post)<- c("EVI")
betas<- betas.post %>% 
  pivot_longer(., cols = c(EVI), names_to = "coeff", values_to = "value") %>% 
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
  # scale_x_discrete(labels = c("EVI")) +
  scale_color_fish_d("", option = "Scarus_tricolor") +
  coord_flip() +
  theme_bw() +
  labs(x="", y="") +
  theme(axis.text = element_text(size = 14),
        panel.grid = element_blank())


##################################
### Viz partial response plots ###
##################################

## EVI

#Generate sequence along green
rango1<- dat %>% 
  filter(selected == 1) %>% 
  dplyr::select(cov1) %>% 
  range()
seq.evi<- seq(rango1[1], rango1[2], length.out = 100)


#Create design matrix where 0s added for all other vars besides green
design.mat<- cbind(seq.evi)

# Take cross-product of design matrix with betas and exponentiate to calc response
y.mu<- exp(design.mat %*% betas$mean)  
y.low<- exp(design.mat %*% betas$lower)
y.up<- exp(design.mat %*% betas$upper)


# Add results to data frame
evi.mu.df<- data.frame(x = seq.evi,
                     y = y.mu,
                     ymin = y.low,
                     ymax = y.up)


# Plot relationship
ggplot(data = evi.mu.df) +
  geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), fill = "forestgreen", alpha =  0.3) +
  geom_line(aes(x, y), color = "forestgreen", size = 1) +
  labs(x = "\nStandardized EVI", y = "Habitat Preference\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12))






######################
### Export Results ###
######################

# write.csv(betas, "SSF coefficients.csv", row.names = F)

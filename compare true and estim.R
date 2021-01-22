plot(mod1$llk,type='l')
seq1=100:length(mod1$llk)
plot(mod1$llk[seq1],type='l')

par(mfrow=c(3,1))
for (i in 1:3) plot(mod1$betas[,i],type='l')

apply(mod1$betas,2,mean)
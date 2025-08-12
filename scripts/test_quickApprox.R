library(RatingShiftHappens)
# Verify type-1 error level
nSim=100
res=rep(NA,nSim)
n=100
for (i in 1:nSim){
  obs=c(rnorm(n),0.+1.*rnorm(n))
  foo=Segmentation(obs,varShift=TRUE,alpha=0.1)
  res[i]=foo$nS
}
plot(obs)
mean(res==2)

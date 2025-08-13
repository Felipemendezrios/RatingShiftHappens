library(RatingShiftHappens)
set.seed(54321)
#---------------------------
# Verify type-1 error level. Since the quick approximation is based on a formal test,
# it should be roughly respected
nSim=100
res=rep(NA,nSim)
n=80
for (i in 1:nSim){
  obs=rnorm(n)
  foo=Segmentation(obs,alpha=0.1)
  res[i]=foo$nS
}
mean(res==2) # rejection rate

#---------------------------
# Compare quick approx and full treatment on contrasted cases
compare <- function(obs,time=1:length(obs),u=0*obs,doRecursive=FALSE,varShift=FALSE){
  if(doRecursive){
    quick=Recursive_Segmentation(obs,time,u,varShift=varShift)
    full=Recursive_Segmentation(obs,time,u,doQuickApprox=FALSE,stout=NULL)
  } else {
    quick=Segmentation(obs,time,u,varShift=varShift)
    full=Segmentation(obs,time,u,doQuickApprox=FALSE,stout=NULL)
  }
  message('Detected shifts [quick]:')
  print(quick$summary$shift)
  message('Detected shifts [full]:')
  print(full$summary$shift)
  par(mfrow=c(1,2))
  plot(obs)
  if(doRecursive){
    plot(full$summary$shift$id_iteration,full$summary$shift$tau,
         ylim=range(full$summary$data$time),
         xlim=c(0,1+max(c(max(full$summary$shift$id_iteration),max(full$summary$shift$id_iteration)))))
    points(full$summary$shift$id_iteration,full$summary$shift$I95_lower,pch=3)
    points(full$summary$shift$id_iteration,full$summary$shift$I95_upper,pch=3)
    points(quick$summary$shift$id_iteration,quick$summary$shift$tau,col='red')
    points(quick$summary$shift$id_iteration,quick$summary$shift$I95_lower,pch=3,col='red')
    points(quick$summary$shift$id_iteration,quick$summary$shift$I95_upper,pch=3,col='red')
  } else {
    hist(full$results[[2]]$mcmc$tau1,50,freq=FALSE,main='',xlab='tau')
    hist(quick$results[[2]]$mcmc$tau1,50,freq=FALSE,col=rgb(1,0,0,0.5),add=TRUE,alpha=0.5)
  }
}

# middle change
obs=c(rnorm(50),1+rnorm(50))
compare(obs)
# early change
obs=c(rnorm(10),2+rnorm(50))
compare(obs)
# many changes
obs=c(rnorm(10),3+rnorm(10),rnorm(10),-4+rnorm(10),-2+rnorm(10),rnorm(10))
compare(obs,doRecursive=TRUE)
# Detect varying variance? (only available with quick-approx)
obs=c(rnorm(50),2*rnorm(50))
compare(obs,varShift=TRUE)
# What if the shifting variance corresponds to more observation uncertainty?
compare(obs,u=c(rep(1,50),rep(2,50)),varShift=TRUE)

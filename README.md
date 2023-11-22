Rating Shift Happens
================
Felipe MENDEZ and Benjamin RENARD (INRAE, RiverLy and RECOVER).
December 2023

## Introduction

Darienzo (2021) developed
[BayDERS](https://github.com/MatteoDarienzo/BayDERS) for detecting and
estimating stage-discharge rating shifts for retrospective and real-time
streamflow quantification. In 2022, a sensibility analysis was conducted
to assess the stability and performance of the tool, leading to multiple
corrections. Nevertheless, debugging the tool proved difficult due to
the integration of different models such as rating curve and recessions
estimation. Hence, a dissociation procedure has been implemented to
extract essential functions, enabling the automatic segmentation of a
random variable. Theses functions have been organized and coded in a
package format, which will be explained below.

The goal of `RatingShiftHappens` package is to create a tools package
for detecting, visualizing and estimating rating shifts. This
documentation provides the description of several functions available to
the segmentation process.

Three fundamental functions are available so far :

1.  segmentation.engine
2.  segmentation
3.  recursive.segmentation

## Installation

You can install the development version of
[RatingShiftHappens](https://github.com/Felipemendezrios/RatingShiftHappens)
from GitHub using the following command. Please note that the
[RBaM](https://github.com/BaM-tools/RBaM) package is also required to
run this package.

``` r
# install.packages("devtools")
devtools::install_github("Felipemendezrios/RatingShiftHappens")
devtools::install_github('BaM-tools/RBaM') 

library(RatingShiftHappens)
```

Functions will be explained more precisely below along with a example.

## segmentation.engine function

Segmentation procedure for a **known** given number of segments

This is a basic example which shows you how to segment a random
variable:

``` r
 # Set random generation
 set.seed(1)

 # Create observation vector
 obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
 nS.user=2

 # Run segmentation engine function
 res <- segmentation.engine(obs=obs,nS=nS.user)

 # Estimated shift time
 res$tau
#> [1] 25.1351

 # intervals defined by time shifts
 if(nS.user!=1){
  intervals.time.shift=c(res$data.p$time[[1]][1],res$tau,rev(res$data.p$time[[nS.user]])[1])
 }else{
  intervals.time.shift=list(res$data.p$time[1],rev(res$data.p$time)[1])
 }

 # Maximum a posterior value per segment indexed by the list number
 res$segments
#> [[1]]
#>  [1] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> [10] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> [19] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> 
#> [[2]]
#>  [1] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559
#> [10] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559
#> [19] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559

 # Uncertainty in shift time
 if(nS.user!=1){
    Shift=res$mcmc$tau1
    hist(Shift,
         main='Histogram of first shift')

    uncertainty95_shift <- list()
    for(i in 1:(nS.user-1)){
     uncertainty95_shift[[i]] = stats::quantile(res$mcmc[,nS.user+i],probs=c(0.025,0.975))
    }
 }
```

<img src="man/readme/README-segmentation.engine-1.png" width="100%" />

``` r

 # Uncertainty in segment estimation
 mu.seg.1.unc=res$mcmc$mu1
 graphics::hist(mu.seg.1.unc,
                xlab='obs',
                main='Histogram of first segment of observation')
```

<img src="man/readme/README-segmentation.engine-2.png" width="100%" />

``` r
 uncertainty95_segment <- list()
 for(i in 1:nS.user){
    uncertainty95_segment [[i]] = stats::quantile(res$mcmc[,i],probs=c(0.025,0.975))
 }

 # Separate and assign information by identified stable period
 res$data.p
#> $obs.p
#> $obs.p[[1]]
#>  [1] -0.62645381  0.18364332 -0.83562861  1.59528080  0.32950777 -0.82046838
#>  [7]  0.48742905  0.73832471  0.57578135 -0.30538839  1.51178117  0.38984324
#> [13] -0.62124058 -2.21469989  1.12493092 -0.04493361 -0.01619026  0.94383621
#> [19]  0.82122120  0.59390132  0.91897737  0.78213630  0.07456498 -1.98935170
#> [25]  0.61982575
#> 
#> $obs.p[[2]]
#>  [1] 1.9438713 1.8442045 0.5292476 1.5218499 2.4179416 3.3586796 1.8972123
#>  [8] 2.3876716 1.9461950 0.6229404 1.5850054 1.6057100 1.9406866 3.1000254
#> [15] 2.7631757 1.8354764 1.7466383 2.6969634 2.5566632 1.3112443 1.2925048
#> [22] 2.3645820 2.7685329 1.8876538 2.8811077
#> 
#> 
#> $time.p
#> $time.p[[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> 
#> $time.p[[2]]
#>  [1] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
#> 
#> 
#> $u.p
#> $u.p[[1]]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $u.p[[2]]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
 
 # DIC estimation
 res$DIC
#> [1] 131.985

 # Setting transparency for plotting
 alpha = 125

 # Set color plot
 color_customized_rect <- function(alpha){
    color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
                   rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
                   rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
                   rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
                   rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
    return(color)
 }

 # Assign period to data
 obs_segmented <- data.frame()

 # Conditional to separate non segmentation case
 if(typeof(res$data.p$obs.p)=='list'){
   for(i in 1:length(res$data.p$obs.p)){
    obs_segmented_temp=cbind(obs=res$data.p$obs.p[[i]],period=i)
    obs_segmented=rbind(obs_segmented,obs_segmented_temp)
   }
 }else{
  obs_segmented=data.frame(obs=res$data.p$obs.p,period=1)
 }

 # Plot observations
 plot(x=obs_segmented$obs,
      col=factor(obs_segmented$period),
      pch=16,
      main='Final segmentation',
      ylab='obs',
      xlab='time')

 # Plot segments
 for(i in 1:nS.user){
   segments(x0=intervals.time.shift[[i]],
            x1=intervals.time.shift[[i+1]],
            y0=res$segments[[i]],
            y1=res$segments[[i]],
            col='blue')
   rect(xleft=intervals.time.shift[[i]],
        xright=intervals.time.shift[[i+1]],
        ybottom=uncertainty95_segment[[i]][1],
        ytop=uncertainty95_segment[[i]][2],
        col= rgb(0,0,255,max=255,alpha=125,names='blue'),
        border = 'transparent')
 }

 # Plot shifts
 if(nS.user!=1){
  for(i in 1:(nS.user-1)){
   abline(v=res$tau[i],col=color_customized_rect(255)[[i]], lwd=2)
   rect(xleft=uncertainty95_shift[[i]][1],
        xright=rev(uncertainty95_shift[[i]])[1],
        ybottom=min(obs)*2,
        ytop=max(obs)*2,
        col= color_customized_rect(125)[[i]],
        border = 'transparent')
   }
 }
```

<img src="man/readme/README-segmentation.engine-3.png" width="100%" />

## segmentation function

Segmentation procedure for a **unknown** given number of segments

This is a basic example which shows you how to segment a variable with a
**unknown** number of segments :

``` r
 # Set random generation
 set.seed(1)

 # Create observation vector
 obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
 
 # Set the maximum number of segments 
 nSmax.user=3

 # Run segmentation function
 res <- segmentation(obs=obs,nSmax=nSmax.user)

 # Optimal number of segments nSopt
 nSopt <- res$nS
 nSopt
#> [1] 2

 # Estimated shift time
 res$results[[nSopt]]$tau
#> [1] 25.1351

 # intervals defined by time shifts
 if(nSopt!=1){
  intervals.time.shift=c(res$results[[nSopt]]$data.p$time[[1]][1],
                         res$results[[nSopt]]$tau,
                         rev(res$results[[nSopt]]$data.p$time[[nSopt]])[1])
 }else{
  intervals.time.shift=list(res$results[[nSopt]]$data.p$time[1],
                            rev(res$results[[nSopt]]$data.p$time)[1])
 }

 # Maximum a posterior value per segment indexed by the list number
 res$results[[nSopt]]$segments
#> [[1]]
#>  [1] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> [10] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> [19] 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329 0.18329
#> 
#> [[2]]
#>  [1] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559
#> [10] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559
#> [19] 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559 2.06559

 # Uncertainty in shift time
 if(nSopt!=1){
    Shift=res$results[[nSopt]]$mcmc$tau1
    hist(Shift,
         main='Histogram of first shift')

    uncertainty95_shift <- list()
    for(i in 1:(nSopt-1)){
     uncertainty95_shift[[i]] = stats::quantile(res$results[[nSopt]]$mcmc[,nSopt+i],
                                                probs=c(0.025,0.975))
    }
 }
```

<img src="man/readme/README-segmentation-1.png" width="100%" />

``` r
 # Uncertainty in segment estimation
 mu.seg.1.unc=res$results[[nSopt]]$mcmc$mu1
 graphics::hist(mu.seg.1.unc,
                xlab='obs',
                main='Histogram of first segment of observation')
```

<img src="man/readme/README-segmentation-2.png" width="100%" />

``` r
 uncertainty95_segment <- list()
 for(i in 1:nSopt){
    uncertainty95_segment [[i]] = stats::quantile(res$results[[nSopt]]$mcmc[,i],
                                                  probs=c(0.025,0.975))
 }
 # Separate and assign information by identified stable period
 res$results[[nSopt]]$data.p
#> $obs.p
#> $obs.p[[1]]
#>  [1] -0.62645381  0.18364332 -0.83562861  1.59528080  0.32950777 -0.82046838
#>  [7]  0.48742905  0.73832471  0.57578135 -0.30538839  1.51178117  0.38984324
#> [13] -0.62124058 -2.21469989  1.12493092 -0.04493361 -0.01619026  0.94383621
#> [19]  0.82122120  0.59390132  0.91897737  0.78213630  0.07456498 -1.98935170
#> [25]  0.61982575
#> 
#> $obs.p[[2]]
#>  [1] 1.9438713 1.8442045 0.5292476 1.5218499 2.4179416 3.3586796 1.8972123
#>  [8] 2.3876716 1.9461950 0.6229404 1.5850054 1.6057100 1.9406866 3.1000254
#> [15] 2.7631757 1.8354764 1.7466383 2.6969634 2.5566632 1.3112443 1.2925048
#> [22] 2.3645820 2.7685329 1.8876538 2.8811077
#> 
#> 
#> $time.p
#> $time.p[[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> 
#> $time.p[[2]]
#>  [1] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
#> 
#> 
#> $u.p
#> $u.p[[1]]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> 
#> $u.p[[2]]
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

 # DIC estimation
 res$results[[nSopt]]$DIC
#> [1] 131.985

 # Setting plot

 # Transparency
 alpha <- 125

 # Set color plot
 color_customized_rect <- function(alpha){
    color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
                   rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
                   rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
                   rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
                   rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
    return(color)
 }

 # Assign period to data
 obs_segmented <- data.frame()

 # Conditional to separate non segmentation case
 if(typeof(res$results[[nSopt]]$data.p$obs.p)=='list'){
   for(i in 1:length(res$results[[nSopt]]$data.p$obs.p)){
    obs_segmented_temp=cbind(obs=res$results[[nSopt]]$data.p$obs.p[[i]],period=i)
    obs_segmented=rbind(obs_segmented,obs_segmented_temp)
   }
 }else{
  obs_segmented=data.frame(obs=res$results[[nSopt]]$data.p$obs.p,period=1)
 }

 # Plot observations
 plot(x=obs_segmented$obs,
      col=factor(obs_segmented$period),
      pch=16,
      main='Final segmentation',
      ylab='obs',
      xlab='time')

 # Plot segments
 for(i in 1:nSopt){
   segments(x0=intervals.time.shift[[i]],
            x1=intervals.time.shift[[i+1]],
            y0=res$results[[nSopt]]$segments[[i]],
            y1=res$results[[nSopt]]$segments[[i]],
            col='blue')
   rect(xleft=intervals.time.shift[[i]],
        xright=intervals.time.shift[[i+1]],
        ybottom=uncertainty95_segment[[i]][1],
        ytop=uncertainty95_segment[[i]][2],
        col= rgb(0,0,255,max=255,alpha=125,names='blue'),
        border = 'transparent')
 }

 # Plot shifts
 if(nSopt!=1){
  for(i in 1:(nSopt-1)){
   abline(v=res$results[[nSopt]]$tau[i],col=color_customized_rect(255)[[i]], lwd=2)
   rect(xleft=uncertainty95_shift[[i]][1],
        xright=rev(uncertainty95_shift[[i]])[1],
        ybottom=min(obs)*2,
        ytop=max(obs)*2,
        col= color_customized_rect(125)[[i]],
        border = 'transparent')
   }
 }
```

<img src="man/readme/README-segmentation-3.png" width="100%" />

## recursive.segmentation function

Recursive segmentation procedure for a **unknown** given number of
segments

This is a basic example which shows you how to segment a variable with a
**unknown** number of segments using a recursive process:

``` r
# Set random generation
set.seed(1)
# Create series to be segmented
# obs=c(rnorm(30,mean=0,sd=1),rnorm(30,mean=2,sd=1))
# obs=c(rnorm(30,mean=0,sd=1),
#       rnorm(30,mean=2,sd=1),
#       rnorm(30,mean=1,sd=1),
#       rnorm(30,mean=2,sd=1))

obs=c(rnorm(30,mean=0.5,sd=1),
      rnorm(30,mean=1.5,sd=1),
      rnorm(30,mean=1,sd=1),
      rnorm(30,mean=2,sd=1))

time=1:length(obs)

# Assign a maximum number of user-defined segments to be assessed
nSmax.user=3

# Apply recursive segmentation
results=recursive.segmentation(obs, nSmax=nSmax.user)

# Have a look at recursion tree
results$tree
#>   indx level parent nS
#> 1    1     1      0  3
#> 2    2     2      1  1
#> 3    3     2      1  2
#> 4    4     2      1  1
#> 5    5     3      3  1
#> 6    6     3      3  1

# Visualize tree with data.tree package
# if(NROW(results$tree)>1){
#   tree <- data.tree::as.Node(results$tree[-1,c(3,1)],mode = "network")
#   plot(tree)
# } else { # No segmentation took place, make a dummy plot
#   tree <- data.tree::as.Node(data.frame(1,2),mode = "network")
#   plot(tree)
# }

# Get terminal nodes
terminal=which(results$tree$nS==1)

# Get node with time shifts
nodes.shift.time=which(results$tree$nS!=1)

nodes.shift.time
#> [1] 1 3

# Estimated shift time along with uncertainties
shift.time.list <- c()
for(i in 1:length(nodes.shift.time)){
  nSopt.p = results$res[[nodes.shift.time[[i]]]]$nS
  results.p=results$res[[nodes.shift.time[[i]]]]$results
  # shift.time.p=results.p[[nSopt.p]]$tau
  shift.time.p=cbind(c(results.p[[nSopt.p]]$tau))

  for(j in 1:(nSopt.p-1)){

    shift.time.p.unc=data.frame(tau=as.numeric(shift.time.p[j,]),
                                u2.5=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                     probs=c(0.025)),
                                u97.5=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                      probs=c(0.975)))
    shift.time.list <- rbind(shift.time.list,
                         shift.time.p.unc)
  }
}

all.shift.time <- shift.time.list[order(shift.time.list$tau),]
all.shift.time
#>           tau     u2.5     u97.5
#> 2.5%  29.5616 18.93388  39.98611
#> 2.5%2 63.2496 32.13913  88.11706
#> 2.5%1 90.1033 81.52526 105.78840

# intervals defined by time shifts
intervals.time.shift=c(results$res[[1]]$results[[1]]$data.p$time.p[1],
                       all.shift.time$tau,
                       rev(results$res[[1]]$results[[1]]$data.p$time.p)[1])

# Get stable periods by adding information as period, segment,
data.stable <- c()
for(i in 1:length(terminal)){
  data.stable.p=results$res[[terminal[i]]]$results[[1]]   #Save data from stable period

  node = data.frame(obs=data.stable.p$data.p$obs.p,
                    time=data.stable.p$data.p$time.p,
                    u=data.stable.p$data.p$u.p,
                    period = rep(i,length(data.stable.p$data.p$obs.p)))
  data.stable = rbind(data.stable,node)
}

# Setting plot

# Transparency
alpha <- 125

# Set color plot
color_customized_rect <- function(alpha){
  color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
                 rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
                 rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
                 rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
                 rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
  return(color)
}

# Plot observations
plot(x=data.stable$time,
y=data.stable$obs,
col=data.stable$period,
xlab='time',
ylab='obs',
main='Final segmentation'
)

# Plot shifts
for(i in 1:nrow(all.shift.time)){
  abline(v=all.shift.time$tau[i],
         col=color_customized_rect(255)[[i]])
  rect(xleft=all.shift.time$u2.5[[i]],
       xright=all.shift.time$u97.5[[i]],
       ybottom=min(obs)*2,
       ytop=max(obs)*2,
       col= color_customized_rect(50)[[i]],
       border = 'transparent')
}

# Pre-treatment data to extract segments from corresponding node
all.global.data <- c()
for(i in 1:length(nodes.shift.time)){
  nSopt.p = results$res[[nodes.shift.time[[i]]]]$nS
  results.p=results$res[[nodes.shift.time[[i]]]]$results

  relation.parent.children.p=results$tree[which(results$tree$parent==nodes.shift.time[[i]]),]
  id.segment.p=which(relation.parent.children.p$nS==1)

  for(j in 1:length(id.segment.p)){
    segments.p=results.p[[nSopt.p]]$segments[[id.segment.p[[j]]]]
    obs.p=results.p[[nSopt.p]]$data.p$obs.p[[id.segment.p[[j]]]]
    time.p=results.p[[nSopt.p]]$data.p$time.p[[id.segment.p[[j]]]]
    u.p=results.p[[nSopt.p]]$data.p$u.p[[id.segment.p[[j]]]]
    unc.segment.p=stats::quantile(results.p[[nSopt.p]]$mcmc[,id.segment.p[[j]]],
                                  probs=c(0.025,0.975))

    data.p.temp=data.frame(mu=segments.p,
                           obs=obs.p,
                           u=u.p,
                           time=time.p,
                           u2.5=rep(unc.segment.p[1],length(segments.p)),
                           u97.5=rep(unc.segment.p[2],length(segments.p)))

    all.global.data=rbind(all.global.data,data.p.temp)
  }
}

all.global.data=all.global.data[order(all.global.data$time),]

final.data.segmented=unname((split(all.global.data,all.global.data$mu)))

# Plot segments
for(i in 1:(nrow(all.shift.time)+1)){
  text(final.data.segmented[[i]]$time,
       final.data.segmented[[i]]$obs,
       terminal[i],pos=3,
       cex=0.8,
       col='black')
  segments(x0=intervals.time.shift[[i]],
           x1=intervals.time.shift[[i+1]],
           y0=final.data.segmented[[i]]$mu,
           y1=final.data.segmented[[i]]$mu,
           col='blue')
  rect(xleft=intervals.time.shift[[i]],
       xright=intervals.time.shift[[i+1]],
       ybottom=final.data.segmented[[i]]$u2.5,
       ytop=final.data.segmented[[i]]$u97.5,
       col= rgb(0,0,255,max=255,alpha=5,names='blue'),
       border = 'transparent')

}
```

<img src="man/readme/README-unnamed-chunk-3-1.png" width="100%" />

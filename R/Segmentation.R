#' Segmentation engine
#'
#' Segmentation procedure for a \strong{known} given number of segments
#'
#' @param obs real vector, observations
#' @param time real vector, time
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nS integer, number of segments
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#'
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#' }
#' @examples
#' # Set random generation
#' set.seed(1)
#'
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Assign a number of user-defined segments to be assessed
#' nS.user=2
#'
#' # Run segmentation engine function
#' res <- segmentation.engine(obs=obs,nS=nS.user)
#'
#' # Estimated shift time
#' res$tau
#'
#' # intervals defined by time shifts
#' if(nS.user!=1){
#'  intervals.time.shift=c(res$data.p$time[[1]][1],res$tau,rev(res$data.p$time[[nS.user]])[1])
#' }else{
#'  intervals.time.shift=list(res$data.p$time[1],rev(res$data.p$time)[1])
#' }
#'
#' # Maximum a posterior value per segment indexed by the list number
#' res$segments
#'
#' # Uncertainty in shift time
#' if(nS.user!=1){
#'    Shift=res$mcmc$tau1
#'    graphics::hist(Shift,
#'                   main='Histogram of first shift')
#'
#'    uncertainty95.shift <- list()
#'    for(i in 1:(nS.user-1)){
#'     uncertainty95.shift[[i]] = stats::quantile(res$mcmc[,nS.user+i],probs=c(0.025,0.975))
#'    }
#' }
#'
#'
#' # Uncertainty in segment estimation
#' mu.seg.1.unc=res$mcmc$mu1
#' graphics::hist(mu.seg.1.unc,
#'                xlab='obs',
#'                main='Histogram of first segment of observation')
#'
#' uncertainty95.segment <- list()
#' for(i in 1:nS.user){
#'    uncertainty95.segment [[i]] = stats::quantile(res$mcmc[,i],probs=c(0.025,0.975))
#' }
#'
#' # Separate and assign information by identified stable period
#' res$data.p
#' # DIC estimation
#' res$DIC
#'
#' # Setting plot
#'
#' # Transparency
#' alpha <- 125
#'
#' # Set color plot
#' color_customized_rect <- function(alpha){
#'    color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
#'                   rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
#'                   rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
#'                   rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
#'                   rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
#'    return(color)
#' }
#'
#' # Assign period to data
#' obs_segmented <- data.frame()
#'
#' # Conditional to separate non segmentation case
#' if(typeof(res$data.p$obs.p)=='list'){
#'   for(i in 1:length(res$data.p$obs.p)){
#'    obs_segmented_temp=cbind(obs=res$data.p$obs.p[[i]],period=i)
#'    obs_segmented=rbind(obs_segmented,obs_segmented_temp)
#'   }
#' }else{
#'  obs_segmented=data.frame(obs=res$data.p$obs.p,period=1)
#' }
#'
#' # Plot observations
#' plot(x=obs_segmented$obs,
#'      col=factor(obs_segmented$period),
#'      pch=16,
#'      main='Final segmentation',
#'      ylab='obs',
#'      xlab='time')
#'
#' # Plot segments
#' for(i in 1:nS.user){
#'   segments(x0=intervals.time.shift[[i]],
#'            x1=intervals.time.shift[[i+1]],
#'            y0=res$segments[[i]],
#'            y1=res$segments[[i]],
#'            col='blue')
#'   rect(xleft=intervals.time.shift[[i]],
#'        xright=intervals.time.shift[[i+1]],
#'        ybottom=uncertainty95.segment[[i]][1],
#'        ytop=uncertainty95.segment[[i]][2],
#'        col= rgb(0,0,255,max=255,alpha=125,names='blue'),
#'        border = 'transparent')
#' }
#'
#' # Plot shifts
#' if(nS.user!=1){
#'  for(i in 1:(nS.user-1)){
#'   abline(v=res$tau[i],col=color_customized_rect(255)[[i]], lwd=2)
#'   rect(xleft=uncertainty95.shift[[i]][1],
#'        xright=rev(uncertainty95.shift[[i]])[1],
#'        ybottom=min(obs)*2,
#'        ytop=max(obs)*2,
#'        col= color_customized_rect(125)[[i]],
#'        border = 'transparent')
#'   }
#' }
#' @export
#' @importFrom RBaM parameter xtraModelInfo model dataset mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom stats quantile sd
#' @importFrom utils read.table
segmentation.engine <- function(obs,
                                time=1:length(obs),
                                u=0*obs,
                                nS=2,
                                nMin= 1,
                                nCycles=100,
                                burn=0.5,
                                nSlim=max(nCycles/10,1),
                                temp.folder=file.path(tempdir(),'BaM')){


  if(length(obs)<nS){
    stop('Number of observations is lower than the number of segments',call.=FALSE)
  }
  if(any(is.na(obs)) | any(is.na(time)) | any(is.na(u))){
    stop('Missing values not allowed in observation, time and uncertainty')
  }
  if(nS<=0){
    stop('Maximum number of segments should be larger than 0',call.=FALSE)
  }

  # sort data frame case time not ascending
  DF.order <- data.frame(obs=obs,
                         time=time,
                         u=u)
  DF.order.f <- DF.order[order(DF.order$time),]

  obs <- DF.order.f$obs
  time <- DF.order.f$time
  u <- DF.order.f$u
  npar = nS + nS - 1

  priors <- vector(mode = 'list',length = npar)

  for(i in 1:nS){
    priors [[i]] <- RBaM::parameter(name=paste0('mu',i),
                              init=mean(obs),
                              prior.dist = 'FlatPrior' ,
                              prior.par = NULL)
  }

  prior_tau_init <- as.numeric(stats::quantile(time,probs = seq(1,nS-1)/nS))

  if(i>1){
    for(i in 1:(nS-1)){
      priors [[nS+i]] <- RBaM::parameter(name=paste0('tau',i),
                                         init= prior_tau_init[i],
                                         prior.dist = 'FlatPrior' ,
                                         prior.par = NULL)
    }
  }
  # Config_Xtra
  xtra=RBaM::xtraModelInfo(object=c(nS=nS,
                                    tmin_xtra=0,
                                    nmin_xtra=nMin,
                                    option_xtra=1)
  )
  # Model
  mod=RBaM::model(
            ID='Segmentation',
            nX=1,
            nY=1,
            par=priors,
            xtra=xtra)

  # dataset object
  data=RBaM::dataset(X=data.frame(time),
                     Y=data.frame(obs),
                     Yu=data.frame(u),
                     data.dir=temp.folder)


  mcmc_temp=RBaM::mcmcOptions(nCycles=nCycles)

  cook_temp=RBaM::mcmcCooking(burn=burn,
                        nSlim=nSlim)

  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Constant",
                                         par = list(RBaM::parameter(name="gamma1",
                                                              init=stats::sd(obs) ,
                                                              prior.dist = "FlatPrior+"))))

  # Run BaM executable
  RBaM::BaM(mod=mod,
            data=data,
            workspace=temp.folder,
            mcmc=mcmc_temp,
            cook = cook_temp,
            dir.exe = file.path(find.package("RBaM"), "bin"),
            remnant = remnant_prior
  )

  mcmc.segm    <- utils::read.table(file=file.path(temp.folder,"Results_Cooking.txt"),header=TRUE)
  mcmc.DIC     <- utils::read.table(file=file.path(temp.folder,"Results_DIC.txt"),header=FALSE)
  resid.segm   <- utils::read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE)

  # unlink(temp.folder, recursive=TRUE)

  colnames(mcmc.segm)[ncol(mcmc.segm)-1] <- "structural_sd"

  simulation.MAP <- resid.segm$Y1_sim

  data = data.frame(time=time,
                    obs=obs,
                    u=u,
                    I95_lower=obs+stats::qnorm(0.025)*u,
                    I95_upper=obs+stats::qnorm(0.975)*u,
                    period = 1)
  if(nS==1){
    obss=obs # Subseries = whole series
    segments.MAP=simulation.MAP # Sub series = whole series
    times=time
    us=u

    shift=data.frame(tau=numeric(0), # no shift time
                     I95_lower=numeric(0),
                     I95_upper=numeric(0))
    tau.MAP=shift$tau
  } else {

    shift <- c()
    for(j in 1:(nS-1)){

        shift.time.p.unc=data.frame(tau=mcmc.segm[which.max(mcmc.segm$LogPost),
                                                  (nS+j)],
                                    I95_lower=stats::quantile(mcmc.segm[,(nS+j)],
                                                              probs=c(0.025)),
                                    I95_upper=stats::quantile(mcmc.segm[,(nS+j)],
                                                              probs=c(0.975)))
        shift <- rbind(shift,
                       shift.time.p.unc)
    }

    shift <- shift[order(shift$tau),]

    tau.MAP <- shift$tau
    # Store sub series into a list
    obss=segments.MAP=times=us=periods=vector(mode='list',length=nS)
    intervals.time.shift=c(-Inf,tau.MAP,Inf) # intervals defined by time shifts
    # intervals.time.shift=c(time[1],tau.MAP,rev(time)[1]) # intervals defined by time shifts

    # data_update=data

    for(i in 1:nS){
      # position.ti.p <- which((data_update$time-intervals.time.shift[[i]])>=0)[1]
      # position.tf.p <- rev(which((data_update$time-intervals.time.shift[[i+1]])<=0))[1]
      position.ti.p <- which((time-intervals.time.shift[[i]])>=0)[1]
      position.tf.p <- rev(which((time-intervals.time.shift[[i+1]])<0))[1]
      # obss[[i]]=data_update$obs[position.ti.p:position.tf.p]
      obss[[i]]=obs[position.ti.p:position.tf.p]
      segments.MAP[[i]]=simulation.MAP[position.ti.p:position.tf.p]
      times[[i]]=time[position.ti.p:position.tf.p]
      # times[[i]]=data_update$time[position.ti.p:position.tf.p]
      us[[i]]=u[position.ti.p:position.tf.p]
      # us[[i]]=data_update$u[position.ti.p:position.tf.p]
      periods[[i]]=rep(i,length(obss[[i]]))

      # data_update <- data_update[-(position.ti.p:position.tf.p),]
    }
    data$period = unlist(periods)
  }

  return(list(summary = list(data=data,
                             shift=shift),
              tau=tau.MAP,
              segments=segments.MAP,
              mcmc=mcmc.segm,
              data.p = list(obs.p=obss,time.p=times,u.p=us),
              DIC=mcmc.DIC[1,2]))
}

#' Segmentation
#'
#' Segmentation procedure for a \strong{unknown} given number of segments
#'
#' @param obs real vector, observations
#' @param time real vector, time
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#'
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#'   \item nS: integer, optimal number of segments following DIC criterion
#' }
#' @examples
#' # Set random generation
#' set.seed(1)
#'
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Assign a maximum number of user-defined segments to be assessed
#' nSmax.user=3
#'
#' # Run segmentation function
#' res <- segmentation(obs=obs,nSmax=nSmax.user)
#'
#' # Optimal number of segments nSopt
#' nSopt <- res$nS
#' nSopt
#'
#' # Estimated shift time
#' res$results[[nSopt]]$tau
#'
#' # intervals defined by time shifts
#' if(nSopt!=1){
#'  intervals.time.shift=c(res$results[[nSopt]]$data.p$time[[1]][1],
#'                         res$results[[nSopt]]$tau,
#'                         rev(res$results[[nSopt]]$data.p$time[[nSopt]])[1])
#' }else{
#'  intervals.time.shift=list(res$results[[nSopt]]$data.p$time[1],
#'                            rev(res$results[[nSopt]]$data.p$time)[1])
#' }
#'
#' # Maximum a posterior value per segment indexed by the list number
#' res$results[[nSopt]]$segments
#'
#' # Uncertainty in shift time
#' if(nSopt!=1){
#'    Shift=res$results[[nSopt]]$mcmc$tau1
#'    graphics::hist(Shift,
#'                   main='Histogram of first shift')
#'
#'    uncertainty95.shift <- list()
#'    for(i in 1:(nSopt-1)){
#'     uncertainty95.shift[[i]] = stats::quantile(res$results[[nSopt]]$mcmc[,nSopt+i],
#'                                                probs=c(0.025,0.975))
#'    }
#' }
#'
#' # Uncertainty in segment estimation
#' mu.seg.1.unc=res$results[[nSopt]]$mcmc$mu1
#' graphics::hist(mu.seg.1.unc,
#'                xlab='obs',
#'                main='Histogram of first segment of observation')
#'
#' uncertainty95.segment <- list()
#' for(i in 1:nSopt){
#'    uncertainty95.segment [[i]] = stats::quantile(res$results[[nSopt]]$mcmc[,i],
#'                                                  probs=c(0.025,0.975))
#' }
#' # Separate and assign information by identified stable period
#' res$results[[nSopt]]$data.p
#'
#' # DIC estimation
#' res$results[[nSopt]]$DIC
#'
#' # Setting plot
#'
#' # Transparency
#' alpha <- 125
#'
#' # Set color plot
#' color_customized_rect <- function(alpha){
#'    color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
#'                   rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
#'                   rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
#'                   rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
#'                   rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
#'    return(color)
#' }
#'
#' # Assign period to data
#' obs_segmented <- data.frame()
#'
#' # Conditional to separate non segmentation case
#' if(typeof(res$results[[nSopt]]$data.p$obs.p)=='list'){
#'   for(i in 1:length(res$results[[nSopt]]$data.p$obs.p)){
#'    obs_segmented_temp=cbind(obs=res$results[[nSopt]]$data.p$obs.p[[i]],period=i)
#'    obs_segmented=rbind(obs_segmented,obs_segmented_temp)
#'   }
#' }else{
#'  obs_segmented=data.frame(obs=res$results[[nSopt]]$data.p$obs.p,period=1)
#' }
#'
#' # Plot observations
#' plot(x=obs_segmented$obs,
#'      col=factor(obs_segmented$period),
#'      pch=16,
#'      main='Final segmentation',
#'      ylab='obs',
#'      xlab='time')
#'
#' # Plot segments
#' for(i in 1:nSopt){
#'   segments(x0=intervals.time.shift[[i]],
#'            x1=intervals.time.shift[[i+1]],
#'            y0=res$results[[nSopt]]$segments[[i]],
#'            y1=res$results[[nSopt]]$segments[[i]],
#'            col='blue')
#'   rect(xleft=intervals.time.shift[[i]],
#'        xright=intervals.time.shift[[i+1]],
#'        ybottom=uncertainty95.segment[[i]][1],
#'        ytop=uncertainty95.segment[[i]][2],
#'        col= rgb(0,0,255,max=255,alpha=125,names='blue'),
#'        border = 'transparent')
#' }
#'
#' # Plot shifts
#' if(nSopt!=1){
#'  for(i in 1:(nSopt-1)){
#'   abline(v=res$results[[nSopt]]$tau[i],col=color_customized_rect(255)[[i]], lwd=2)
#'   rect(xleft=uncertainty95.shift[[i]][1],
#'        xright=rev(uncertainty95.shift[[i]])[1],
#'        ybottom=min(obs)*2,
#'        ytop=max(obs)*2,
#'        col= color_customized_rect(125)[[i]],
#'        border = 'transparent')
#'   }
#' }
#' @export
segmentation <- function(obs,
                         time=1:length(obs),
                         u=0*obs,
                         nSmax=2,
                         nMin= 1,
                         nCycles=100,
                         burn=0.5,
                         nSlim=max(nCycles/10,1),
                         temp.folder=file.path(tempdir(),'BaM')){


  if(nSmax<=0){
    stop('Maximum number of segments should be larger than 0',call.=FALSE)
  }


  res=vector(mode = 'list',length = nSmax)
  DICs <- rep(NA,nSmax)
  for(i in (1:nSmax)){
    nS <- i
    if(length(obs)<nS){
      DICs [i] <- NA
    }else{
      res[[i]] <- segmentation.engine(obs,time,u,nS,nMin,nCycles,burn,nSlim,temp.folder)
      DICs [i] <- res[[i]]$DIC
    }
  }
  nS=which.min(DICs)
  summary<- res[[nS]]$summary
  return(list(summary=summary,results=res,nS=nS))
}

#' Recursive segmentation
#'
#' Recursive segmentation procedure for a \strong{unknown} given number of segments
#'
#' @param obs real vector, observations
#' @param time real vector, time
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles.
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#'
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#'   \item nS: integer, optimal number of segments following DIC criterion
#'   \item tree : data frame, table for tree structure after segmentation
#' }
#' @examples
#'
#' # Set random generation
#' set.seed(1)
#'
#' # Create series to be segmented
#' obs=c(rnorm(30,mean=0.5,sd=1),
#'       rnorm(30,mean=1.5,sd=1),
#'       rnorm(30,mean=1,sd=1),
#'       rnorm(30,mean=2,sd=1))
#'
#' time=1:length(obs)
#'
#' # Assign a maximum number of user-defined segments to be assessed
#' nSmax.user=3
#'
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs, nSmax=nSmax.user)
#'
#' # Have a look at recursion tree
#' results$tree
#'
#' # Plot tree
#' plotTree(results$tree)
#'
#' # Plot segmentation
#' plotSegmentation(data=results$summary$data,
#'                  shift=results$summary$shift)
#' # Get terminal nodes
#' terminal=which(results$tree$nS==1)
#'
#' # Get node with time shifts
#' nodes.shift.time=which(results$tree$nS!=1)
#'
#' nodes.shift.time
#'
#' # Estimated shift time along with uncertainties
#' shift.time.list <- c()
#' for(i in 1:length(nodes.shift.time)){
#'   nSopt.p = results$res[[nodes.shift.time[[i]]]]$nS
#'   results.p=results$res[[nodes.shift.time[[i]]]]$results
#'   # shift.time.p=results.p[[nSopt.p]]$tau
#'   shift.time.p=cbind(c(results.p[[nSopt.p]]$tau))
#'
#'   for(j in 1:(nSopt.p-1)){
#'
#'     shift.time.p.unc=data.frame(tau=as.numeric(shift.time.p[j,]),
#'                                 u2.5=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
#'                                                      probs=c(0.025)),
#'                                 u97.5=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
#'                                                       probs=c(0.975)))
#'     shift.time.list <- rbind(shift.time.list,
#'                              shift.time.p.unc)
#'   }
#' }
#'
#' all.shift.time <- shift.time.list[order(shift.time.list$tau),]
#' all.shift.time
#'
#' # intervals defined by time shifts
#' intervals.time.shift=c(results$res[[1]]$results[[1]]$data.p$time.p[1],
#'                        all.shift.time$tau,
#'                        rev(results$res[[1]]$results[[1]]$data.p$time.p)[1])
#'
#' # Get stable periods by adding information as period, segment,
#' data.stable <- c()
#' for(i in 1:length(terminal)){
#'   data.stable.p=results$res[[terminal[i]]]$results[[1]]   #Save data from stable period
#'
#'   node = data.frame(obs=data.stable.p$data.p$obs.p,
#'                     time=data.stable.p$data.p$time.p,
#'                     u=data.stable.p$data.p$u.p,
#'                     period = rep(i,length(data.stable.p$data.p$obs.p)))
#'   data.stable = rbind(data.stable,node)
#' }
#'
#' # Setting plot
#'
#' # Transparency
#' alpha <- 125
#'
#' # Set color plot
#' color_customized_rect <- function(alpha){
#'   color <-  list(rgb(0, 255, 170, max = 255, alpha = alpha, names ='green'),
#'                  rgb(0, 221, 255, max = 255, alpha = alpha, names='sky blue'),
#'                  rgb(255, 0, 255, max = 255, alpha = alpha, names='purple'),
#'                  rgb(255, 157, 0, max = 255, alpha = alpha, names='orange'),
#'                  rgb(255, 0, 212, max = 255, alpha = alpha, names='magenta' ))
#'   return(color)
#' }
#'
#' # Plot observations
#' plot(x=data.stable$time,
#'      y=data.stable$obs,
#'      col=data.stable$period,
#'      xlab='time',
#'      ylab='obs',
#'      main='Final segmentation'
#' )
#'
#' # Plot shifts
#' for(i in 1:nrow(all.shift.time)){
#'   abline(v=all.shift.time$tau[i],
#'          col=color_customized_rect(255)[[i]])
#'   rect(xleft=all.shift.time$u2.5[[i]],
#'        xright=all.shift.time$u97.5[[i]],
#'        ybottom=min(obs)*2,
#'        ytop=max(obs)*2,
#'        col= color_customized_rect(50)[[i]],
#'        border = 'transparent')
#' }
#'
#' # Pre-treatment data to extract segments from corresponding node
#' all.global.data <- c()
#' for(i in 1:length(nodes.shift.time)){
#'   nSopt.p = results$res[[nodes.shift.time[[i]]]]$nS
#'   results.p=results$res[[nodes.shift.time[[i]]]]$results
#'
#' relation.parent.children.p=results$tree[which(results$tree$parent==nodes.shift.time[[i]]),]
#'   id.segment.p=which(relation.parent.children.p$nS==1)
#'
#'   for(j in 1:length(id.segment.p)){
#'     segments.p=results.p[[nSopt.p]]$segments[[id.segment.p[[j]]]]
#'     obs.p=results.p[[nSopt.p]]$data.p$obs.p[[id.segment.p[[j]]]]
#'     time.p=results.p[[nSopt.p]]$data.p$time.p[[id.segment.p[[j]]]]
#'     u.p=results.p[[nSopt.p]]$data.p$u.p[[id.segment.p[[j]]]]
#'     unc.segment.p=stats::quantile(results.p[[nSopt.p]]$mcmc[,id.segment.p[[j]]],
#'                                   probs=c(0.025,0.975))
#'
#'
#'     data.p.temp=data.frame(mu=segments.p,
#'                            obs=obs.p,
#'                            u=u.p,
#'                            time=time.p,
#'                            u2.5=rep(unc.segment.p[1],length(segments.p)),
#'                            u97.5=rep(unc.segment.p[2],length(segments.p)))
#'
#'     all.global.data=rbind(all.global.data,data.p.temp)
#'   }
#' }
#'
#' all.global.data=all.global.data[order(all.global.data$time),]
#'
#' final.data.segmented=unname((split(all.global.data,all.global.data$mu)))
#'
#' # Plot segments
#' for(i in 1:(nrow(all.shift.time)+1)){
#'   text(final.data.segmented[[i]]$time,
#'        final.data.segmented[[i]]$obs,
#'        terminal[i],pos=3,
#'        cex=0.8,
#'        col='black')
#'   segments(x0=intervals.time.shift[[i]],
#'            x1=intervals.time.shift[[i+1]],
#'            y0=final.data.segmented[[i]]$mu,
#'            y1=final.data.segmented[[i]]$mu,
#'            col='blue')
#'   rect(xleft=intervals.time.shift[[i]],
#'        xright=intervals.time.shift[[i+1]],
#'        ybottom=final.data.segmented[[i]]$u2.5,
#'        ytop=final.data.segmented[[i]]$u97.5,
#'        col= rgb(0,0,255,max=255,alpha=5,names='blue'),
#'        border = 'transparent')
#' }
#' @export
recursive.segmentation <- function(obs,
                                   time=1:length(obs),
                                   u=0*obs,
                                   nSmax=2,
                                   nMin= 1,
                                   nCycles=100,
                                   burn=0.5,
                                   nSlim=max(nCycles/10,1),
                                   temp.folder=file.path(tempdir(),'BaM')){
  # Initialization
  allRes=list() # store segmentation results for all nodes in a sequential list
  k=0 # Main counter used to control indices in allRes
  tree=data.frame() # store tree structure (parents - children relationship)
  p=1 # Auxiliary counter needed to keep track of children / parents indices
  level=0 # Recursion level. The tree is created level-by-level rather than branch-by-branch
  X=list(obs) # List of all nodes (each corresponding to a subseries of x) to be segmented at this level. Start with a unique node corresponding to the whole series
  TIME=list(time) # List of corresponding times
  U=list(u) # List of corresponding uncertainties
  indices=c(1) # Vector containing the indices of each node - same size as X
  parents=c(0) # Vector containing the indices of the parents of each node - same size as X
  continue=TRUE # Logical determining whether recursion should continue

  while(continue){
    level=level+1 # Increment recursion level
    nX=length(X) # Number of nodes at this level
    keepgoing=rep(NA,nX) # Should recursion continue for each node?
    newX=newTIME=newU=newIndices=newParents=c() # Will be used to update subseries, indices and parents at the end of each recursion level
    m=0 # Local counter used to control indices in the 4 vectors above => reset to 0 at each new level of the recursion
    for(j in 1:nX){ # Loop on each node
      k=k+1 # Increment main counter
      partial.segmentation=segmentation(obs=X[[j]],time=TIME[[j]],u=U[[j]],
                                        nSmax,nMin,nCycles,burn,nSlim,temp.folder) # Apply segmentation to subseries stored in node X[[j]]
      # Save results for this node
      allRes[[k]]=partial.segmentation
      # Save optimal number of segments
      nSopt=partial.segmentation$nS
      # Update recursion tree
      tree=rbind(tree,data.frame(indx=k,level=level,parent=parents[j],nS=nSopt))
      # This was the trickiest part: keeping track of indices and parents
      keepgoing[j]=nSopt>1 # if nS=1, segmentation will not continue for this node which is hence terminal
      if(keepgoing[j]){ # Save results for segmentation at next level
        for(i in 1:nSopt){ # Loop on each segment detected for the current node
          p=p+1 # Increment auxiliary counter
          m=m+1 # Increment local counter
          newX[[m]]=partial.segmentation$results[[nSopt]]$data.p$obs.p[[i]] # Save ith segment (on a total of nS)
          newTIME[[m]]=partial.segmentation$results[[nSopt]]$data.p$time.p[[i]] # Save corresponding times
          newU[[m]]=partial.segmentation$results[[nSopt]]$data.p$u.p[[i]] # Save corresponding uncertainty
          newParents[m]=indices[j] # At next level, the parent of this segment will be the index of current node
          newIndices[m]=p # At next level, the index of this segment will be p
        }
      }
    }
    # Check if recursion should continue at all, i.e. if at least one node is not terminal
    if(all(keepgoing==FALSE)) continue=FALSE
    # Update list of nodes to be further segmented at next level + parents and indices
    X=newX
    TIME=newTIME
    U=newU
    parents=newParents
    indices=newIndices
  }
  # Get terminal nodes
  terminal=which(tree$nS==1)

  # Get stable periods by adding information as period, segment,
  data <- c()
  for(i in 1:length(terminal)){
    data.stable.p=allRes[[terminal[i]]]$results[[1]]   #Save data from stable period

    node = data.frame(time=data.stable.p$data.p$time.p,
                      obs=data.stable.p$data.p$obs.p,
                      u=data.stable.p$data.p$u.p,
                      I95_lower=data.stable.p$data.p$obs.p+stats::qnorm(0.025)*data.stable.p$data.p$u.p,
                      I95_upper=data.stable.p$data.p$obs.p+stats::qnorm(0.975)*data.stable.p$data.p$u.p,
                      id = rep(i,length(data.stable.p$data.p$obs.p)))

    data = rbind(data,node)
  }

  data=data[order(data$time),]
  data$period <- rep(NA,nrow(data))

  period.counter=1
  for(i in 1:length(terminal)){
    data$period[which(unique(data$id)[i]==data$id)] <- period.counter
    period.counter=period.counter+1
  }

  data = data[,-which('id'==colnames(data))]

  # Get node with time shifts
  nodes.shift.time=which(tree$nS!=1)

  shift <- c()
  for(i in 1:length(nodes.shift.time)){
    nSopt.p = allRes[[nodes.shift.time[[i]]]]$nS
    results.p = allRes[[nodes.shift.time[[i]]]]$results

    shift.time.p=cbind(c(results.p[[nSopt.p]]$tau))

    for(j in 1:(nSopt.p-1)){

      shift.time.p.unc=data.frame(tau=as.numeric(shift.time.p[j,]),
                                  I95_lower=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                            probs=c(0.025)),
                                  I95_upper=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                            probs=c(0.975)))
      shift <- rbind(shift,
                     shift.time.p.unc)
    }
  }

  shift <- shift[order(shift$tau),]

  return(list(summary=list(data=data,
                           shift=shift),
              res=allRes,tree=tree))
}





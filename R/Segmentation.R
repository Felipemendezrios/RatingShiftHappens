#' Segmentation engine
#'
#' Segmentation procedure for a \strong{known} given number of segments
#'
#' @param obs real vector, observations
#' @param time vector, time in POSIXct, string or numeric format
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
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'        \item data: data frame, all data with their respective periods after segmentation
#'        \item shift: data frame, all detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item tau: real vector, estimated shift times in numeric or POSIXct format in UTC
#'   \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#' }
#' @examples
#' # Run segmentation engine function at two segments
#' # for data set : RhoneRiver (further details on ?RhoneRiver)
#' res=segmentation.engine(obs=RhoneRiver$H,
#'                         time=RhoneRiver$Year,
#'                         u=RhoneRiver$uH,nS=2)
#' # Data information
#' knitr::kable(head(res$summary$data),
#'              align = 'c',row.names = FALSE)
#' # Shift information
#' knitr::kable(head(res$summary$shift),
#'              align = 'c',row.names = FALSE)
#' # Plot segmentation
#' plotSegmentation(res$summary)
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
  if(trunc(length(obs)/nS)<nMin){
    stop(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                   ') and the number of segments (',nS,')'))
  }
  # sort data frame case time not ascending
  DF.order <- data.frame(obs=obs,
                         time=time,
                         u=u)
  DF.order.f <- DF.order[order(DF.order$time),]

  # Check time format
  numeric.check=TRUE
  # Save class of time from input data
  class.time <- class(time)[1]
  if(!is.numeric(DF.order.f$time)){
    DateTransformed <- DateFormatTransform(date=DF.order.f$time)
    DF.order.f$time <- DateTransformed$time
    origin.date <- DateTransformed$origin
    numeric.check=FALSE
  }

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

  remnantInit=stats::sd(obs)
  if(is.na(remnantInit)){ # happens when nObs=1
    remnantInit=abs(mean(obs))
  }
  if(remnantInit==0){remnantInit=1}
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Constant",
                                         par = list(RBaM::parameter(name="gamma1",
                                                              init=remnantInit,
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

  if(numeric.check!=TRUE){
    time = NumericFormatTransform(numeric.date = time,
                                  class = class.time,
                                  origin.date = origin.date)
  }
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
    rownames(shift) <- NULL
    shift <- shift[order(shift$tau),]

    if(numeric.check!=TRUE &  all(shift$tau!= 0)){
      # Transform all time in POSIXct format
      shift <- data.frame(lapply(shift, function(column) {
        NumericFormatTransform(numeric.date = column,
                               class = class.time,
                               origin.date = origin.date)
      }))

    }
    tau.MAP <- shift$tau

    # Store sub series into a list
    obss=segments.MAP=times=us=periods=vector(mode='list',length=nS)
    intervals.time.shift=c(-Inf,tau.MAP,Inf) # intervals defined by time shifts

    for(i in 1:nS){
      position.ti.p <- which((time-intervals.time.shift[[i]])>=0)[1]
      position.tf.p <- rev(which((time-intervals.time.shift[[i+1]])<0))[1]

      obss[[i]]=obs[position.ti.p:position.tf.p]
      segments.MAP[[i]]=simulation.MAP[position.ti.p:position.tf.p]
      times[[i]]=time[position.ti.p:position.tf.p]
      us[[i]]=u[position.ti.p:position.tf.p]
      periods[[i]]=rep(i,length(obss[[i]]))
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
#' Segmentation procedure for an \strong{unknown} given number of segments
#'
#' @param obs real vector, observations
#' @param time real vector, time in POSIXct, string or numeric format
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
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'       \item data: data frame, all data with their respective periods after segmentation
#'       \item shift: data frame, all detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item res: list, provide all the information of the periods from tree structure
#'   \itemize{
#'       \item tau: real vector, estimated shift times in numeric or POSIXct format in UTC
#'       \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'       \item mcmc: data frame, MCMC simulation
#'       \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'       \item DIC: real, DIC estimation
#'   }
#'   \item nS: integer, optimal number of segments following DIC criterion
#' }
#' @examples
#' # Run segmentation engine function at two segments
#' res=segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH,nSmax=3)
#'
#' # Get lower DIC value and optimal number of segments (to define optimal solution)
#' DIC.df = data.frame(nS=c(1:3),DIC=c(res$results[[1]]$DIC,res$results[[2]]$DIC,res$results[[3]]$DIC))
#' nSopt=res$nS
#'
#' ggplot2::ggplot(DIC.df,ggplot2::aes(x=nS,y=DIC,col=factor(nS)))+
#'   ggplot2::geom_point(size=3,show.legend = FALSE)+
#'   ggplot2::geom_segment(ggplot2::aes(x=nSopt,y=min(DIC)*1.03,xend=nSopt,yend=min(DIC)*1.005),
#'                         arrow=ggplot2::arrow(length=ggplot2::unit(0.5,'cm')),
#'                         color='BLACK',lwd=1, show.legend = FALSE)+
#'   ggplot2::theme_bw()
#' # Data information
#' knitr::kable(head(res$results[[nSopt]]$summary$data),
#'              align = 'c',row.names = FALSE)
#' # Shift information
#' knitr::kable(head(res$results[[nSopt]]$summary$shift),
#'              align = 'c',row.names = FALSE)
#' # Plot segmentation
#' plotSegmentation(res$summary)
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
      warning(paste0('NA was returned because the number of observations (',length(obs),
                     ') is lower than the number of segments (',nS,')'))
      DICs [i] <- NA
    }else if(trunc(length(obs)/nS)<nMin){
      warning(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                     ') and the number of segments (',nS,')'))
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
#' Recursive segmentation procedure for an \strong{unknown} number of segments
#'
#' @param obs real vector, observations
#' @param time real vector, time in POSIXct, string or numeric format
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
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'       \item data: data frame, all data of (H,Q and uQ) with their respective periods after segmentation
#'       \item shift: data frame, all detected shift time
#'    }
#'    \item res: list, provide all the information of the periods from tree structure
#'    \itemize{
#'       \item tau: real vector, estimated shift times
#'       \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'       \item mcmc: data frame, MCMC simulation
#'       \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'       \item DIC: real, DIC estimation
#'       \item nS: integer, optimal number of segments following DIC criterion
#'    }
#'    \item tree: data frame, provide tree structure
#'   }
#' @examples
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH,nSmax=3)
#'
#' # Data information
#' knitr::kable(head(results$summary$data),
#'              align = 'c',row.names = FALSE)
#' # Shift information
#' knitr::kable(head(results$summary$shift),
#'              align = 'c',row.names = FALSE)
#' # Have a look at recursion tree
#' results$tree
#'
#' # Visualize tree structure
#' plotTree(results$tree)
#' # Plot segmentation
#' plotSegmentation(summary=results$summary)
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
  if(any(which(tree$nS!=1))){
    nodes.shift.time=which(tree$nS!=1)

    shift <- c()
    for(i in 1:length(nodes.shift.time)){
      nSopt.p = allRes[[nodes.shift.time[[i]]]]$nS
      results.p = allRes[[nodes.shift.time[[i]]]]$results

      shift.time.p=data.frame(c(results.p[[nSopt.p]]$tau))

      for(j in 1:(nSopt.p-1)){

        shift.time.p.unc=data.frame(tau=shift.time.p[j,],
                                    I95_lower=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                              probs=c(0.025)),
                                    I95_upper=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                              probs=c(0.975)))
        shift <- rbind(shift,
                       shift.time.p.unc)
      }
    }

    rownames(shift) <- NULL
    shift <- shift[order(shift$tau),]

    # Transform uncertainty on the shift in POSIXct format
    if(all(is.numeric(shift$tau)!=TRUE)){
      shift[,c(2,3)] <- data.frame(lapply(shift[,c(2,3)], function(column) {
        NumericFormatTransform(numeric.date = column,
                               class = class(data$time)[1],
                               origin.date = min(data$time))
      }))
    }
  }else{
    shift=NULL
  }

  return(list(summary=list(data=data,
                           shift=shift),
              res=allRes,tree=tree))
}





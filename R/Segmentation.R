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
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment mean value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#' }
#' @examples
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Run segmentation engine function
#' res <- segmentation.engine(obs=obs)
#'
#' # Estimated shift time
#' res$tau
#'
#' # mean value per segment indexed by the list number
#' res$segments
#'
#' # Uncertainty in shift time
#' hist(res$mcmc$tau1)
#'
#' # Separate and assign information by identified stable period
#' res$data.p
#'
#' # DIC estimation
#' res$DIC
#'
#' # Plot
#' plot(obs)
#' lines(x=res$data.p$time.p[[1]],y=res$segments[[1]],col='blue')
#' lines(x=res$data.p$time.p[[2]],y=res$segments[[2]],col='blue')
#' abline(v=res$tau,col='green')
#' abline(v=quantile(res$mcmc$tau1,probs=c(0.025,0.975)),col='green',lty=2)
#' @export
#' @importFrom RBaM parameter xtraModelInfo model dataset mcmcOptions mcmcCooking remnantErrorModel BaM
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
   stop('Number of observations is lower than number of segments',call.=FALSE)
  }
  if(any(is.na(obs)) | any(is.na(time)) | any(is.na(u))){
    stop('Missing values not allowed in observation, time and uncertainty')
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

  prior_tau_init <- as.numeric(quantile(time,probs = seq(1,nS-1)/nS))

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
                                                              init=sd(obs) ,
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

  mcmc.segm    <- read.table(file=file.path(temp.folder,"Results_Cooking.txt"),header=TRUE)
  mcmc.DIC     <- read.table(file=file.path(temp.folder,"Results_DIC.txt"),header=FALSE)
  resid.segm   <- read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE)

  # unlink(temp.folder, recursive=TRUE)

  colnames(mcmc.segm)[ncol(mcmc.segm)-1] <- "structural_sd"

  simulation.mean <- resid.segm$Y1_sim

  if(nS==1){
    obss=obs # Subseries = whole series
    segments=simulation.mean # Subseries = whole series
    times=time
    us=u
    tau.MAP=NULL # no shift time
  } else {
    tau.MAP <- mcmc.segm[which.max(mcmc.segm$LogPost),
                         ((nS+1):(nS+nS-1))]
    # Store subseries into a list
    obss=segments=times=us=vector(mode='list',length=nS)
    intervals.time.shift=c(time[1],tau.MAP,rev(time)[1]) # intervals defined by time shifts

    for(i in 1:nS){

      position.ti.p <- which((time-intervals.time.shift[i])>=0)[1]
      position.tf.p <- rev(which((time-intervals.time.shift[i+1])<=0))[1]

      obss[[i]]=obs[position.ti.p:position.tf.p]
      segments[[i]]=simulation.mean[position.ti.p:position.tf.p]
      times[[i]]=time[position.ti.p:position.tf.p]
      us[[i]]=u[position.ti.p:position.tf.p]
    }
  }

  # Clever way to identify the segments index associated with each observation
  # index <- cumsum(diff(c(Inf,simulation))!=0)

  return(list(tau=tau.MAP,
       segments=segments,
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
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment mean value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#'   \item nS: integer, optimal number of segments following DIC criterion
#' }
#' @examples
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Run segmentation function
#' res <- segmentation(obs=obs)
#'
#' # Optimal number of segments nSopt
#' nSopt <- res$nS
#' nSopt
#'
#' # Estimated shift time
#' res$results[[nSopt]]$tau
#'
#' # Uncertainty in shift time
#' hist(res$results[[nSopt]]$mcmc$tau)
#'
#' # Separate and assign information by identified stable period
#' res$results[[nSopt]]$data.p
#'
#' # DIC estimation
#' res$results[[nSopt]]$DIC
#'
#' # Plot
#' plot(obs)
#' lines(x=res$results[[res$nS]]$data.p$time.p[[1]],y=res$results[[res$nS]]$segments[[1]],col='blue')
#' lines(x=res$results[[res$nS]]$data.p$time.p[[2]],y=res$results[[res$nS]]$segments[[2]],col='blue')
#' abline(v=res$results[[nSopt]]$tau,col='green')
#' abline(v=quantile(res$results[[nSopt]]$mcmc$tau1,probs=c(0.025,0.975)),col='green',lty=2)
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
    res[[i]] <- segmentation.engine(obs,time,u,nS,nMin,nCycles,burn,nSlim,temp.folder)
    DICs [i] <- res[[i]]$DIC
  }

  return(list(results=res,nS=which.min(DICs)))
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
#' @return List with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: list, segment mean value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#'   \item nS: integer, optimal number of segments following DIC criterion
#'   \item tree : data frame, table for tree structure after segmentation
#' }
#' @examples
#' #------------------------------------------------------
#' # First example:
#' #------------------------------------------------------
#'
#' set.seed(2023)
#' # Create series to be segmented
#' obs=c(rnorm(30,mean=0,sd=1),rnorm(30,mean=2,sd=1))
#' time=1:length(obs)
#'
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs)
#'
#' # Have a look at recursion tree
#' results$tree
#'
#' # Get terminal nodes
#' terminal=which(results$tree$nS==1)
#'
#' # Plot original series and terminal nodes defining final segments
#' X11();plot(time,obs)
#' for(i in 1:length(terminal)){
#'data.stable.p=results$res[[terminal[i]]]$results[[1]]   #Save data from stable period
#'node=list(obs=data.stable.p$data.p$obs.p,
#'          times=data.stable.p$data.p$time.p,
#'          u=data.stable.p$data.p$u.p)
#'points(node$times,node$obs,col=i)
#'text(node$times,node$obs,terminal[i],pos=3,col=i)
#'}
#'
#'# Get time shifts
# time.shifts=which(results$tree$nS!=1)
# for(i in 1:length(time.shifts)){
# nSopt.p = results$res[[time.shifts[[i]]]]$nS
# abline(v=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$tau,col='green')
# abline(v=quantile(results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$mcmc$tau1,probs=c(0.025,0.975)),col='green',lty=2)
# segments(
# x0=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$data.p$time.p[[1]][1],
# x1=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$tau,
# y0=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$segments[[1]],
# y1=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$segments[[1]],
# col='blue')
# segments(
# x0=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$tau,
# x1=rev(results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$data.p$time.p[[2]])[1],
# y0=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$segments[[2]],
# y1=results$res[[time.shifts[[i]]]]$results[[nSopt.p]]$segments[[2]],
# col='blue')
#' #'}
#'
#'# Visualize tree with data.tree package
#'if(NROW(results$tree)>1){
#'  tree <- data.tree::as.Node(results$tree[-1,c(3,1)],mode = "network")
#'  plot(tree)
#'} else { # No segmentation took place, make a dummy plot
#'  tree <- data.tree::as.Node(data.frame(1,2),mode = "network")
#'  plot(tree)
#'}
#'
#' #------------------------------------------------------
#' # Second example:
#' #------------------------------------------------------
#' set.seed(2023)
#' obs=c(rnorm(30,mean=0,sd=1),rnorm(30,mean=2,sd=1),rnorm(10,mean=5,sd=1),rnorm(5,mean=1,sd=1))
#' time=1:length(obs)
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs)
#'
#' # Have a look at recursion tree
#' results$tree
#'
#' # Get terminal nodes
#' terminal=which(results$tree$nS==1)
#'
#' # Plot original series and terminal nodes defining final segments
#' X11();plot(time,obs)
#' for(i in 1:length(terminal)){
#' data.stable.p=results$res[[terminal[i]]]$results[[1]]   #Save data from stable period
#' node=list(obs=data.stable.p$data.p$obs.p,
#'          times=data.stable.p$data.p$time.p,
#'          u=data.stable.p$data.p$u.p)
#' points(node$times,node$obs,col=i)
#' text(node$times,node$obs,terminal[i],pos=3,col=i)
#' }
#'
#'# Get time shifts
#' time.shifts=which(results$tree$nS!=1)
#' nSopt.p <- c()
#'
#' for(i in 1:length(time.shifts)){
#' nSopt.p [i] = results$res[[time.shifts[[i]]]]$nS
#'
#' }
#'
#'
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
  continue=TRUE

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
  return(list(res=allRes,tree=tree))
}




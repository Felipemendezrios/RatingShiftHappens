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
#' @param mu_prior list, object describing prior knowledge about residual between the rating curve and observation if user-defined (see details)
#' @param doQuickApprox logical, use quick approximation? see ?Segmentation_quickApprox
#' @param varShift logical, allow for a shifting variance? Only used when doQuickApprox=TRUE.
#' @param alpha real in (0;1), type-I error level of the underlying step-change test. Only used when doQuickApprox=TRUE.
#' @param ... other arguments passed to RBaM::BaM.
#'
#' @return list with the following components:
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'        \item data: data frame, all data with their respective periods after segmentation
#'        \item shift: data frame, all detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item plot : list, data formatted to use as input for some plot functions
#'   \itemize{
#'        \item density.tau: data frame, a table with three columns. The first column indicates the specific shift being analyzed.
#'        The second column contains the values assessed during the MCMC exploration. The last columns shows the probability density associated
#'        with each tested value
#'        \item density.inc.tau: data frame, all information about the 95% credibility interval and the Maximum a posterior (MAP) estimation
#'        for shift times with their associated probability densities
#'   }
#'   \item tau: real vector, estimated shift times in numeric or POSIXct format in UTC
#'   \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'   \item mcmc: data frame, MCMC simulation
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, DIC estimation
#'   \item origin.date.p: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#' }
#' @details
#' The residuals between the gaugings and the RC are defined as follows:
#' \deqn{r_i = \tilde{Q_i} - \hat{Q_i} \quad i = 1, \dots, N}
#' where \eqn{\tilde{Q_i}}  is the gauged discharge, \eqn{\hat{Q_i}} is the RC-estimated discharge and N is the number of gaugings.
#' Non-informative prior knowledge has been provided by default. If the user wants to modify this, it is advised to keep
#' in mind that residual must tend towards zero to obtain satisfactory results. Besides, a same prior will be assigned for all calculation.
#' Additionally, please ensure that the prior has been created using the `parameter` function from `RBaM` package;
#' otherwise, an error message will appear as shown below:
#' - Error in mu_list$name : the $ operator is invalid for atomic vectors
#'
#' More information about prior knowledge on mu parameter available in `sources`.
#'
#' @source \url{https://theses.hal.science/tel-03211343}
#' @source \url{https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2020WR028607}
#'
#' @examples
#' # Run segmentation engine function at two segments
#' # for data set : RhoneRiverAMAX (further details on ?RhoneRiverAMAX)
#' res=Segmentation_Engine(obs=RhoneRiverAMAX$H,
#'                         time=RhoneRiverAMAX$Year,
#'                         u=RhoneRiverAMAX$uH,nS=2,
#'                         mu_prior = list(RBaM::parameter(name=paste0('mu'),
#'                                         init=c(0),
#'                                         prior.dist = 'Gaussian',
#'                                         prior.par = c(0,50))))
#' # Data information
#' knitr::kable(head(res$summary$data),
#'              align = 'c',row.names = FALSE)
#' # Shift information
#' knitr::kable(head(res$summary$shift),
#'              align = 'c',row.names = FALSE)
#' # Plot segmentation
#' PlotSegmentation(res$summary,
#'                  res$plot)
#' @export
#' @import dplyr
#' @importFrom RBaM parameter xtraModelInfo model dataset mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom stats quantile sd approx
#' @importFrom utils read.table
#' @importFrom tidyr gather
Segmentation_Engine <- function(obs,
                                time=1:length(obs),
                                u=0*obs,
                                nS=2,
                                nMin=ifelse(doQuickApprox,3,1),
                                nCycles=100,
                                burn=0.5,
                                nSlim=max(nCycles/10,1),
                                temp.folder=file.path(tempdir(),'BaM'),
                                mu_prior=list(NULL),doQuickApprox=TRUE,
                                varShift=FALSE,alpha=0.1,...){

  if(length(obs)<2)stop('At least 2 observations are required',call.=FALSE)
  if(length(obs)<nS)stop('Number of observations is lower than the number of segments',call.=FALSE)
  if(any(is.na(obs)) | any(is.na(time)) | any(is.na(u)))stop('Missing values not allowed in observation, time and uncertainty')

  if(nS<=0)stop('Maximum number of segments should be larger than 0',call.=FALSE)
  if(trunc(length(obs)/nS)<nMin)stop(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                                            ') and the number of segments (',nS,')'))
  if(is.null(check_vector_lengths(obs,time,u)))stop('The observations, time and uncertainty have not the same length')

  ## Think how to check if mu_prior has been properly provided
  # if(!is.null(mu_prior[[1]])){
  #   if(length(mu_prior)!=1)stop('Please provide a value for `mu_prior`; the same prior is mandatory for calculations')
  #   if(class(mu_prior[[1]]!='parameter'))stop('Please ensure `mu_prior` have been created beforehand using the parameter function from the RBaM package')
  # }

  # sort data frame case time not ascending
  DF.order <- data.frame(obs=obs,
                         time=time,
                         u=u)
  DF.order <- DF.order[order(DF.order$time),]

  # Check time format
  numeric.check=TRUE

  # Get origin date of the segmentation
  origin.date <- min(DF.order$time)

  # Date transformation function to passe to numeric format if necessary
  if(!is.numeric(DF.order$time)){
    DateTransformed <- TransformDateFormat(date=DF.order$time)
    DF.order$time <- DateTransformed$time
    origin.date <- DateTransformed$origin
    numeric.check=FALSE
  }

  obs <- DF.order$obs
  time <- DF.order$time
  u <- DF.order$u

  if(doQuickApprox){
    hasPrior=!sapply(mu_prior,is.null)
    if(hasPrior){
      warning('Quick approximation does not handle prior information. The provided priors will be ignored.')
    }
    nSim=as.integer(100*nCycles*burn/nSlim)
    out=Segmentation_quickApprox(obs=obs,time=time,u=u,nS=nS,nMin=nMin,
                                 nSim=nSim,varShift=varShift,alpha=alpha)
  } else {
    npar = nS + nS - 1

    priors <- vector(mode = 'list',length = npar)
    # Extract mu parameter from if user-defined

    if(is.null(mu_prior[[1]])){
      # Set default for mu_prior if no mu_prior are provided
      for(i in 1:nS){
        priors [[i]] <- RBaM::parameter(name=paste0('mu',i),
                                        init=mean(obs),
                                        prior.dist = 'FlatPrior' ,
                                        prior.par = NULL)
      }
    }else{
      # Use the provided mu parameter from mu_prior
      for(i in 1:nS){
        mu_list <- mu_prior[[1]]
        mu_list$name <- paste0(mu_list$name,i)
        priors [[i]] <- mu_list
      }
    }

    if(i>1){
        # Set default for tau_prior
        prior_tau_init <- as.numeric(stats::quantile(time,probs = seq(1,nS-1)/nS))

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
              remnant = remnant_prior,...)

    mcmc.segm    <- utils::read.table(file=file.path(temp.folder,"Results_Cooking.txt"),header=TRUE)
    mcmc.DIC     <- utils::read.table(file=file.path(temp.folder,"Results_DIC.txt"),header=FALSE)
    resid.segm   <- utils::read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE)

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
      rownames(shift) <- NULL
      shift <- shift[order(shift$tau),]

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

    # Save only MCMC related to shift
    if(nS!=1){

      if(nS==2){
        MCMC.shift.time = data.frame(tau1=mcmc.segm[,c((nS+1):(nS*2-1))])
      }else{
        MCMC.shift.time = mcmc.segm[,c((nS+1):(nS*2-1))]
      }

      MCMC.shift.time.plot <- MCMC.shift.time %>%
        tidyr::gather(key = "Shift", value = "Value")

      # Credibility interval at 95% by shift
      IC_merge=merge(data.frame(MCMC.shift.time.plot%>%
                                  group_by(Shift)%>%
                                  summarize(
                                    Lower_inc = stats::quantile(Value, probs=0.025)
                                  )),
                     data.frame(MCMC.shift.time.plot%>%
                                  group_by(Shift)%>%
                                  summarize(
                                    Upper_inc = stats::quantile(Value,probs= 0.975)
                                  )
                     ))

      # Get density values for each shift time
      density_data<-c()
      density_data.p <- list()
      for(i in 1:(nS-1)){
        density.MCMC.shift.time.p=density(MCMC.shift.time[,i])
        density_data.p[[i]]=data.frame(Shift=IC_merge$Shift[i],
                                       Value=density.MCMC.shift.time.p$x,
                                       Density=density.MCMC.shift.time.p$y)

        density_data=rbind(density_data,density_data.p[[i]])
      }

      density_inc_95 <- c()
      for (i in 1:nrow(IC_merge)) {
        linear.interpolation=stats::approx(x=density_data.p[[i]]$Value,
                                           y=density_data.p[[i]]$Density,
                                           xout=c(as.double(IC_merge[i,-1]),tau.MAP[i]))

        local.res.interpolation <- data.frame(Shift=IC_merge$Shift[i],
                                              tau_lower_inc=linear.interpolation$x[1],
                                              density_tau_lower_inc=linear.interpolation$y[1],
                                              tau_upper_inc=linear.interpolation$x[2],
                                              density_tau_upper_inc=linear.interpolation$y[2],
                                              taU_MAP=linear.interpolation$x[3],
                                              density_taU_MAP=linear.interpolation$y[3])

        density_inc_95 = rbind(density_inc_95,local.res.interpolation)
      }
    }else{
      density_data=NULL
      density_inc_95=NULL
    }
    # Assemble output object
    out=list(summary = list(data=data,
                            shift=shift),
             plot = list(density.tau = density_data,
                         density.inc.tau = density_inc_95),
             tau=tau.MAP,
             segments=segments.MAP,
             mcmc=mcmc.segm,
             data.p = list(obs.p=obss,time.p=times,u.p=us),
             DIC=mcmc.DIC[1,2],
             origin.date.p=origin.date)
  }
  # Transform back all time into original units
  if(!numeric.check){out=transformDatesInOutput(out,origin.date)}
  return(out)
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
#' @param mu_prior list, object describing prior knowledge about residual between the rating curve and observation if user-defined (see details)
#' @param doQuickApprox logical, use quick approximation? see ?Segmentation_quickApprox
#' @param varShift logical, allow for a shifting variance? Only used when doQuickApprox=TRUE.
#' @param alpha real in (0;1), type-I error level of the underlying step-change test. Only used when doQuickApprox=TRUE.
#' @param ... other arguments passed to RBaM::BaM.
#'
#' @return list with the following components:
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'       \item data: data frame, all data with their respective periods after segmentation
#'       \item shift: data frame, all detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item plot : list, data formatted to use as input for some plot functions
#'   \item res: list, provide all the information of the periods from tree structure
#'   \itemize{
#'       \item tau: real vector, estimated shift times in numeric or POSIXct format in UTC
#'       \item segments: list, segment maximum a posterior (MAP) value indexed by the list number
#'       \item mcmc: data frame, MCMC simulation
#'       \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'       \item DIC: real, DIC estimation
#'        \item origin.date: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#'   }
#'   \item nS: integer, optimal number of segments following DIC criterion
#' }
#' @details
#' User may enter prior knowledge about the mu parameter (see `?Segmentation_Engine`) instead of default values.
#' This information must be provided using the parameter function in the RBaM package, as shown in the example.
#'
#' @examples
#' # Run segmentation engine function at two segments
#' res=Segmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#'
#' # Verify optimal nS corresponds to lowest DIC value
#' res$nS
#' res$results[[1]]$DIC
#' res$results[[2]]$DIC
#'
#' # Data information
#' knitr::kable(head(res$results[[res$nS]]$summary$data),
#'              align = 'c',row.names = FALSE)
#' # Shift information
#' knitr::kable(head(res$results[[res$nS]]$summary$shift),
#'              align = 'c',row.names = FALSE)
#' # Plot segmentation
#' PlotSegmentation(res$summary,res$plot)
#' @export
#' @importFrom rlang is_empty
Segmentation <- function(obs,
                         time=1:length(obs),
                         u=0*obs,
                         nSmax=2,
                         nMin=ifelse(doQuickApprox,3,1),
                         nCycles=100,
                         burn=0.5,
                         nSlim=max(nCycles/10,1),
                         temp.folder=file.path(tempdir(),'BaM'),
                         mu_prior = list(),doQuickApprox=TRUE,
                         varShift=FALSE,alpha=0.1,...){

  if(nSmax<=0){
    stop('Maximum number of segments should be larger than 0',call.=FALSE)
  }
  if(length(obs)<2){
    stop('At least 2 observations are required',call.=FALSE)
  }

  res=vector(mode = 'list',length = nSmax)
  DICs <- rep(NA,nSmax)
  for(i in (1:nSmax)){
    nS <- i
    quick=ifelse(nS<3,doQuickApprox,FALSE)
    if(length(obs)<nS){
      warning(paste0('NA was returned because the number of observations (',length(obs),
                     ') is lower than the number of segments (',nS,')'))
      DICs [i] <- NA
    }else if(trunc(length(obs)/nS)<nMin){
      warning(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                     ') and the number of segments (',nS,')'))
      DICs [i] <- NA
    }else{
      # Check if prior knowledge has been provided:
      if(!rlang::is_empty(mu_prior)){
        mu_args <- mu_prior
      }else{
        mu_args <- list(NULL)
      }
      res[[i]] <- Segmentation_Engine(obs,time,u,nS,nMin,nCycles,burn,nSlim,temp.folder,mu_prior=mu_args,
                                      doQuickApprox=quick,varShift=varShift,alpha=alpha,...)
      DICs [i] <- res[[i]]$DIC
    }
  }
  nS=which.min(DICs)
  summary<- res[[nS]]$summary

  return(list(summary=summary,
              plot = list(density.tau = res[[nS]]$plot$density.tau,
                          density.inc.tau = res[[nS]]$plot$density.inc.tau),
              results=res,
              nS=nS,
              origin.date=res[[nS]]$origin.date.p))
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
#' @param mu_prior list, object describing prior knowledge about residual between the rating curve and observation if user-defined (see details)
#' @param doQuickApprox logical, use quick approximation? see ?Segmentation_quickApprox
#' @param varShift logical, allow for a shifting variance? Only used when doQuickApprox=TRUE.
#' @param alpha real in (0;1), type-I error level of the underlying step-change test. Only used when doQuickApprox=TRUE.
#' @param ... other arguments passed to RBaM::BaM.
#'
#' @return list with the following components:
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'       \item data: data frame, all data of (H,Q and uQ) with their respective periods after segmentation
#'       \item shift: data frame, all detected shift time
#'    }
#'    \item plot : list, data formatted to use as input for some plot functions
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
#'    \item origin.date: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#'   }
#' @details
#' User may enter prior knowledge about the mu parameter (see `?Segmentation_Engine`) instead of default values.
#' This information must be provided using the parameter function in the RBaM package, as shown in the example.
#'
#' @examples
#' # Apply recursive segmentation
#' results=Recursive_Segmentation(obs=RhoneRiverAMAX$H,
#'                                time=RhoneRiverAMAX$Year,
#'                                u=RhoneRiverAMAX$uH)
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
#' PlotTree(results$tree)
#' # Plot segmentation
#' PlotSegmentation(summary=results$summary,
#'                  plot_summary=results$plot)
#' @export
Recursive_Segmentation <- function(obs,
                                   time=1:length(obs),
                                   u=0*obs,
                                   nSmax=2,
                                   nMin=ifelse(doQuickApprox,3,1),
                                   nCycles=100,
                                   burn=0.5,
                                   nSlim=max(nCycles/10,1),
                                   temp.folder=file.path(tempdir(),'BaM'),
                                   mu_prior = list(),doQuickApprox=TRUE,
                                   varShift=FALSE,alpha=0.1,...){
  if(length(obs)<2){
    stop('At least 2 observations are required',call.=FALSE)
  }
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
      if(NROW(X[[j]])<nSmax*nMin){ # Can't segment, return default result
        partial.segmentation=getOutputList(time=TIME[[j]],obs=X[[j]],u=U[[j]])
      } else {
        partial.segmentation=Segmentation(obs=X[[j]],time=TIME[[j]],u=U[[j]],
                                          nSmax,nMin,nCycles,burn,nSlim,temp.folder,
                                          mu_prior=mu_prior,doQuickApprox=doQuickApprox,
                                          varShift=varShift,alpha=alpha,...) # Apply segmentation to subseries stored in node X[[j]]
      }
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
    DF.origin.date <-  list()
    density.tau <- c()
    density.inc.tau <- c()
    counter <- 0
    for(i in 1:length(nodes.shift.time)){
      nSopt.p = allRes[[nodes.shift.time[[i]]]]$nS
      results.p = allRes[[nodes.shift.time[[i]]]]$results

      shift.time.p=data.frame(c(results.p[[nSopt.p]]$tau))

      # Add id to the MCMC to estimate density of shift time
      density.tau.p = cbind(allRes[[nodes.shift.time[[i]]]]$plot$density.tau,
                            id_iteration=nodes.shift.time[[i]])
      density.tau = rbind(density.tau,density.tau.p)

      density.inc.tau.p = cbind(allRes[[nodes.shift.time[[i]]]]$plot$density.inc.tau,
                                id_iteration=nodes.shift.time[[i]])
      density.inc.tau = rbind(density.inc.tau,
                              density.inc.tau.p)

      for(j in 1:(nSopt.p-1)){

        shift.time.p.unc=data.frame(tau=shift.time.p[j,],
                                    I95_lower=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                              probs=c(0.025)),
                                    I95_upper=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],
                                                              probs=c(0.975)),
                                    id_iteration=nodes.shift.time[[i]] )
        shift <- rbind(shift,
                       shift.time.p.unc)

        counter=counter+1
        # Get origin date of each segmentation
        DF.origin.date [[counter]] = results.p[[nSopt.p]]$origin.date.p
      }
    }

    rownames(shift) <- NULL


    # Transform uncertainty on the shift in POSIXct format
    if(all(is.numeric(shift$tau)!=TRUE)){
      transformed.IC.shift <- c()
      for(i in 1:length(DF.origin.date)){

        transformed.IC.shift.p <- data.frame(lapply(shift[i,c(2,3)], function(column) {
          NumericFormatTransform(numeric.date = column,
                                 origin.date = DF.origin.date[[i]])
        }))

        transformed.IC.shift=rbind(transformed.IC.shift,transformed.IC.shift.p)
      }
      shift[,c(2,3)]=transformed.IC.shift
    }
    shift <- shift[order(shift$tau),]
  }else{
    shift=density.tau=density.inc.tau=DF.origin.date=NULL
  }

  return(list(summary=list(data=data,
                           shift=shift),
              plot = list(density.tau = density.tau,
                          density.inc.tau = density.inc.tau),
              res=allRes,
              tree=tree,
              origin.date=DF.origin.date))
}

#' Segmentation engine - quick approximation algorithm
#'
#' A quick approximation to the segmentation procedure for either one or two segments, and
#' no prior information. The approximation is based on a likelihood-ratio test for a single step change
#' at an unknown position. DIC0 is hence replaced with the deviance of a no-shift model M0,
#' while DIC1 is replaced by the maximum deviance of a single-shift model M1 + the critical value of the test.
#' This way, a change will be detected when DIC1<DIC0, i.e. when deviance1+critical value<deviance0,
#' i.e. when deviance0-deviance1>critical value, which precisely corresponds to the outcome of the likelihood ratio test.
#' \cr
#' Critical values are computed using the approximations suggested by Gombay and Horvath (1994, 1996a, 1996b, 1997),
#' as detailled in Renard (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' \cr
#' The uncertainty in tau is still estimated by plugging-in point-estimates of mu1, mu2 and sigma
#' for each possible value of tau, and computing the associated likelihood. This likelihood can then be normalized to
#' estimate the posterior pdf of tau, which can then be numerically integrated to give the posterior cdf of tau.
#' Samples of tau values can then be simulated from this cdf using uniform sampling then inverse-cdf transformation.
#' The resulting algorithm is MCMC-free and independent of RBaM.
#'
#' @param obs real vector, observations
#' @param time vector, time in POSIXct, string or numeric format
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nS integer, number of segments, either 1 or 2
#' @param nMin integer, minimum number of observations by segment
#' @param nSim integer, number of simulated tau values
#' @param varShift logical, allow for a shifting variance?
#' @param alpha real in (0;1), type-I error level of the underlying step-change test.
#' @return list with the following components:
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'        \item data: data frame, data augmented with a column denoting the period after segmentation
#'        \item shift: data frame, detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item plot : list, data formatted to use as input for some plot functions
#'   \itemize{
#'        \item density.tau: data frame, a table with three columns. The first column indicates the specific shift being analyzed.
#'        The second column contains the value of the shift time tau. The last columns shows the probability density associated
#'        with each tau value
#'        \item density.inc.tau: data frame, all information about the 95% credibility interval and the Maximum a posterior (MAP) estimation
#'        for shift times with their associated probability densities
#'   }
#'   \item tau: real, estimated shift time in numeric or POSIXct format in UTC
#'   \item segments: list, segment mean value indexed by the list number
#'   \item mcmc: data frame, Monte-Carlo simulations. Note that values for mu's and sigma's are fixed to the values corresponding
#'       to the maxpost tau.
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, pseudo-DIC (see description)
#'   \item origin.date.p: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#' }
#' @examples
#' # Segmentation into two segments for the RhoneRiverAMAX data set (details in ?RhoneRiverAMAX)
#' res=Segmentation_quickApprox(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH,nS=2)
#' res$summary$shift
#' PlotSegmentation(res$summary,res$plot)
#' @export
#' @importFrom stats quantile sd approxfun runif optim
Segmentation_quickApprox <- function(obs,time=1:length(obs),u=0*obs,nS=2,
                                     nMin=3,nSim=500,varShift=FALSE,alpha=0.1){
  n=NROW(obs)
  cll=rep(-Inf,n) # conditional log-likelihoods (for each tau)
  foo=sort.int(time,index.return=TRUE) # make sure everything is sorted chronologically
  obs=obs[foo$ix];time=time[foo$ix];u=u[foo$ix]
  ix=nMin:(n-nMin) # indices of candidate taus
  if(any(diff(ix)<0)){
    mess=paste0('Not enough data (n=',n,') to detect a shift with nMin=',nMin,
                ' values on each side of the shift')
    stop(mess)
  }
  if(nMin<3){stop('nMin cannot be smaller than 3')}
  if(nS>2){stop('Quick-approximation algorithm is only available for nS<3 (i.e. single-shift model at most)')}

  # Start preparing output list
  out=getOutputList_Engine(time,obs,u)
  # no-change model
  start=c(mean(obs),sd(obs))
  if(all(u==0)){ # no need to optimize, estimate is explicit
    theta0=start
    f0=quickApprox_llfunk0(start,obs,u)
  } else { # maximize log-lkh
    opt=optim(start,quickApprox_llfunk0,obs=obs,u=u,control=list(fnscale=-1))
    theta0=opt$par
    f0=opt$value
  }
  if(nS==1){ # stop here and return
    out$segments=rep(theta0[1],n)
    out$mcmc=data.frame(mu1=rep(theta0[1],nSim),
                        structural_sd=rep(theta0[2],nSim),
                        LogPost=rep(f0,nSim))
    out$data.p=list(obs.p=obs,time.p=time,u.p=u)
    out$DIC=-2*f0
    return(out)
  }

  # Single-shift model
  if(varShift){ # mu1,mu2,sig1,sig2
    start=c(theta0[1],theta0[1],theta0[2],theta0[2])
  } else {  # mu1,mu2,sig
    start=c(theta0[1],theta0[1],theta0[2])
  }
  taus=time[ix]
  pars=vector('list',length(time))
  for (i in 1:length(taus)){ # for each tau, find mus/sds by max. likelihood
    foo=quickApprox_getTheta(tau=taus[i],time=time,obs=obs,u=u,varShift=varShift)
    pars[[ix[i]]]=foo$par
    cll[ix[i]]=foo$value
  }
  # Approximate cdf by rectangle integration
  post=exp(cll-max(cll)) # avoid numerical zeros, remains proportional to pdf
  w=c(diff(time),0) # basis of each rectangle
  cdf0=c(0,cumsum(w*post))[1:n] # unnormalized cdf
  cdf=cdf0/max(cdf0) # normalized cdf
  post=post/max(cdf0) # normalized pdf
  # compute pseudo-DIC (see Description)
  d=ifelse(varShift,2,1)
  crit=getCriticalValue(d=d,alpha=alpha,n=n)
  DIC1=-2*max(cll)+crit
  # simulate taus from posterior
  unif=runif(nSim)
  foo=approxfun(x=cdf,y=time,ties='ordered')
  sim=sapply(unif,foo)
  # ready to return
  imax=which.max(post)
  tau=time[imax]
  out$summary$data$period[time>tau]=2
  n1=sum(out$summary$data$period==1)
  n2=sum(out$summary$data$period==2)
  out$summary$shift[1,]=c(tau,quantile(sim,probs=c(0.025,0.975)))
  out$plot$density.tau=data.frame(Shift='tau1',Value=time,Density=post)
  foo=approxfun(x=time,y=post)
  out$plot$density.inc.tau=data.frame(Shift='tau1',
                                      tau_lower_inc=out$summary$shift$I95_lower,
                                      density_tau_lower_inc=foo(out$summary$shift$I95_lower),
                                      tau_upper_inc=out$summary$shift$I95_upper,
                                      density_tau_upper_inc=foo(out$summary$shift$I95_upper),
                                      taU_MAP=tau,
                                      density_taU_MAP=foo(tau))
  out$tau=tau
  out$segments=list(rep(pars[[imax]][1],n1),rep(pars[[imax]][2],n2))
  out$mcmc=data.frame(mu1=pars[[imax]][1],mu2=pars[[imax]][2],tau1=sim)
  if(varShift){
    out$mcmc$sig1=pars[[imax]][3]
    out$mcmc$sig2=pars[[imax]][4]
  } else {
    out$mcmc$structural_sd=pars[[imax]][3]
  }
  foo=approxfun(x=time,y=post)
  out$mcmc$LogPost=foo(sim)
  mask=(time<=tau)
  out$data.p$obs.p=list(obs[mask],obs[!mask])
  out$data.p$time.p=list(time[mask],time[!mask])
  out$data.p$u.p=list(u[mask],u[!mask])
  out$DIC=DIC1
  return(out)
}

#' No-shift likelihood
#'
#' No-shift likelihood function for the quick approximation algorithm.
#'
#' @param theta real vector, parameter vector (mu,sigma)
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a numeric value equal to the log-likelihood of the no-shift model
#' @keywords internal
#' @importFrom stats dnorm
quickApprox_llfunk0 <- function(theta,obs,u){
  out=sum(dnorm(obs,mean=theta[1],sd=sqrt(theta[2]^2+u^2),log=TRUE))
  return(out)
}

#' Single-shift likelihood
#'
#' Single-shift likelihood function used in the quick approximation algorithm.
#'
#' @param theta real vector, parameter vector: either (mu1,mu2,sigma) (treated as a fixed-var model) or
#'     (mu1,mu2,sigma1,sigma2) (treated as a shifting-var model)
#' @param tau real, shift time
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a numeric value equal to the log-likelihood of the single-shift model
#' @keywords internal
#' @importFrom stats dnorm
quickApprox_llfunk1 <- function(theta,tau,time,obs,u){
  mask=(time<=tau)
  if(length(theta)==4){ # shifting var
    s1=theta[3];s2=theta[4]
  } else {
    s1=s2=theta[3]
  }
  l1=sum(dnorm(obs[mask],mean=theta[1],sd=sqrt(s1^2+u[mask]^2),log=TRUE))
  l2=sum(dnorm(obs[!mask],mean=theta[2],sd=sqrt(s2^2+u[!mask]^2),log=TRUE))
  out=l1+l2
  return(out)
}

#' Compute theta given tau
#'
#' Compute theta (mu's and sigma's) given tau (shift time) by maximizing the log-likelihood,
#' or using explicit estimates if available (which is the case when u==0). Also returns the
#' associated log-likelihood.
#'
#' @param tau real, shift time
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param varShift logical, allow for a shifting variance?
#' @return A list with the following components:
#' \enumerate{
#'   \item par: numeric vector, parameter vector: either (mu1,mu2,sigma) (if varShift is FALSE) or
#'     (mu1,mu2,sigma1,sigma2) (if varShift is TRUE)
#'   \item value: real, corresponding log-likelihood
#'  }
#' @keywords internal
#' @importFrom stats optim
quickApprox_getTheta <- function(tau,time,obs,u,varShift){
  mask=(time<=tau)
  start=c()
  start[1]=mean(obs[mask])
  start[2]=mean(obs[!mask])
  if(varShift){
    start[3]=sd(obs[mask])
    start[4]=sd(obs[!mask])
  } else {
    z=obs
    z[mask]=obs[mask]-start[1]
    z[!mask]=obs[!mask]-start[2]
    start[3]=sd(z)
  }
  if(all(u==0)){ # no need to optimize, explicit estimates
    pars=start
    val=quickApprox_llfunk1(theta=start,tau=tau,time=time,obs=obs,u=u)
  } else { # maximize log-lkh
    opt=optim(start,quickApprox_llfunk1,tau=tau,time=time,obs=obs,u=u,control=list(fnscale=-1))
    pars=opt$par
    val=opt$value
  }
  out=list(par=pars,value=val)
  return(out)
}

#' Compute theta given tau
#'
#' Compute theta (mu's and sigma's) given tau (shift time) by maximizing the log-likelihood,
#' or using explicit estimates if available (which is the case when u==0).
#'
#' @inheritParams quickApprox_getTheta
#' @return A numeric vector: either (mu1,mu2,sigma) (if varShift is FALSE) or
#'     (mu1,mu2,sigma1,sigma2) (if varShift is TRUE)
#' @keywords internal
quickApprox_getThetaOnly <- function(tau,time,obs,u,varShift){
  out=quickApprox_getTheta(tau,time,obs,u,varShift)
  return(out$par)
}

#' Default output list
#'
#' Get the default output list returned by Segmentation_Engine when failing.
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a list, see ?Recursive_Segmentation for details.
#' @keywords internal
getOutputList_Engine <- function(time,obs,u){
  out=list()
  out$summary=list(data=data.frame(time=time,obs=obs,u=u,
                                   I95_lower=obs-1.96*u,
                                   I95_upper=obs+1.96*u,
                                   period=1),
                   shift=data.frame(tau=numeric(0),I95_lower=numeric(0),I95_upper=numeric(0)))
  out$plot=list(density.tau=NULL,density.inc.tau=NULL)
  out$tau=numeric(0)
  out$segments=numeric(0)
  out$mcmc=data.frame()
  out$data.p=list()
  out$DIC=NA
  out$origin.date.p=min(time)
  return(out)
}

#' Default output list
#'
#' Get the default output list returned by Segmentation when failing.
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a list, see ?Segmentation for details.
#' @keywords internal
getOutputList <- function(time,obs,u){
  res=getOutputList_Engine(time,obs,u)
  res$data.p=list(obs.p=res$summary$data$obs,time.p=res$summary$data$time,u.p=res$summary$data$u)
  out=list(summary=res$summary,plot=res$plot,results=list(res),nS=1,
           origin.date=res$origin.date.p)
  return(out)
}

#' Date transformer
#'
#' Transform the dates contained in a raw output list back into their original format.
#'
#' @param out list, output from Segmentation (see ?Segmentation for details)
#' @param origin.date date, origin used to transform date back
#' @return a list, see ?Segmentation for details. It is the same as the list `out`,
#'     but with all dates reformatted.
#' @keywords internal
transformDatesInOutput <- function(out,origin.date){
  # Data time (summary)
  out$summary$data$time = NumericFormatTransform(numeric.date = out$summary$data$time,
                                                 origin.date = origin.date)
  # time series indexed by number of segments identified
  if(is.list(out$data.p$time.p)==T){
    for(i in 1:length(out$data.p$time.p)){
      out$data.p$time.p[[i]] =  NumericFormatTransform(numeric.date = out$data.p$time.p[[i]],
                                                       origin.date = origin.date)
    }
  }else{
    out$data.p$time.p = NumericFormatTransform(numeric.date = out$data.p$time.p,
                                               origin.date = origin.date)
  }
  # Rating shift time summary
  if(all(out$summary$shift$tau!= 0)){
    # Transform all time in POSIXct format
    out$summary$shift <- data.frame(lapply(out$summary$shift, function(column) {
      NumericFormatTransform(numeric.date = column,
                             origin.date = origin.date)
    }))
  }
  # Tau in input date format
  out$tau=out$summary$shift$tau
  # Transform density data in POSIXct format if necessary
  if(NROW(out$summary$shift)>0){
    out$plot$density.tau$Value <- NumericFormatTransform(numeric.date = out$plot$density.tau$Value,
                                                         origin.date = origin.date)
    out$plot$density.inc.tau$tau_lower_inc <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$tau_lower_inc,
                                                                     origin.date = origin.date)
    out$plot$density.inc.tau$tau_upper_inc <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$tau_upper_inc,
                                                                     origin.date = origin.date)
    out$plot$density.inc.tau$taU_MAP <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$taU_MAP,
                                                               origin.date = origin.date)
  }
  out$origin.date.p=origin.date
  return(out)
}

#' Compute critical value
#'
#' Compute the critical value of a likelihood-ratio test for a single step change
#' at an unknown position. Based on the approximations suggested by Gombay and Horvath (1994, 1996a, 1996b, 1997),
#' as detailled in Renard (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' @param d integer>1, difference of dimension between M0 (no-change) and M1 (single-change)
#' @param alpha real in (0;1), type-I error level of the test.
#' @param n integer, sample size
#' @return a real number equal to the critical value of the test.
#' @keywords internal
#' @importFrom stats uniroot
getCriticalValue <- function(d,alpha,n){
  logn=log(n)
  if(d==1){
    h=(logn^1.5)/n
    S=log(((1-h)^2)/(h^2))
    foo=uniroot(BrownianFunk,interval=c(2,10),d=d,h=h,S=S,alpha=alpha)
    out=(foo$root)^2
  } else {
    crit=-log(-0.5*log(1-alpha))
    a=sqrt(2*log(logn))
    b=2*log(logn)+0.5*d*log(log(logn))-lgamma(0.5*d)
    out=((crit+b)/a)^2
  }
  return(out)
}

#' Gombay and Horvath's Brownian function
#'
#' Function to be nullified to compute the critical value when d=1.
#' See Gombay and Horvath (1994, 1996a, 1996b, 1997), and Renard
#' (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' @param z real, critical value
#' @param d integer>1, difference of dimension between M0 (no-change) and M1 (single-change)
#' @param h real, h variable (see suggested references)
#' @param S integer, S variable (noted T in suggested references)
#' @param alpha real in (0;1), type-I error level of the test.
#' @return a real number equal to the evaluated Brownian function.
#' @keywords internal
BrownianFunk <- function(z,d,h,S,alpha){
  out=(1/((2^(0.5*d))*gamma(0.5*d)))*(z^d)*exp(-0.5*z^2)*(S-(d/(z^2))*S+(4/(z^2)))-alpha
  return(out)
}

#' Segmentation engine
#'
#' Segmentation procedure for a \strong{known} given number of segments
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nS integer, number of segments
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles.
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#' @return list with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: data frame, segment index and mean value
#'   \item mcmc: data frame, MCMC simulation
#'   \item DIC: real, DIC estimation
#' }
#' @examples
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Run segmentation engine function
#' res <- segmentation_engine(obs=obs)
#'
#' # Estimated shift time
#' res$tau
#'
#' # Uncertainty in shift time
#' hist(res$mcmc$tau1)
#'
#' # DIC estimation
#' res$DIC
#'
#' # Plot
#' plot(obs)
#' lines(res$segments$mean)
#'
#' @export
#' @importFrom RBaM parameter xtraModelInfo model dataset mcmcOptions mcmcCooking remnantErrorModel BaM
segmentation_engine <- function(time=1:length(obs),
                         obs,
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

  npar = nS + nS - 1

  priors <- vector(mode = 'list',length = npar)

  for(j in 1:nS){
    priors [[j]] <- RBaM::parameter(name=paste0('mu',j),
                              init=mean(obs),
                              prior.dist = 'FlatPrior' ,
                              prior.par = NULL)
  }

  prior_tau_init <- as.numeric(quantile(time,probs = seq(1,nS-1)/nS))

  if(j>1){
    for(j in 1:(nS-1)){
      priors [[nS+j]] <- RBaM::parameter(name=paste0('tau',j),
                                         init= prior_tau_init[j],
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

  # mcmcSummary()

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

  if(nS>1){
    tau_MAP <- mcmc.segm[which.max(mcmc.segm$LogPost),
                         ((nS+1):(nS+nS-1))]
  }else{
    tau_MAP=NULL
  }
  simulation <- resid.segm$Y1_sim

  # Clever way to identify the segments index associated with each observation
  index <- cumsum(diff(c(Inf,simulation))!=0)

  return(list(tau=tau_MAP,
       segments=data.frame(index=index,mean=simulation),
       mcmc=mcmc.segm,
       DIC=mcmc.DIC[1,2]))
}

#' Segmentation
#'
#' Segmentation procedure for a \strong{unknown} given number of segments
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles.
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#' @return list with the following components :
#' \enumerate{
#'   \item tau: real vector, estimated shift times
#'   \item segments: data frame, segment index and mean value
#'   \item mcmc: data frame, MCMC simulation
#'   \item DIC: real, DIC estimation
#' }
#' @examples
#' # Create observation vector
#' obs=c(rnorm(25,mean=0,sd=1),rnorm(25,mean=2,sd=1))
#'
#' # Run segmentation function
#' res <- segmentation(obs=obs)
#'
#' # Estimated shift time
#' res$tau
#'
#' # Uncertainty in shift time
#' hist(res$mcmc$tau1)
#'
#' # DIC estimation
#' res$DIC
#'
#' # Plot
#' plot(obs)
#' lines(res$segments$mean)
#'
#' @export
segmentation <- function(time=1:length(obs),
                                obs,
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
    res[[i]] <- segmentation_engine(time,obs,u,nS,nMin,nCycles,burn,nSlim,temp.folder)
    DICs [i] <- res[[i]]$DIC
  }

  return(list(results=res,nS=which.min(DICs)))
}




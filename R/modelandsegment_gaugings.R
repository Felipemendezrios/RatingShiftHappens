#' Recursive modelling  and segmentation
#'
#' Recursive procedure for an \strong{unknown} number of segments, updating the rating curve at each iteration before segmentation
#'
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param time vector, time in POSIXct, string or numeric format
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param nSmax integer, maximum number of segments to assess
#' @param nMin integer, minimum number of observations by segment
#' @param nCycles integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles.
#' @param burn real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlim integer, MCMC slim step
#' @param temp.folder directory, temporary directory to write computations
#' @param funk the function for estimating the rating curve to be applied: see ‘Details’
#' @param ... optional arguments to funk
#' @param mu_prior list, object describing prior knowledge about residual between the rating curve and observation if user-defined (see details)
#'
#' @details
#' Some functions for estimating the rating curve are available in this package.
#' Use `GetCatalog()$models` to discover the supported models. More information in `?GetCatalog()`.
#' User may enter prior knowledge about the mu parameter (see `?segmentation.engine`, ) instead of default values.
#' This information must be provided using the parameter function in the RBaM package, as shown in the example of
#' creating an prior knowledge on "a" for computing the rating curve
#' @return List with the following components :
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'       \item data: data frame, all data of (H,Q and uQ) with their respective periods after segmentation
#'       \item shift: data frame, all detected shift time
#'       \item param.equation: data frame, parameters estimation
#'       }
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
#' }
#' @export
#' @examples
#' # Define fit model for rating curve
#' fit=fitRC_exponential
#' equation=EquationRC_Exponential
#'
#' # Apply recursive model and segmentation
#' results=recursive.ModelAndSegmentation(H=ArdecheRiverGaugings$H,
#'                                        Q=ArdecheRiverGaugings$Q,
#'                                        time=ArdecheRiverGaugings$Date,
#'                                        uQ=ArdecheRiverGaugings$uQ,
#'                                        nSmax=3,nMin=2,funk=fit)
#'
#' # Data information
#' knitr::kable(head(results$summary$data),
#'              align = 'c',row.names = FALSE)
#'
#' # Shift information
#' knitr::kable(head(results$summary$shift),
#'              align = 'c',row.names = FALSE)
#'
#' # Parameters estimation of the rating curve
#' results$summary$param.equation
#'
#' # Have a look at recursion tree
#' results$tree
#'
#' # Terminal nodes
#' terminal = results$tree$indx[which(results$tree$nS==1)]
#' terminal
#'
#' # Visualize tree structure
#' plotTree(results$tree)
#'
#' # See the arguments of the specified fit model for the rating curve
#' args(equation)
#'
#' # See parameters estimates for each rating curve
#' results$summary$param.equation
#'
#' # It must be pass as input data a and b as parameters of the specified fit model in this case
#' a=results$summary$param.equation$a
#' b=results$summary$param.equation$b
#'
#' # Plot the rating curve after segmentation following a regression exponential
#' plotRC_ModelAndSegmentation(summary=results$summary,
#'                             equation=equation,
#'                             a=a,
#'                             b=b)
#'
#' # Plot the rating curves after segmentation with zoom user-defined
#' plotRC_ModelAndSegmentation(summary=results$summary,
#'                             equation = equation,
#'                             autoscale = FALSE,
#'                             Hmin_user = 1,
#'                             Hmax_user = 2,
#'                             H_step_discretization = 0.01,
#'                             a=a,b=b)
#'
#' # Plot the rating curves after segmentation in log scale
#' plotRC_ModelAndSegmentation(summary=results$summary,
#'                             logscale=TRUE,
#'                             equation = equation,
#'                             a=a,
#'                             b=b)
#'
#' # Plot the rating curves after segmentation in log scale with zoom
#' plotRC_ModelAndSegmentation(summary=results$summary,
#'                             a=a,
#'                             b=b,
#'                             logscale=TRUE,
#'                             equation = equation,
#'                             autoscale = FALSE,
#'                             Hmin_user = 0.5,
#'                             Hmax_user = 2,
#'                             H_step_discretization = 0.01)
#'
#' # Plot shift times in stage record
#' plot_H_ModelAndSegmentation(summary=results$summary,
#'                             plot_summary=results$plot)
#'
#' # Plot shift times in discharge observations
#' plot_Q_ModelAndSegmentation(summary=results$summary,
#'                             plot_summary=results$plot)
#'
#' # Plot residual
#' plotResidual_ModelAndSegmentation(summary=results$summary,
#'                                   plot_summary=results$plot)
#'
#'
#' # example with BaRatin method (k-a-c)
#' fit=fitRC_BaRatinKAC
#'
#' # Hydraulic matrix control is also required.
#' # The `control_matrix_builder` was developed to help the user to create this matrix.
#' # In this case, 3 hydraulic control has been set
#'
#' controlMatrix=matrix(c(1,0,0,0,1,1,0,0,1),ncol=3,nrow=3)
#'
#' # Prior information about input data is required.
#' # The `prior_infor_param_builder` function was developed to help the user to create these objects.
#' #' # Prior information for Ardeche River at Meyras
#'
#' a1=RBaM::parameter(name='a1',init=14.17,prior.dist='Gaussian',prior.par=c(14.17,3.65))
#' k1=RBaM::parameter(name='k1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.6,0.5))
#' c1=RBaM::parameter(name='c1',init=1.5,prior.dist='Gaussian',prior.par=c(1.5,0.025))
#' a2=RBaM::parameter(name='a2',init=26.5,prior.dist='Gaussian',prior.par=c(26.5,8.4))
#' k2=RBaM::parameter(name='k2',init=0,prior.dist='Gaussian',prior.par=c(0,0.5))
#' c2=RBaM::parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
#' a3=RBaM::parameter(name='a3',init=31.82,prior.dist='Gaussian',prior.par=c(31.8,10.9))
#' k3=RBaM::parameter(name='k3',init=1.2,prior.dist='Gaussian',prior.par=c(1.2,0.2))
#' c3=RBaM::parameter(name='c3',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
#'
#' a.object=list(a1,a2,a3)
#' k.object=list(k1,k2,k3)
#' c.object=list(c1,c2,c3)
#'
#' resultsBaRatin=recursive.ModelAndSegmentation(H=ArdecheRiverGaugings$H,
#'                                               Q=ArdecheRiverGaugings$Q,
#'                                               time=ArdecheRiverGaugings$Date,
#'                                               uQ=ArdecheRiverGaugings$uQ,
#'                                               nSmax=3,
#'                                               nMin=2,
#'                                               funk=fit,
#'                                               a.object=a.object,
#'                                               k.object=k.object,
#'                                               c.object=c.object,
#'                                               controlMatrix=controlMatrix
#'                                               )
#'
#' # Visualize tree structure
#' plotTree(resultsBaRatin$tree)
#'
#' # Terminal nodes
#' terminal = resultsBaRatin$tree$indx[which(resultsBaRatin$tree$nS==1)]
#' terminal
#'
#' # Plot the rating curves after using BaRatin method. Function specially created to this method
#' plots_prediction=plotRCPrediction(Hgrid=data.frame(seq(-1,2,by=0.01)),
#'                                   autoscale=FALSE,
#'                                   temp.folder=file.path(tempdir(),'BaM'),
#'                                   CalibrationData='CalibrationData.txt',
#'                                   allnodes=FALSE,
#'                                   nodes=terminal)
#'
#' # Plot rating curve of the first terminal node
#' plots_prediction$PlotRCPred[[1]]
#'
#' # Get data of the first terminal node
#' plots_prediction$allData[[1]]
#'
#' # Plot shift times in stage record
#' plot_H_ModelAndSegmentation(summary=resultsBaRatin$summary,
#'                             plot_summary=resultsBaRatin$plot)
#'
#' # Plot shift times in discharge observations
#' plot_Q_ModelAndSegmentation(summary=resultsBaRatin$summary,
#'                             plot_summary=resultsBaRatin$plot)
#'
#' # Plot residual
#' plotResidual_ModelAndSegmentation(summary=resultsBaRatin$summary,
#'                                   plot_summary=resultsBaRatin$plot)
#'
#' # Plot all gaugings after segmentation in a plot H-Q
#' plotGaugingsSegmented(summary=resultsBaRatin$summary)
recursive.ModelAndSegmentation <- function(H,
                                           Q,
                                           time=1:length(H),
                                           uQ=0*Q,
                                           uH=0*H,
                                           nSmax=2,
                                           nMin= 1,
                                           nCycles=100,
                                           burn=0.5,
                                           nSlim=max(nCycles/10,1),
                                           temp.folder=file.path(tempdir(),'BaM'),
                                           funk=fitRC_exponential,
                                           ...,
                                           mu_prior = list()
                                           ){
  # Initialization
  allRes=list() # store segmentation results for all nodes in a sequential list
  k=0 # Main counter used to control indices in allRes
  tree=data.frame() # store tree structure (parents - children relationship)
  p=1 # Auxiliary counter needed to keep track of children / parents indices
  level=0 # Recursion level. The tree is created level-by-level rather than branch-by-branch
  indices=c(1) # Vector containing the indices of each node - same size as X
  parents=c(0) # Vector containing the indices of the parents of each node - same size as X
  continue=TRUE # Logical determining whether recursion should continue

  if(any(Q<0))stop('Dishcarge cannot be negative, verify information')
  if(any(is.na(time)) | any(is.na(H)) | any(is.na(Q)) | any(is.na(uQ) | any(is.na(uH))))stop('Missing values not allowed in time, stage, discharge or uncertainties')

  check <- check_vector_lengths(time,H,Q,uQ,uH)
  if(is.null(check))stop('time, hauteur, discharge or uncertainties do not have the same length')

  DF.order <- data.frame(H=H,
                         time=time,
                         Q=Q,
                         uQ=uQ,
                         uH=uH)

  DF.order <- DF.order[order(DF.order$time),]


  # Save results from first prediction using the grid for plotting rating curve
  if(identical(funk,fitRC_SimplifiedBaRatin)||identical(funk,fitRC_SimplifiedBaRatinWithPrior)||identical(funk,fitRC_BaRatinKAC)||identical(funk,fitRC_BaRatinBAC)){

    residualsData.all <- list(funk(time=DF.order$time,H=DF.order$H,Q=DF.order$Q,uQ=DF.order$uQ,uH=DF.order$uH,
                                   temp.folder.RC=file.path(temp.folder,'RC'),...)) # initialize first residual data to be segmented

    dir.destination=file.path(temp.folder,'it_1/')
    # Ensure the destination directory exists
    if (!dir.exists(dir.destination)) {
      dir.create(dir.destination)
    }else{
      # erase old files
      unlink(file.path(dir.destination,"*"))
    }
    invisible(copy_files_to_folder(dir.source=file.path(temp.folder,'RC'),
                                   dir.destination=dir.destination))
  }else{
    residualsData.all <- list(funk(time=DF.order$time,H=DF.order$H,Q=DF.order$Q,uQ=DF.order$uQ,uH=DF.order$uH,...)) # initialize first residual data to be segmented

  }
  residualsData <- list(cbind(residualsData.all[[1]][[1]],id=seq(1,nrow(residualsData.all[[1]]$ResultsResiduals))))
  param.equation.p <- list(residualsData.all[[1]][[2]])

  if(any(is.na(residualsData[[1]]))){
    stop('There is not enough data to run the segmentation model')
  }
  residuals=list(residualsData[[1]]$Q_res) # List of all nodes (each corresponding to a subseries of residual) to be segmented at this level. Start with a unique node corresponding to the whole series
  TIME=list(residualsData[[1]]$time) # List of corresponding times
  u_residuals=list(sqrt(residualsData[[1]]$uQ_sim^2+residualsData[[1]]$uQ_obs^2)) # List of corresponding uncertainties

  while(continue){
    level=level+1 # Increment recursion level
    n.residuals=length(residuals) # Number of nodes at this level
    keepgoing=rep(NA,n.residuals) # Should recursion continue for each node?
    new_residuals=newTIME=new_u_residuals=newIndices=newParents=c() # Will be used to update subseries, indices and parents at the end of each recursion level
    m=0 # Local counter used to control indices in the 4 vectors above => reset to 0 at each new level of the recursion
    for(j in 1:n.residuals){ # Loop on each node
      k=k+1 # Increment main counter
      # When the residual cannot be computed by any limitation of the rating curve model, hence terminal node.
      if(any(is.na(residuals[[j]]))){
        nSopt=1
        ##
      }else{
        # Apply segmentation to subseries stored in node residuals[[j]]
        partial.segmentation=segmentation(obs=residuals[[j]],
                                          time=TIME[[j]],
                                          u=u_residuals[[j]],
                                          nSmax,nMin,nCycles,burn,nSlim,
                                          temp.folder=temp.folder,
                                          mu_prior=mu_prior)  # Pass mu_prior to segmentation if user-defined
        # Save results for this node
        allRes[[k]]=partial.segmentation
        # Save optimal number of segments
        nSopt=partial.segmentation$nS
      }
      # Update recursion tree
      tree=rbind(tree,data.frame(indx=k,level=level,parent=parents[j],nS=nSopt))
      # This was the trickiest part: keeping track of indices and parents
      keepgoing[j]=nSopt>1 # if nS=1, segmentation will not continue for this node which is hence terminal
      if(keepgoing[j]){ # Save results for segmentation at next level
        for(i in 1:nSopt){ # Loop on each segment detected for the current node
          p=p+1 # Increment auxiliary counter
          m=m+1 # Increment local counter
          data.local.segmented = partial.segmentation$results[[nSopt]]$summary$data # data summary of local node
          id.lines=which(data.local.segmented$period==i) # identify the ids with a stable period
          newTIME[[m]]=data.local.segmented$time[id.lines] # Save corresponding times
          newParents[m]=indices[j] # At next level, the parent of this segment will be the index of current node
          newIndices[m]=p # At next level, the index of this segment will be p
          NewH=residualsData[[newParents[m]]]$H[match(id.lines,residualsData[[newParents[m]]]$id)] # match ids according to id indexed in the table to target a integer value rather than time (format issues)
          NewQ=residualsData[[newParents[m]]]$Q_obs[match(id.lines,residualsData[[newParents[m]]]$id)]
          NewuQ=residualsData[[newParents[m]]]$uQ_obs[match(id.lines,residualsData[[newParents[m]]]$id)]
          NewuH=residualsData[[newParents[m]]]$uH[match(id.lines,residualsData[[newParents[m]]]$id)]

          # Save results from first prediction using the grid for plotting rating curve
          if(identical(funk,fitRC_SimplifiedBaRatin)||identical(funk,fitRC_SimplifiedBaRatinWithPrior)||identical(funk,fitRC_BaRatinKAC)||identical(funk,fitRC_BaRatinBAC)){

            # Update rating curve estimation
            residualsData.all[[p]] <- funk(time=newTIME[[m]],H=NewH,Q=NewQ,uQ=NewuQ,uH=NewuH,
                                           temp.folder.RC=file.path(temp.folder,'RC'),...)

            dir.destination=file.path(temp.folder,paste0('it_',p,'/'))
            # Ensure the destination directory exists
            if (!dir.exists(dir.destination)) {
              dir.create(dir.destination)
            }else{
              # erase old files
              unlink(file.path(dir.destination,"*"))
            }

            invisible(copy_files_to_folder(dir.source=file.path(temp.folder,'RC'),
                                           dir.destination=dir.destination))
          }else{
            # Update rating curve estimation
            residualsData.all[[p]] <- funk(time=newTIME[[m]],H=NewH,Q=NewQ,uQ=NewuQ,uH=NewuH,...)
          }
          residualsData[[p]] <- data.frame(residualsData.all[[p]][[1]],id=seq(1,nrow(residualsData.all[[p]][[1]])))
          param.equation.p[[p]] <- residualsData.all[[p]][[2]]

          if(any(is.na(residualsData[[p]]))){
            new_residuals[[m]]=NA # Save ith segment (on a total of nS)
            new_u_residuals[[m]]=NA # Save corresponding uncertainty
            residualsData[[p]] =data.frame(time=newTIME[[m]],     # Save data of terminal
                                           H=NewH,
                                           uH=NewuQ,
                                           Q_obs=NewQ,
                                           Q_sim=NA,
                                           Q_res=NA,
                                           uQ_obs=NewuQ,
                                           uQ_sim=NA,
                                           id=seq(1,length(newTIME[[m]])))
            param.equation.p[[p]] = NA
          }else{
          # update residual of new rating curve
          new_residuals[[m]]=residualsData[[p]]$Q_res # Save ith segment (on a total of nS)
          new_u_residuals[[m]]=sqrt(residualsData[[p]]$uQ_sim^2+residualsData[[p]]$uQ_obs^2) # Save corresponding uncertainty
          }
        }
      }
    }
    # Check if recursion should continue at all, i.e. if at least one node is not terminal
    if(all(keepgoing==FALSE)) continue=FALSE
    # Update list of nodes to be further segmented at next level + parents and indices
    residuals=new_residuals
    TIME=newTIME
    u_residuals=new_u_residuals
    parents=newParents
    indices=newIndices
  }
  # Get terminal nodes
  terminal=which(tree$nS==1)

  # Get stable periods by adding information about information to be returned
  data <- c()
  param.equation <- list()
  for(i in 1:length(terminal)){
    data.stable.p=residualsData[[terminal[[i]]]] # Save data from stable period
    node = data.frame(time=data.stable.p$time,
                      H=data.stable.p$H,
                      uH=data.stable.p$uH,
                      Q=data.stable.p$Q_obs,
                      uQ=data.stable.p$uQ_obs,
                      Q_I95_lower=data.stable.p$Q_obs+stats::qnorm(0.025)*data.stable.p$uQ_obs,
                      Q_I95_upper=data.stable.p$Q_obs+stats::qnorm(0.975)*data.stable.p$uQ_obs,
                      Qsim=data.stable.p$Q_sim,
                      uQ_sim=data.stable.p$uQ_sim,
                      Qsim_I95_lower=data.stable.p$Q_sim+stats::qnorm(0.025)*data.stable.p$uQ_sim,
                      Qsim_I95_upper=data.stable.p$Q_sim+stats::qnorm(0.975)*data.stable.p$uQ_sim,
                      Qres=data.stable.p$Q_res,
                      id = rep(i,length(data.stable.p$Q_obs)))
    data = rbind(data,node)

    # Get parameters of rating curve
    param.equation [[i]] = param.equation.p[[terminal[[i]]]]
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
                                    id_iteration=nodes.shift.time[[i]])
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

  # Convert a list to a data frame if possible
  param.equation <- List_to_DF(param.equation)

  return(list(summary=list(data=data,
                           shift=shift,
                           param.equation=param.equation),
              plot = list(density.tau = density.tau,
                          density.inc.tau = density.inc.tau),
              res=allRes,
              tree=tree,
              origin.date=DF.origin.date))
}

#' Exponential recession model with four parameters
#'
#' @param raw.data data frame, all stage record
#' @param data.object object, model object output by the RBaM function
#' @param current.rec.model character, equation selected to modelling recessions
#' @param nCyclesrec integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @return list, recession estimation with their uncertainties
#' @export
#' @importFrom RBaM parameter parameter_VAR xtraModelInfo mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom utils read.table
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_extract
Estimation_Recession_BR1 <- function(raw.data,
                                     data.object,
                                     current.rec.model,
                                     temp.folder.Recession,
                                     nCyclesrec=100,
                                     burnrec=0.5,
                                     nSlimrec=max(nCyclesrec/10,1)){

  rec.data=data.object$data

  # Define prior information :
  # Parameters
  alpha=RBaM::parameter(name='alpha',
                        init=max(rec.data$hrec)-min(rec.data$hrec),
                        prior.dist='FlatPrior+')
  lambda=RBaM::parameter(name='lambda',
                         init=1/stats::median(rec.data$time_rec),          # need to be checked : several lambda to get half-life of a recession
                         prior.dist='FlatPrior+')
  c=RBaM::parameter(name='c',
                    init=1,
                    prior.dist='FlatPrior+')

  Ncurves=max(rec.data$indx)
  beta_k=RBaM::parameter_VAR(name='beta_k',
                             index='indx',
                             d=data.object,
                             # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                             init=rep(min(rec.data$hrec),Ncurves), # first guesses
                             prior.dist=rep('FlatPrior',Ncurves),
                             prior.par =rep(list(NULL), Ncurves)) # prior distributions

  # use xtraModelInfo to pass the names of the inputs and the formulas
  xtra=RBaM::xtraModelInfo(object=list(inputs=c('t'),formulas=current.rec.model))
  # model
  mod=RBaM::model(ID='TextFile',nX=1,nY=1,par=list(alpha,lambda,c,beta_k),xtra=xtra)

  # Cooking
  mcmc_temp=RBaM::mcmcOptions(nCycles=nCyclesrec)
  cook_temp=RBaM::mcmcCooking(burn=burnrec,
                              nSlim=nSlimrec)
  # Error model
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))

  RBaM::BaM(mod=mod,
            data=data.object,
            workspace=temp.folder.Recession,
            mcmc=mcmc_temp,
            cook = cook_temp,
            dir.exe = file.path(find.package("RBaM"), "bin"),
            remnant = remnant_prior)

  residuals   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  residuals.indx=data.frame(residuals,uH_obs=rec.data$uHrec,indx=rec.data$indx)
  residuals.indx$weight = NA

  # Mean max value of all recession
  max.time.rec.local =  data.frame(rec.data %>%
                                     dplyr::group_by(indx)%>%
                                     dplyr::summarize(max=max(time_rec)))

  mean.max.val <- mean(max.time.rec.local$max)

  # Weighting for each recession data item to give greater importance to the long duration of the recession
  for(i in 1:Ncurves){
    rec.duration = rec.data$time_rec[which(rec.data$indx==i)]

    min.val <- min(rec.duration)
    max.val <- max(rec.duration)
    if(max.val < mean.max.val){
      max.val = mean.max.val
    }
    # Define the minimum and maximum weights
    min.weight <- 0.1
    max.weight <- 1 # It could be better than fixed values

    # Calculate normalized weights
    weights <- min.weight + (rec.duration - min.val) / (max.val - min.val) * (max.weight - min.weight)

    residuals.indx$weight[which(residuals.indx$indx==i)] <- weights
  }

  # RMSE is in the same units that input data -> meters
  rmse_curve=data.frame(residuals.indx %>%
                          dplyr::group_by(indx) %>%
                          dplyr::summarize(value_rmse=rmse.weighted(Y1_obs,Y1_sim,weight)))

  # criticize recession modeled curve by curve
  rmse_threshold = 0.3 # tolerance threshold (could be better)

  table_comparison = cbind(rmse_curve,rmse_threshold)

  recession_accepted = table_comparison$indx[which(table_comparison$value_rmse <= table_comparison$rmse_threshold)]

  rejected.curves = unique(rec.data$indx)[!unique(rec.data$indx) %in% recession_accepted]

  # Detect curve to reject in summary.MCMC
  if(length(rejected.curves)!=0){
    summary.MCMC=summary.MCMC[,-match(paste0('beta_k_',rejected.curves),colnames(summary.MCMC))]
  }

  # residuals.all.info helps to plot all recession in a same plot (simulated and observed and a status accepted or jerected)
  # Paste rmse calculation
  residuals.all.info.p=merge(residuals.indx,rmse_curve,by='indx')

  residuals.all.info=data.frame(residuals.all.info.p,
                                rmse_threshold=rmse_threshold,
                                status=NA)

  residuals.all.info$status[residuals.all.info$indx %in% rejected.curves] = "Rejected"
  residuals.all.info$status[residuals.all.info$indx %in% recession_accepted] = "Accepted"

  # Focus in asymptotic height (b) estimated :
  indx_MAP=which(c('MaxPost')==rownames(summary.MCMC))
  indx_st_dev=which(c('St.Dev.')==rownames(summary.MCMC))

  table.asymp.height=(summary.MCMC[c(indx_MAP,indx_st_dev),grep(paste0('beta_k_'),colnames(summary.MCMC))])


  b_obs=data.frame(table.asymp.height ['MaxPost',] %>%
                     tidyr::pivot_longer(
                       cols = dplyr::starts_with("beta_k"),   # Pivot all columns
                       names_to = "b",                        # Create 'b' column for each beta variable
                       values_to = "b_estimated"                    # Store the values in 'value'
                     )) %>%
    dplyr::mutate(indx = stringr::str_extract(b, "\\d+$") %>% as.numeric()) %>%
    dplyr::select(-b)

  ub_obs=data.frame(table.asymp.height ['St.Dev.',] %>%
                      tidyr::pivot_longer(
                        cols = dplyr::starts_with("beta_k"),   # Pivot all columns
                        names_to = "b",                        # Create 'b' column for each beta variable
                        values_to = "ub_estimated"                    # Store the values in 'value'
                      )) %>%
    dplyr::mutate(indx = stringr::str_extract(b, "\\d+$") %>% as.numeric()) %>%
    dplyr::select(-b)

  # Get recession data h(t) for plotting
  rec.data.plot.h.dt = data.frame(raw.data, status=NA)
  rec.data.plot.h.dt$status[raw.data$indx %in% rejected.curves] = "Rejected"
  rec.data.plot.h.dt$status[raw.data$indx %in% recession_accepted] = "Accepted"

  # Get data frame to segment
  df.to.segm.p=merge(b_obs,ub_obs)

  last.date.record=data.frame(raw.data %>%
                                group_by(indx) %>%
                                summarize(date=max(daterec)))

  df.to.segm=merge(df.to.segm.p,last.date.record)

  return(list(df.to.segm=df.to.segm,
              rec.data.plot.h.dt=rec.data.plot.h.dt,
              residuals.all.info=residuals.all.info))
}



#' Catalog of recession models available
#'
#' @param matched.model character, model selected by user
#' @param temp.folder.Recession directory, temporal folder for recession computation
#' @param ... optional arguments specific to matched.model defined by user
#'
#' @return list, recession estimation with their uncertainties
#' @export
Estimation_Recession <- function(matched.model,
                                 temp.folder.Recession= file.path(tempdir(),'BaM','Recession'),
                                 ...){
  # Create named variables for model actions
  model_actions <- list(
    'alpha*exp(-lambda*t^c)+beta_k' = Estimation_Recession_BR1(current.rec.model=matched.model,
                                                               temp.folder.Recession=temp.folder.Recession,
                                                               ...)
    # 'alpha_1_k * exp(-lambda_1*t) + alpha_2_k * exp(-lambda_2*t) + beta_k' = Estimation_Recession_M3(current.rec.model=matched.model,
    #                                                                                                  temp.folder.Recession=temp.folder.Recession,
    #                                                                                                  ...)
  )

  action <- model_actions[[matched.model]]
  return(action)
}

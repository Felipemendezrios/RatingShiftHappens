#' Fit rating curve LOESS model
#'
#' LOESS model used to estimate a simple rating curve
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return list with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data.frame
#'   \itemize{
#'        \item NULL, because parameters cannot be estimated
#'        }
#' }
#' @export
#'
#' @examples
#' # Dataset
#' subset = RhoneRiver[1:20,]
#'
#' # LOESS regression to estimate a simple rating curve
#' fit.funk=fitRC_loess(time=subset$Time,H=subset$H,Q=subset$Q,uQ=subset$uQ)
#' fit=fit.funk$ResultsResiduals
#'
#' # plot rating curve with model results
#' plot(fit$H, fit$Q_obs, ylim = range(c(fit$Q_obs - fit$uQ_obs, fit$Q_obs + fit$uQ_obs)),
#'      xlab = "H", ylab = "Q_obs")
#'
#' # Add error bars following uQ_obs
#' arrows(fit$H, fit$Q_obs - fit$uQ_obs, fit$H, fit$Q_obs + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#'
#' # Loess model regression
#' lines(x=fit$H,y=fit$Q_sim, col="blue",lty=3, lwd=2)
#'
#' # plot residuals
#' plot(x=fit$time,y=fit$Q_res, ylim = range(c(fit$Q_res - fit$uQ_obs, fit$Q_res + fit$uQ_obs)),
#'      xlab='H', ylab='Residual')
#' arrows(fit$time, fit$Q_res - fit$uQ_obs, fit$time, fit$Q_res + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#' abline(h=0, col='red')
fitRC_loess<-function(time,H,Q,uQ){
  if(length(time)<=2){ # because second degree polynomial by default (loess function)
    warning('NA was returned because it not possible to perform LOESS regression with only two points.
    A second degree polynomial requires at least three points for prediction')
    return(list(NA,NA))
  }
  data=data.frame(time=time,H=H,Q=Q)

  # LOESS model for estimating rating curve
  mod <- stats::loess(Q ~ H, data = data)

  # Prediction of discharge Q_rc from the water level (H)
  qsim <- stats::predict(mod, newdata = data$H)

  # Residual calculation
  residuals <- data$Q - qsim

  # Residual standard deviation
  residual_sd <- stats::sd(residuals)  # incorrect interpreted as predictive uncertainty

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
  )
  return(list(ResultsResiduals=ResultsResiduals,
              parameters=NA))
}

#' Fit rating curve using a linear regression
#'
#' Linear regression used to estimate a simple rating curve. Formula : \deqn{Q(h)= a \cdot h + b}
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the linear regression expressed \deqn{Q(h)= a \cdot h + b}
#'   \itemize{
#'        \item a : real value, parameter of the linear regression
#'        \item b : real value, parameter of the linear regression
#'        }
#' }
#' @export
#'
#' @examples
#' # Dataset
#' subset = RhoneRiver[1:20,]
#'
#' # Linear regression to estimate a simple rating curve
#' fit.funk=fitRC_LinearRegression(time=subset$Time,H=subset$H,Q=subset$Q,uQ=subset$uQ)
#' fit=fit.funk$ResultsResiduals
#'
#' # Parameters of the linear regression
#' fit.funk$parameters
#'
#' # plot rating curve with model results
#' plot(fit$H, fit$Q_obs, ylim = range(c(fit$Q_obs - fit$uQ_obs, fit$Q_obs + fit$uQ_obs)),
#'      xlab = "H", ylab = "Q_obs")
#'
#' # Add error bars following uQ_obs
#' arrows(fit$H, fit$Q_obs - fit$uQ_obs, fit$H, fit$Q_obs + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#'
#' # Linear regression
#' lines(x=fit$H,y=fit$Q_sim, col="blue",lty=3, lwd=2)
#'
#' # plot residuals
#' plot(x=fit$time,y=fit$Q_res, ylim = range(c(fit$Q_res - fit$uQ_obs, fit$Q_res + fit$uQ_obs)),
#'      xlab='H', ylab='Residual')
#' arrows(fit$time, fit$Q_res - fit$uQ_obs, fit$time, fit$Q_res + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#' abline(h=0, col='red')
fitRC_LinearRegression <- function(time,H,Q,uQ){
  if(length(time)<2){ # because second degree polynomial by default (loess function)
    warning('NA was returned because it not possible to perform linear regression with fewer than two points.')
    return(list(NA,NA))
  }
  data=data.frame(time=time,H=H,Q=Q)

  # Linear regression for estimating rating curve
  mod <- stats::lm(Q ~ H, data = data)

  # Parameters
  b= stats::coef(mod)[1]
  a= stats::coef(mod)[2]

  # Prediction of discharge Q_rc from the water level (H)
  qsim <- stats::predict(mod, newdata = data.frame(data$H))

  # Residual calculation
  residuals <- data$Q - qsim

  # Residual standard deviation
  residual_sd <- stats::sd(residuals)  # incorrect interpreted as predictive uncertainty

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
                              )

  parameters=data.frame(a=a,b=b)

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))
}

#' Fit rating curve using an exponential regression
#'
#' Exponential regression expressed as \deqn{Q(h)=a \cdot \exp(b \cdot h)}
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the exponential regression expressed \deqn{Q(h)=Q0 \cdot \exp(mu \cdot h)}
#'   \itemize{
#'        \item a : real value, parameter of exponential model regression
#'        \item b : real value, parameter of exponential model regression
#'        }
#' }
#' @export
#'
#' @examples
#' # Dataset
#' subset = RhoneRiver[1:20,]
#'
#' # Exponential regression to estimate a simple rating curve
#' fit.all=fitRC_exponential(time=subset$Time,H=subset$H,Q=subset$Q,uQ=subset$uQ)
#'
#' fit=fit.all$ResultsResiduals
#' fit.param=fit.all$parameters
#'
#' # Parameters of the exponential regression
#' fit.param
#'
#' # plot rating curve with model results
#' plot(fit$H, fit$Q_obs, ylim = range(c(fit$Q_obs - fit$uQ_obs, fit$Q_obs + fit$uQ_obs)),
#'      xlab = "H", ylab = "Q_obs")
#'
#' # Add error bars following uQ_obs
#' arrows(fit$H, fit$Q_obs - fit$uQ_obs, fit$H, fit$Q_obs + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#'
#' # Exponential model regression
#' curve(expr=fit.param$a*exp(fit.param$b*x),from=min(fit$H),to=max(fit$H),col='blue',lwd=2,add=TRUE)
#'
#' # plot residuals
#' plot(x=fit$time,y=fit$Q_res, ylim = range(c(fit$Q_res - fit$uQ_obs, fit$Q_res + fit$uQ_obs)),
#'      xlab='H', ylab='Residual')
#' arrows(fit$time, fit$Q_res - fit$uQ_obs, fit$time, fit$Q_res + fit$uQ_obs, angle = 90,
#'        code = 3, length = 0.1)
#' abline(h=0, col='red')
fitRC_exponential <- function(time,H,Q,uQ){
  if(length(time)<=2){
    warning('NA was returned because it not possible to perform exponentiel regression with only two points.
    A second degree polynomial requires at least three points for prediction')
    return(list(NA,NA))
  }
  data=data.frame(time=time,H=H,Q=Q)

  # Exponential model Q(h)=Q0*exp(mu*h)
  # Convert exponential regression model to linear using the natural log of Q as the response variable and H as the predictor variable
  # Linear model : ln(Q(h))=a1+a2*h
  # a1 = ln(Q0) ;  a2 = mu
  mod <- stats::lm(log(Q) ~ H, data = data)

  a1=stats::coef(mod)[1]
  a2=stats::coef(mod)[2]

  # Initialization of parameters
  Q0.init=exp(a1)
  mu.init=a2

  # Return to the exponential model
  initialization = list(a=Q0.init, b=mu.init)

  # General formula : Q(h)=a*exp(b*h)
  model.exp=stats::nls(Q ~ a*exp(b*H), data=data, start = initialization)

  Q0=stats::coef(model.exp)[1]
  mu=stats::coef(model.exp)[2]

  # Prediction of discharge Q_rc from the water level (H)
  qsim <- stats::predict(model.exp, newdata = data.frame(data$H))

  # Residual calculation
  residuals <- data$Q - qsim

  # Residual standard deviation
  residual_sd <- stats::sd(residuals)  # incorrect interpreted as predictive uncertainty

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
                              )
  parameters=data.frame(a=Q0,b=mu)

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))
}

#' Fit rating curve using Baratin simplified
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} without any prior about the parameters. Only one hydraulic control is assigned
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param temp.folder.RC.Plot.Grid directory, temporary directory to write computations of rating curve with prediction following Hgrid
#' @param temp.folder.RC.H.obs directory, temporary directory to write computations of rating curve only using observed stages
#' @param Hgrid real data.frame, grid to computing the rating curve
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : real value, parameter of geometry component
#'        \item b : real value, parameter of offset (thalweg or streambed)
#'        \item c : real value, parameter of exponent
#'        }
#' }
fitRC_SimplifiedBaRatin<- function(time,H,Q,uQ,
                                   temp.folder.RC=file.path(tempdir(),'BaM','RC'),
                                   Hgrid=data.frame(Hgrid=seq(floor(min(H)),ceiling(max(H)),by=0.01))){
  # if(length(time)<2){
  #   warning('NA was returned because it not possible to perform linear regression with fewer than two points.')
  #   return(list(NA,NA))
  # }
  if(!is.data.frame(Hgrid))(stop('Hgrid must be a data frame'))

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)

  # Define the calibration dataset by specifying
  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=temp.folder.RC)

  # Extra information to lunch model
  controlMatrix = rbind(c(1))
  hmax_grid = max(Hgrid)

  # Declare a prior information about each parameter
  b1=RBaM::parameter(name='b1',
                     init=min(data$H),
                     prior.dist = "FlatPrior")
  #initial guess a1
  if(mean(data$H)!=min(data$H)){
    a1.init=mean(data$Q)/(mean(data$H)-min(data$H))^(5/3)
  }else if(median(data$H)!=min(data$H)){
    a1.init=median(data$Q)/(median(data$H)-min(data$H))^(5/3)
  }else{
    a1.init=mean(data$Q)/(mean(data$H)*0.1-min(data$H))^(5/3)
  }

  a1=RBaM::parameter(name='a1',
                     init=a1.init,
                     prior.dist = "FlatPrior+")

  c1=RBaM::parameter(name='c1',
                     init=5/3,
                     prior.dist = 'FlatPrior+')

  priors=list(b1,a1,c1)

  # Config_xtra
  xtra=RBaM::xtraModelInfo(object=list(controlMatrix,
                                       hmax_grid)
                           )

  # Stitch it all together into a model object
  M=RBaM::model(ID='BaRatinBAC',
                nX=1,nY=1, # number of input/output variables
                par=priors, # list of model parameters
                xtra=xtra) # use xtraModelInfo() to pass the control matrix

  # Cooking
  nCycles=100
  mcmc_temp=RBaM::mcmcOptions(nCycles=nCycles)

  cook_temp=RBaM::mcmcCooking(burn=0.5,
                              nSlim=10)
  # Error model
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "FlatPrior+"),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "FlatPrior+"))))
  # Run BaM executable
  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             mcmc=mcmc_temp,
             cook = cook_temp,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             remnant = remnant_prior)

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot

  # First prediction : observed data :
  # Define a 'prediction' object for total predictive uncertainty only for observed stages
  # mcmc.segm    <- utils::read.table(file=file.path(temp.folder.RC,"Results_Cooking.txt"),header=TRUE)
  # mcmc.DIC     <- utils::read.table(file=file.path(temp.folder.RC,"Results_DIC.txt"),header=FALSE)
  resid.segm   <- utils::read.table(file=file.path(temp.folder.RC,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.RC,"Results_Summary.txt"),header=TRUE)
  summary.MCMC.MAP <- summary.MCMC[nrow(summary.MCMC),]

  totalU=RBaM::prediction(X=data['H'], # stage values
                          spagFiles='QRC_TotalU.spag', # file where predictions are saved
                          data.dir=temp.folder.RC, # a copy of data files will be saved here
                          doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                          doStructural=TRUE) # propagate structural uncertainty ?

  # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
  paramU=RBaM::prediction(X=data['H'],
                          spagFiles='QRC_ParamU.spag',
                          data.dir=temp.folder.RC,
                          doParametric=TRUE,
                          doStructural=FALSE)

  # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
  maxpost=RBaM::prediction(X=data['H'],
                           spagFiles='QRC_Maxpost.spag',
                           data.dir=temp.folder.RC,
                           doParametric=FALSE,
                           doStructural=FALSE)

  # # Define a 'prediction' object with no uncertainty - this corresponds to the rating curve computed with prior information
  # priorU=RBaM::prediction(X=Hgrid,
  #                         spagFiles='QRC_Prior.spag',
  #                         data.dir=temp.folder.RC,
  #                         priorNsim = nrow(mcmc.segm),
  #                         doParametric=TRUE,
  #                         doStructural=FALSE)

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=list(totalU,paramU,maxpost), # list of predictions
             # pred=priorU, # list of predictions
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)
  # env_QRC_TotalUReplace=replace_negatives_or_zero_values(data_frame=env_QRC_TotalU,
  #                                                 columns='all',
  #                                                 consider_zero = FALSE,
  #                                                 replace=0)
  # Predictive uncertainty propagation
  # env_QRC_ParamU=utils::read.table(file.path(temp.folder.RC,'QRC_ParamU.env'),header=TRUE)
  # env_QRC_ParamUReplace=replace_negatives_or_zero_values(data_frame=env_QRC_ParamU,
  #                                                 columns='all',
  #                                                 consider_zero = FALSE,
  #                                                 replace=0)

  # MaxPost prediction
  # env_QRC_MaxPost=utils::read.table(file.path(temp.folder.RC,'QRC_Maxpost.spag'))
  # env_QRC_MaxPostReplace=replace_negatives_or_zero_values(data_frame=env_QRC_MaxPost,
  #                                                  columns='all',
  #                                                  consider_zero = FALSE,
  #                                                  replace=0)

  # Discharge simulation
  qsim <- resid.segm$Y1_sim

  # Residual calculation
  residuals <- resid.segm$Y1_res

  # Residual standard deviation
  residual_sd <- env_QRC_TotalU$Stdev  # incorrect interpreted as predictive uncertainty

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
  )


  parameters=data.frame(a=summary.MCMC.MAP$a1,
                        b=summary.MCMC.MAP$b1,
                        c=summary.MCMC.MAP$c1,
                        k=summary.MCMC.MAP$k1,
                        gamma1=summary.MCMC.MAP$Y1_gamma1,
                        gamma2=summary.MCMC.MAP$Y1_gamma2)

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.RC,
                       dir.destination=file.path(temp.folder.RC, 'Residual'))

  # Second prediction : Hgrid :
  remove_files(dir.source = temp.folder.RC ,
               files_to_keep=c('Results_Cooking.txt',
                               'Results_Residuals.txt',
                               'Results_Summary.txt',
                               "Residual"))

  totalU=RBaM::prediction(X=Hgrid, # stage values
                          spagFiles='QRC_TotalU.spag', # file where predictions are saved
                          data.dir=temp.folder.RC, # a copy of data files will be saved here
                          doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                          doStructural=TRUE) # propagate structural uncertainty ?

  # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
  paramU=RBaM::prediction(X=Hgrid,
                          spagFiles='QRC_ParamU.spag',
                          data.dir=temp.folder.RC,
                          doParametric=TRUE,
                          doStructural=FALSE)

  # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
  maxpost=RBaM::prediction(X=Hgrid,
                           spagFiles='QRC_Maxpost.spag',
                           data.dir=temp.folder.RC,
                           doParametric=FALSE,
                           doStructural=FALSE)

  # # Define a 'prediction' object with no uncertainty - this corresponds to the rating curve computed with prior information
  # priorU=RBaM::prediction(X=Hgrid,
  #                         spagFiles='QRC_Prior.spag',
  #                         data.dir=temp.folder.RC,
  #                         priorNsim = nrow(mcmc.segm),
  #                         doParametric=TRUE,
  #                         doStructural=FALSE)

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=list(totalU,paramU,maxpost), # list of predictions
             # pred=priorU, # list of predictions
             doCalib=FALSE,
             doPred=TRUE)

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))

}

#' Fit rating curve using rectangular channel BaRatin model
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^5/3}
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param prior_a list, prior information about parameter a, describing geometry properties
#' @param prior_b list, prior information about parameter b, describing talweg
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : real value, parameter of geometry component
#'        \item b : real value, parameter of offset (thalweg or streambed)
#'        \item c : real value, parameter of exponent
#'        }
#' }
fitRC_SimplifiedBaRatin<- function(time,H,Q,uQ, prior_a, prior_b){
  # if(length(time)<2){
  #   warning('NA was returned because it not possible to perform linear regression with fewer than two points.')
  #   return(list(NA,NA))
  # }

  prior_a=list('a1',26.516,'LogNormal',c(3.2797,0.36))
  prior_b=list('b1',-0.6,'Gaussian',c(-0.58,1.49))

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)


  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=file.path(tempdir(),'BaM'))


  # Rectangular channel BaRatin model Q(h)=a*(h-b)^5/3
  controlMatrix=list(rbind(c(1)),5)

  b1=RBaM::parameter(name=prior_b[[1]],
                     init=prior_b[[2]],
                     prior.dist = prior_b[[3]] ,
                     prior.par = prior_b[[4]])

  a1=RBaM::parameter(name=prior_a[[1]],
                     init=prior_a[[2]],
                     prior.dist = prior_a[[3]] ,
                     prior.par = prior_a[[4]])

  c1=RBaM::parameter(name='c1',
                     init=5/3,
                     prior.dist = 'FIX' ,
                     prior.par = NULL)

  # Stitch it all together into a model object
  M=RBaM::model(ID='BaRatinBAC',
          nX=1,nY=1, # number of input/output variables
          par=list(b1,a1,c1), # list of model parameters
          xtra=RBaM::xtraModelInfo(object=controlMatrix)) # use xtraModelInfo() to pass the control matrix

  RBaM:: BaM(mod=M,
             data=D,
             workspace = file.path(tempdir(),'BaM'),
             dir.exe = file.path(find.package("RBaM"), "bin"))

  # Exponential model Q(h)=Q0*exp(mu*h)
  # Convert exponential regression model to linear using the natural log of Q as the response variable and H as the predictor variable
  # Linear model : ln(Q(h))=a1+a2*h
  # a1 = ln(Q0) ;  a2 = mu
  mod <- stats::lm(log(Q) ~ H, data = data)

  a1=stats::coef(mod)[1]
  a2=stats::coef(mod)[2]

  # Initialization of parameters
  Q0.init=exp(a1)
  mu.init=a2

  # Return to the exponential model
  initialization = list(a=Q0.init, b=mu.init)

  # General formula : Q(h)=a*exp(b*h)
  model.exp=stats::nls(Q ~ a*exp(b*H), data=data, start = initialization)

  Q0=stats::coef(model.exp)[1]
  mu=stats::coef(model.exp)[2]

  # Prediction of discharge Q_rc from the water level (H)
  qsim <- stats::predict(model.exp, newdata = data.frame(data$H))

  # Residual calculation
  residuals <- data$Q - qsim

  # Residual standard deviation
  residual_sd <- stats::sd(residuals)  # incorrect interpreted as predictive uncertainty

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
  )
  parameters=data.frame(a=Q0,b=mu)

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))
}


#' Fit rating curve using simplified BaRatin model
#'
#' Rating curve estimated by using BaRatin approach with the parameter c set to (5/3) and single rectangular weir as hydraulic control. The input data must be already ordered.
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return list with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ (5/3)}
#'   \itemize{
#'        \item a : real value, parameter of geometry component
#'        \item b : real value, parameter of offset (thalweg or streambed)
#'        }
#' }
#'#' fitRC_Simplified_BaRatin=function(time,H,Q,uQ,
#'#'                                   b.distr = 'Gaussian' ,a.distr = 'LogNormal',
#'#'                                   b.prior, st_b.prior,
#'  #'                                  Bc.prior, KS.prior, S0.prior,
#'  #'                                  hmax_Xtra=1.5*max(H)){
#'  #'   # Hardcoded variables
#'  #'   hyr_contr_matrix=matrix(c(1),nrow=1,ncol=1)
#'  #'   c = 5/3
#'  #'
#'  #'   if(is.null(check_square_matrix(hyr_contr_matrix)))stop('Matrix must be square. Try another hydraulic matrix control')
#'  #'   if(is.null(check_vector_lengths(hyr_contr_matrix,
#'  #'                                   b.prior, Bc.prior, KS.prior, S0.prior)))stop('Prior information of the parameters does not have the same length')
#'  #'   if(is.null(check_vector_lengths(hyr_contr_matrix,
#'  #'                                   st_b.prior,st_Bc.prior,st_KS.prior, st_S0.prior)))stop('Prior information of the parameters does not have the same length')
#'  #'
#'  #'   #  How to handle uncertainty of the parameters? because depends on the distribution, number of parameters can vary
#'  #'   if(is.null(check_param_distribution(b.distr,st_b.prior)))stop('The number of the parameters does not match with the specified distribution for st_b.prior')
#'  #'   if(is.null(check_param_distribution(a.distr,st_Bc.prior)))stop('The number of the parameters does not match with the specified distribution for st_Bc.prior')
#'  #'   if(is.null(check_param_distribution(a.distr,st_KS.prior)))stop('The number of the parameters does not match with the specified distribution for st_KS.prior')
#'  #'   if(is.null(check_param_distribution(a.distr,st_S0.prior)))stop('The number of the parameters does not match with the specified distribution for st_S0.prior')
#'  #'   if(is.null(check_param_distribution(g1.distr.type,g1.prior)))stop('The number of the parameters does not match with the specified distribution for g1.prior')
#'  #'   if(is.null(check_param_distribution(g2.distr.type,g2.prior)))stop('The number of the parameters does not match with the specified distribution for g2.prior')
#'  #'
#'  #'   npar = 3*ncol(hyr_contr_matrix) #'  Three parameters by control
#'  #'   priors <- vector(mode = 'list',length = npar)
#'  #'
#'  #'   for(i in 1 :ncol(hyr_contr_matrix)){
#'  #'     b1=RBaM::parameter(name=paste0('b',i),
#'  #'                        init=b.prior[i],
#'  #'                        prior.dist = b.distr ,
#'  #'                        prior.par = st_b.prior)
#'  #'       parameter(name='k1',init=-0.5,prior.dist='Uniform',prior.par=c(-1.5,0))
#'  #'
#'  #'   }
#'  #'   #'  put all parameters in priors
#'  #'   data=data.frame(time=time,H=H,Q=Q)
#'  #'
#'  #'   ## ending
#'  #'   #  residual data frame
#'  #'   ResultsResiduals=data.frame(time=data$time,
#'  #'                               H=data$H,
#'  #'                               Q_obs=data$Q,
#'  #'                               Q_sim=qsim,
#'  #'                               Q_res=residuals,
#'  #'                               uQ_obs=uQ,
#'  #'                               uQ_sim=residual_sd
#'  #'   )
#'  #'   parameters=data.frame(a=Q0,b=mu)
#'  #'
#'  #'   return(list(ResultsResiduals=ResultsResiduals,
#'  #'               parameters=parameters))
#'  #' }
#'  #'
#'  #'
#'  #'
#'


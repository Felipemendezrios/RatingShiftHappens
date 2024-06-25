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

#' Fit rating curve using simplified BaRatin method with non-informative prior information with only a rectangular hydraulic control
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} non-informative prior information on the parameters estimated automatically, using stage record and gauging information. Only a rectangular hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data
#' @param temp.folder.RC directory, temporary directory to write computations of rating curve using observed stages and grid for plotting rating curve
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
#'        \item c : real value, parameter of exponent describing the type of hydraulic control
#'        }
#' }
#' @export
fitRC_SimplifiedBaRatin<- function(time,H,Q,uQ,HmaxGrid,
                                   temp.folder.RC=file.path(tempdir(),'BaM','RC')){

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)

  # Define the calibration dataset by specifying
  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=temp.folder.RC)

  # Extra information to lunch model
  controlMatrix = rbind(c(1))
  hmax_grid = abs(HmaxGrid)*10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

  # Declare a prior information about each parameter
  b1=RBaM::parameter(name='b1',
                     init=min(data$H),
                     prior.dist = "Uniform",
                     prior.par=c(1.5*min(data$H),max(data$H)))
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
                     prior.dist = "Gaussian",
                     prior.par=c(5/3,0.05))

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

  # Save data object and model object
  save(D,file = file.path(temp.folder.RC,'DataObject.RData'))
  save(M,file = file.path(temp.folder.RC,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

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

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=totalU, # list of predictions
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

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


  parameters=data.frame(a=summary.MCMC.MAP[,which(a1$name==colnames(summary.MCMC.MAP))],
                        b=summary.MCMC.MAP[,which(b1$name==colnames(summary.MCMC.MAP))],
                        c=summary.MCMC.MAP[,which(c1$name==colnames(summary.MCMC.MAP))],
                        k=summary.MCMC.MAP$k1,
                        gamma1=summary.MCMC.MAP$Y1_gamma1,
                        gamma2=summary.MCMC.MAP$Y1_gamma2)

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.RC,
                       dir.destination=file.path(temp.folder.RC, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.RC ,
                         files_to_keep=c('Results_Cooking.txt',
                                         'Results_Residuals.txt',
                                         'Results_Summary.txt',
                                         'CalibrationData.txt',
                                         'ModelObject.RData',
                                         'DataObject.RData',
                                         'Residual')))

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))

}

#' Fit rating curve using simplified BaRatin method with prior information with only a hydraulic control
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} with prior information on the parameters. Only one hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data
#' @param temp.folder.RC directory, temporary directory to write computations of rating curve using observed stages and grid for plotting rating curve
#' @param a.object object, created by `prior_infor_param_builder` for describing prior information about the geometry properties
#' @param b.object object, created by `prior_infor_param_builder` for describing prior information about the offset (thalweg or streambed)
#' @param c.object object, created by `prior_infor_param_builder` for describing prior information about the type of hydraulic control
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
#'        \item c : real value, parameter of exponent describing the type of hydraulic control
#'        }
#' }
#' @export
fitRC_SimplifiedBaRatinWithPrior<- function(time,H,Q,uQ,
                                            HmaxGrid,
                                            a.object,
                                            b.object,
                                            c.object,
                                            temp.folder.RC=file.path(tempdir(),'BaM','RC')){
  # if(length(time)<2){
  #   warning('NA was returned because it not possible to perform linear regression with fewer than two points.')
  #   return(list(NA,NA))
  # }

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)

  # Define the calibration dataset by specifying
  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=temp.folder.RC)

  # Extra information to lunch model
  controlMatrix = rbind(c(1))
  hmax_grid = abs(HmaxGrid)*10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

  # Declare a prior information about each parameter
  b1=b.object

  #initial guess a1
  a1=a.object

  c1=c.object

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

  # Save data object and model object
  save(D,file = file.path(temp.folder.RC,'DataObject.RData'))
  save(M,file = file.path(temp.folder.RC,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

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

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=totalU, # list of predictions
             # pred=list(totalU,paramU,maxpost), # list of predictions
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

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

  parameters=data.frame(a=summary.MCMC.MAP[,which(a.object$name==colnames(summary.MCMC.MAP))],
                        b=summary.MCMC.MAP[,which(b.object$name==colnames(summary.MCMC.MAP))],
                        c=summary.MCMC.MAP[,which(c.object$name==colnames(summary.MCMC.MAP))],
                        k=summary.MCMC.MAP$k1,
                        gamma1=summary.MCMC.MAP$Y1_gamma1,
                        gamma2=summary.MCMC.MAP$Y1_gamma2)

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.RC,
                       dir.destination=file.path(temp.folder.RC, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.RC ,
                         files_to_keep=c('Results_Cooking.txt',
                                         'Results_Residuals.txt',
                                         'Results_Summary.txt',
                                         'CalibrationData.txt',
                                         'ModelObject.RData',
                                         'DataObject.RData',
                                         'Residual')))

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))

}

#' Fit rating curve using BaRatin method with prior information multi-hydraulic controls
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} with prior information on the parameters. Multi-hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data
#' @param temp.folder.RC directory, temporary directory to write computations of rating curve using observed stages and grid for plotting rating curve
#' @param a.object list of object, created by `prior_infor_param_builder` for describing prior information about the geometry properties for each hydraulic control
#' @param b.object list of object, created by `prior_infor_param_builder` for describing prior information about the offset (thalweg or streambed) for each hydraulic control
#' @param c.object list of object, created by `prior_infor_param_builder` for describing prior information about the type of hydraulic control for each hydraulic control
#' @param controlMatrix matrix, hydraulic control. The function `control_matrix_builder` was developed to help the user to created this control matrix
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
#'        \item a : list, real value of the parameter of geometry component for each hydraulic control
#'        \item b : list, real value of the parameter of offset (thalweg or streambed) for each hydraulic control
#'        \item c : list, real value of the parameter of exponent describing the type of hydraulic control for each hydraulic control
#'        }
#' }
#' @export
fitRC_BaRatin<- function(time,H,Q,uQ,
                         HmaxGrid,
                         a.object,
                         b.object,
                         c.object,
                         controlMatrix,
                         temp.folder.RC=file.path(tempdir(),'BaM','RC')){

  if(is.null(check_vector_lengths(a.object,b.object,c.object)))stop('It must be specified three times the number of hydraulic controls')
  if(!is.matrix(controlMatrix))stop('Control matrix must be a matrix')
  if(ncol(controlMatrix)!=length(a.object))stop('The number of the columns in the control matrix must be identical to the number of objects describing the parameters')
  if(is.null(check_square_matrix(controlMatrix)))stop('Hydraulic control must be a square matrix')
  if(any(controlMatrix!=0 & controlMatrix!=1))stop('Hydraulic control must be filled by 1 (active) and 0 (inactive) for describing hydraulic controls')

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)

  # Define the calibration dataset by specifying
  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=temp.folder.RC)

  # Extra information to lunch model
  hmax_grid = abs(HmaxGrid)*10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

  ncontrols=ncol(controlMatrix)
  # Initialize an empty list for the final result
  priors <- list()

  # Put the lists in a list to iterate over them
  lists <- list(b.object,a.object, c.object)

  # Get the length of the sublists (since all lists have the same size)
  list_length <- length(a.object)

  # Use a loop to add the elements to the final list
  for(i in 1:ncontrols){
    for (lst in lists) {
      priors <- append(priors, lst[i])
    }
  }

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

  # Save data object and model object
  save(D,file = file.path(temp.folder.RC,'DataObject.RData'))
  save(M,file = file.path(temp.folder.RC,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

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

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=totalU, # list of predictions
             # pred=list(totalU,paramU,maxpost), # list of predictions
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

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


  for(i in 1:ncontrols){
    local.parameters.temp=data.frame(summary.MCMC.MAP[,which(a.object[[i]]$name==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(b.object[[i]]$name==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(c.object[[i]]$name==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(paste0('k',i)==colnames(summary.MCMC.MAP))])

    colnames(local.parameters.temp) <- c(paste0('a',i),paste0('b',i),paste0('c',i),paste0('k',i))
    if(i==1){
      local.parameters <- local.parameters.temp
    }else{
      local.parameters <- cbind(local.parameters,local.parameters.temp)
    }
  }

  parameters=cbind(local.parameters,
                   gamma1=summary.MCMC.MAP$Y1_gamma1,
                   gamma2=summary.MCMC.MAP$Y1_gamma2)

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.RC,
                       dir.destination=file.path(temp.folder.RC, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.RC ,
                         files_to_keep=c('Results_Cooking.txt',
                                         'Results_Residuals.txt',
                                         'Results_Summary.txt',
                                         'CalibrationData.txt',
                                         'ModelObject.RData',
                                         'DataObject.RData',
                                         'Residual')))

  return(list(ResultsResiduals=ResultsResiduals,
              parameters=parameters))

}

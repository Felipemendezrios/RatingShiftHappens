#' Fit rating curve LOESS model
#'
#' LOESS model used to estimate a simple rating curve
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return list with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
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
#' @importFrom stats loess predict sd
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
fitRC_loess<-function(time,H,Q,uQ,uH){
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
                              uH=uH,
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
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
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
#' @importFrom stats lm predict sd coef
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
fitRC_LinearRegression <- function(time,H,Q,uQ,uH){
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
                              uH=uH,
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
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
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
#' @importFrom stats lm nls sd coef predict
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
fitRC_exponential <- function(time,H,Q,uQ,uH){
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
                              uH=uH,
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

#' Fit rating curve using simplified BaRatin method (b-a-c) with non-informative prior information with only a rectangular hydraulic control
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} non-informative prior information on the parameters estimated automatically, using stage record and gauging information. Only a rectangular hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data in the historical record
#' @param temp.folder.RC directory, temporary directory to write computations of rating curve using observed stages and grid for plotting rating curve
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : real value, parameter of the geometry component
#'        \item b : real value, parameter of the offset (thalweg or streambed)
#'        \item c : real value, parameter of the exponent describing the type of hydraulic control
#'        \item k : real value, parameter of the activation stage for each hydraulic control
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @importFrom RBaM dataset parameter xtraModelInfo model mcmcOptions mcmcCooking remnantErrorModel BaM prediction
#' @importFrom stats median
#' @importFrom utils read.table
#' @export
fitRC_SimplifiedBaRatin<- function(time,H,Q,uQ,uH,HmaxGrid,
                                   temp.folder.RC=file.path(tempdir(),'BaM','RC')){

  data=data.frame(time=time,H=H,Q=Q,uQ=uQ)

  # Define the calibration dataset by specifying
  D=RBaM::dataset(X=data['H'],
                  Y=data['Q'],
                  Yu=data['uQ'],
                  data.dir=temp.folder.RC)

  # Extra information to lunch model
  controlMatrix = rbind(c(1))
  hmax_grid = abs(HmaxGrid)+10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

  # Declare a prior information about each parameter
  b1=RBaM::parameter(name='b1',
                     init=min(data$H),
                     prior.dist = "Uniform",
                     prior.par=c(1.5*min(data$H),max(data$H)))
  #initial guess a1
  if(mean(data$H)!=min(data$H)){
    a1.init=mean(data$Q)/(mean(data$H)-min(data$H))^(5/3)
  }else if(stats::median(data$H)!=min(data$H)){
    a1.init=stats::median(data$Q)/(stats::median(data$H)-min(data$H))^(5/3)
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
  remnant_prior <-list(RBaM::remnantErrorModel(funk = "Linear",
                                               par = list(RBaM::parameter(name="gamma1",
                                                                          init=1,
                                                                          prior.dist = "Uniform",
                                                                          prior.par = c(0,1000)),
                                                          RBaM::parameter(name="gamma2",
                                                                          init=0.1,
                                                                          prior.dist = "Uniform",
                                                                          prior.par = c(0,1000)))))
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
             remnant = remnant_prior,
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

  # Discharge simulation
  qsim <- resid.segm$Y1_sim

  # Residual calculation
  residuals <- resid.segm$Y1_res

  # Residual standard deviation
  residual_sd <- env_QRC_TotalU$Stdev

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              uH=uH,
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

#' Fit rating curve using simplified BaRatin method (b-a-c)  with prior information with only a hydraulic control
#'
#' Rating curve estimated by the equation \deqn{Q(h)=a \cdot (h-b)^c} with prior information on the parameters. Only one hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data in the historical record
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
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : real value, parameter of the geometry component
#'        \item b : real value, parameter of the offset (thalweg or streambed)
#'        \item c : real value, parameter of the exponent describing the type of hydraulic control
#'        \item k : real value, parameter of the activation stage for each hydraulic control
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRC_SimplifiedBaRatinWithPrior<- function(time,H,Q,uQ,uH,
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
  hmax_grid = abs(HmaxGrid)+10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

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
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))
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
             remnant = remnant_prior,
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

  # Discharge simulation
  qsim <- resid.segm$Y1_sim

  # Residual calculation
  residuals <- resid.segm$Y1_res

  # Residual standard deviation
  residual_sd <- env_QRC_TotalU$Stdev

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              uH=uH,
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

#' Fit rating curve using BaRatin method (b-a-c)  with prior information multi-hydraulic controls
#'
#' Rating curve estimated by the equation \deqn{Q(h) = a \cdot (h-b)^{c} \quad \text{for } (h>k) \quad (\text{and } Q=0 \quad \text{if } h \leq b)}
#' with prior information on the parameters. Multi-hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param HmaxGrid real value, maximum stage of all data in the historical record
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
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : list, real value of the parameter of the geometry component for each hydraulic control
#'        \item b : list, real value of the parameter of the offset (thalweg or streambed) for each hydraulic control
#'        \item c : list, real value of the parameter of the exponent describing the type of hydraulic control for each hydraulic control
#'        \item k : list, real value of the parameter of the activation stage for each hydraulic control
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRC_BaRatinBAC<- function(time,H,Q,uQ,uH,
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
  hmax_grid = abs(HmaxGrid)+10 # 10 -> arbitrary factor value to avoid some bugs of unfeasible starting point due to this parameters

  ncontrols=ncol(controlMatrix)
  # Initialize an empty list for the final result
  priors <- list()

  # Put the lists in a list to iterate over them
  lists <- list(b.object,a.object, c.object)

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
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))
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

  if(all(uH==0)){
    totalU=RBaM::prediction(X=data['H'], # stage values
                            spagFiles='QRC_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.RC, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE) # propagate structural uncertainty ?

  }else{ # Consider uncertainty on stage for vegetation influences (David Besson, CVL )

    set.seed(08071998)
    # Create 100 random values for ht, with the uncertainty given by the user
    htrep=matrix(rnorm(n=length(uH)*100,mean=H,sd=uH),nrow=length(uH),ncol=100)

    totalU=RBaM::prediction(X=list(htrep), # spaghetti of stage values
                            spagFiles='QRC_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.RC, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE) # propagate structural uncertainty ?
  }

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=totalU, # list of predictions
             remnant = remnant_prior,
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

  # Discharge simulation
  qsim <- resid.segm$Y1_sim

  # Residual calculation
  residuals <- resid.segm$Y1_res

  # Residual standard deviation
  residual_sd <- env_QRC_TotalU$Stdev

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              uH=uH,
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

#' Fit rating curve using BaRatin method (k-a-c) with prior information multi-hydraulic controls
#'
#' Rating curve estimated by the equation
#' \deqn{Q(h) = a \cdot (h-b)^{c} \quad \text{for } (h>k) \quad (\text{and } Q=0 \quad \text{if } h \leq b)}
#' with prior information on the parameters. Multi-hydraulic control is considered
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param uH real vector, uncertainty in stage record in meters (as a standard deviation)
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param temp.folder.RC directory, temporary directory to write computations of rating curve using observed stages and grid for plotting rating curve
#' @param a.object list of object, created by `prior_infor_param_builder` for describing prior information about the geometry properties for each hydraulic control
#' @param k.object list of object, created by `prior_infor_param_builder` for describing prior information about the activation stage for each hydraulic control; when the water level falls below , the control becomes inactive;
#' @param c.object list of object, created by `prior_infor_param_builder` for describing prior information about the type of hydraulic control for each hydraulic control
#' @param controlMatrix matrix, hydraulic control. The function `control_matrix_builder` was developed to help the user to created this control matrix
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : data frame, results after fitting curve
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage
#'        \item uH: real value, uncertainty in stage observed (as a standard deviation)
#'        \item Q_obs: real value, discharge observed
#'        \item Q_sim: real value, discharge simulated
#'        \item Q_res: real value, residual between discharge observed and simulated
#'        \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'        \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the simplified BaRatin model : \deqn{Q(h)=a \cdot (H - b)^ c}
#'   \itemize{
#'        \item a : list, real value of the parameter of the geometry component for each hydraulic control
#'        \item b : list, real value of the parameter of the offset (thalweg or streambed) for each hydraulic control
#'        \item c : list, real value of the parameter of the exponent describing the type of hydraulic control for each hydraulic control
#'        \item k : list, real value of the parameter of the activation stage for each hydraulic control
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRC_BaRatinKAC<- function(time,H,Q,uQ,uH,
                            a.object,
                            k.object,
                            c.object,
                            controlMatrix,
                            temp.folder.RC=file.path(tempdir(),'BaM','RC')){

  if(is.null(check_vector_lengths(a.object,k.object,c.object)))stop('It must be specified three times the number of hydraulic controls')
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

  ncontrols=ncol(controlMatrix)
  # Initialize an empty list for the final result
  priors <- list()

  # Put the lists in a list to iterate over them
  lists <- list(k.object,a.object, c.object)

  # Use a loop to add the elements to the final list
  for(i in 1:ncontrols){
    for (lst in lists) {
      priors <- append(priors, lst[i])
    }
  }

  # Config_xtra
  xtra=RBaM::xtraModelInfo(object=controlMatrix)

  # Stitch it all together into a model object
  M=RBaM::model(ID='BaRatin',
                nX=1,nY=1, # number of input/output variables
                par=priors, # list of model parameters
                xtra=xtra) # use xtraModelInfo() to pass the control matrix

  # Cooking
  nCycles=100
  mcmc_temp=RBaM::mcmcOptions(nCycles=nCycles)


  cook_temp=RBaM::mcmcCooking(burn=0.5,
                              nSlim=10)

  # Error model fixed for all functions in the package
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))
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

  if(all(uH==0)){
    totalU=RBaM::prediction(X=data['H'], # stage values
                            spagFiles='QRC_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.RC, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE) # propagate structural uncertainty ?

  }else{ # Consider uncertainty on stage for vegetation influences (David Besson, CVL )

    set.seed(08071998)
    # Create 100 random values for ht, with the uncertainty given by the user
    htrep=matrix(rnorm(n=length(uH)*100,mean=H,sd=uH),nrow=length(uH),ncol=100)

    totalU=RBaM::prediction(X=list(htrep), # spaghetti of stage values
                            spagFiles='QRC_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.RC, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE) # propagate structural uncertainty ?
  }

  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.RC,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             pred=totalU, # list of predictions
             remnant = remnant_prior,
             doCalib=FALSE,
             doPred=TRUE)

  # Total uncertainty propagation
  env_QRC_TotalU=utils::read.table(file.path(temp.folder.RC,'QRC_TotalU.env'),header=TRUE)

  # Discharge simulation
  qsim <- resid.segm$Y1_sim

  # Residual calculation
  residuals <- resid.segm$Y1_res

  # Residual standard deviation
  residual_sd <- env_QRC_TotalU$Stdev

  # residual data frame
  ResultsResiduals=data.frame(time=data$time,
                              H=data$H,
                              uH=uH,
                              Q_obs=data$Q,
                              Q_sim=qsim,
                              Q_res=residuals,
                              uQ_obs=uQ,
                              uQ_sim=residual_sd
  )


  for(i in 1:ncontrols){
    local.parameters.temp=data.frame(summary.MCMC.MAP[,which(a.object[[i]]$name==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(paste0('b',i)==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(c.object[[i]]$name==colnames(summary.MCMC.MAP))],
                                     summary.MCMC.MAP[,which(k.object[[i]]$name==colnames(summary.MCMC.MAP))])


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

#' Recession model with two exponential and asymptotic
#'
#' Recession modelling following a the exponential function specified as M3 according to (Darienzo, 2022):
#' \deqn{h(t) = \alpha_1(k) \cdot \exp(-\lambda_1 \cdot t) + \alpha_2(k) \cdot \exp(-\lambda_2 \cdot t) + \beta(k)}
#' This model includes three recession-specific parameters: \eqn{\alpha_1}, \eqn{\alpha_2}, and \eqn{\beta}
#' and two stable parameters: \eqn{\lambda_1} and \eqn{\lambda_2}.
#'
#' Default values for this recession are shown in ‘Details’
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param hrec real vector, stage value of the recessions
#' @param uHrec real vector, uncertainty of stage value of the recessions
#' @param indx  integer, factor used to gather the data of a same recession
#' @param alpha1.object object, prior knowledge on the parameter representing the initial stage of the first exponential function
#' @param alpha2.object object, prior knowledge on the parameter representing the initial stage of the second exponential function
#' @param beta.object object, prior knowledge on the parameter representing the asymptotic stage
#' @param lambda1.object object, prior knowledge on the  parameter describing the fast runoff
#' @param lambda2.object object, prior knowledge on the parameter describing the slow emptying of the aquifer
#' @param Dataset.object object, dataset given by user
#' @param nCyclesrec  integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @details
#' By default, prior values are the same for all recession for the parameter alpha1, alpha2 and beta
#'
#' Default values are:
#' \describe{
#'   \item{Starting in 200, we assume that \eqn{\alpha_1 \sim U(0, 1000).}}
#'   \item{Starting in 50, we assume that \eqn{\alpha_2 \sim U(0, 500).}}
#'   \item{Starting in \eqn{\exp(-\log(0.5) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_1 \sim \text{LogNormal}(\mu = -\log(0.5) + \log(\log(2)), \sigma = 1).}}
#'   \item{Starting in \eqn{\exp(-\log(50) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_2 \sim \text{LogNormal}(\mu = -\log(80) + \log(\log(2)), \sigma = 0.5).}}
#'   \item{Starting in 1000, we assume that \eqn{\beta \sim U(-10000, 10000).}}
#' }
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : list, results after fitting the recession model to each recession observed
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage observed
#'        \item H_rec_sim: real value, stage simulated
#'        \item H_res: real value, residual between stage observed and simulated
#'        \item uH_rec: real value, uncertainty in stage observed (as a standard deviation)
#'        \item uH_sim: real value, uncertainty in stage simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the recession model
#'   \itemize{
#'        \item alpha1 : real value, parameter representing the initial stage
#'        \item lambda1 : real value, parameter describing the fast runoff
#'        \item alpha2 : real value, parameter representing the initial stage
#'        \item lambda2 : real value, parameter describing the slow emptying of the aquifer
#'        \item beta : real value, parameter representing the asymptotic stage
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRecession_M3 <- function(time_rec,hrec,uHrec,indx,
                            alpha1.object=NULL,
                            alpha2.object=NULL,
                            beta.object=NULL,
                            lambda1.object=NULL,
                            lambda2.object=NULL,
                            Dataset.object=NULL,
                            nCyclesrec=100,
                            burnrec=0.5,
                            nSlimrec=max(nCyclesrec/10,1),
                            temp.folder.Recession=file.path(tempdir(),'BaM','Recession')){

  # Create data frame to used in all calculation
  data=data.frame(time_rec=time_rec,hrec=hrec,uHrec=uHrec,indx=indx)
  # Read the number of recessions
  Ncurves=max(indx)

  if(any(!is.null(alpha1.object)|
         !is.null(alpha2.object)|
         !is.null(beta.object)|
         !is.null(lambda1.object)|
         !is.null(lambda2.object))&
     is.null(Dataset.object))stop('Prior knowledge about recession models are specified, dataset must be given too using the "dataset" function of RBaM')

  # Check and assign the dataset
  if(!is.null(Dataset.object)){
    if(class(Dataset.object)!='dataset')stop('Dataset.object must be created using the "dataset" function of RBaM before it can be used')
    # Use the calibration dataset given in the input
    colnames(Dataset.object$data) = c('time_rec','hrec','uHrec','indx')
    D=Dataset.object
  }else{
    # Define the calibration dataset by specifying
    D=RBaM::dataset(X=data['time_rec'],
                    Y=data['hrec'],
                    Yu=data['uHrec'],
                    VAR.indx=data['indx'],
                    data.dir=temp.folder.Recession)
  }

  # Check prior recessions
  # Prior knowledge based on user manual of BayDERS (Darienzo, 2022) (Fixed)
  if(!is.null(alpha1.object)){
    if(class(alpha1.object)!='parameter_VAR')stop('alpha1.object must be created using the "parameter_VAR" function of RBaM before it can be used')
    if(alpha1.object$name!='alpha1')stop("name of alpha1.object must be 'alpha1'")
    alpha1=alpha1.object
  }else{
    alpha1=RBaM::parameter_VAR(name='alpha1',
                               index='indx',
                               d=D,
                               init=rep(200,Ncurves), # first guesses
                               prior.dist=rep('Uniform',Ncurves), # prior distributions
                               prior.par=rep(list(c(0,1000)),Ncurves)) # prior parameters

  }

  if(!is.null(alpha2.object)){
    if(class(alpha2.object)!='parameter_VAR')stop('alpha2.object must be created using the "parameter_VAR" function of RBaM before it can be used')
    if(alpha2.object$name!='alpha2')stop("name of alpha2.object must be 'alpha2'")
    alpha2=alpha2.object
  }else{
    alpha2=RBaM::parameter_VAR(name='alpha2',
                               index='indx',
                               d=D,
                               init=rep(50,Ncurves), # first guesses
                               prior.dist=rep('Uniform',Ncurves), # prior distributions
                               prior.par=rep(list(c(0,500)),Ncurves)) # prior parameters
  }

  if(!is.null(beta.object)){
    if(class(beta.object)!='parameter_VAR')stop('beta.object must be created using the "parameter_VAR" function of RBaM before it can be used')
    if(beta.object$name!='beta')stop("name of beta.object must be 'beta'")
    beta=beta.object
  }else{
    beta=RBaM::parameter_VAR(name='beta',
                             index='indx',
                             d=D,
                             init=rep(1000,Ncurves), # first guesses
                             prior.dist=rep('Uniform',Ncurves), # prior distributions
                             prior.par=rep(list(c(-10000 , 10000)),Ncurves)) # prior parameters
  }

  if(!is.null(lambda1.object)){
    if(class(lambda1.object)!='parameter')stop('lambda1.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda1.object$name!='lambda1')stop("name of lambda1.object must be 'lambda1'")
    lambda1=lambda1.object
  }else{
    lambda1=RBaM::parameter(name='lambda1',
                            init=exp(-log(0.5)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c(-log((0.5))+log(log(2)),1))
  }

  if(!is.null(lambda2.object)){
    if(class(lambda2.object)!='parameter')stop('lambda2.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda2.object$name!='lambda2')stop("name of lambda2.object must be 'lambda2'")
    lambda2=lambda2.object
  }else{
    lambda2=RBaM::parameter(name='lambda2',
                            init=exp(-log(50)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c( -log(80)+log(log(2)),0.5))
  }

  priors=list(alpha1,lambda1,alpha2,lambda2,beta)
  # Stitch it all together into a model object
  M=RBaM::model(ID='Recession_h',
                nX=1,nY=1, # number of input/output variables
                par=priors) # list of model parameters

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
                                                                           prior.par = c(0,1000)))))  # Run BaM executable
  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.Recession,
             mcmc=mcmc_temp,
             cook = cook_temp,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             remnant = remnant_prior)

  # Save data object and model object
  save(D,file = file.path(temp.folder.Recession,'DataObject.RData'))
  save(M,file = file.path(temp.folder.Recession,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

  # First prediction : observed data :
  # Define a 'prediction' object for total predictive uncertainty only for observed stages
  MCMC    <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Cooking.txt"),header=TRUE)
  residus   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  indx_MAP=which(c('MaxPost')==rownames(summary.MCMC))
  indx_st_dev=which(c('St.Dev.')==rownames(summary.MCMC))

  # Way to handle prediction for VAR parameters :
  ResultsResiduals = c()

  for( i in 1:Ncurves){
    alpha1_nonVAR=RBaM::parameter(name=alpha1$name,
                                  init=alpha1$init[i],
                                  prior.dist=alpha1$prior[[i]]$dist,
                                  prior.par=alpha1$prior[[i]]$par)

    alpha2_nonVAR=RBaM::parameter(name=alpha2$name,
                                  init=alpha2$init[i],
                                  prior.dist=alpha2$prior[[i]]$dist,
                                  prior.par=alpha2$prior[[i]]$par)

    beta_nonVAR=RBaM::parameter(name=beta$name,
                                init=beta$init[i],
                                prior.dist=beta$prior[[i]]$dist,
                                prior.par=beta$prior[[i]]$par)

    M_nonVAR=RBaM::model(ID='Recession_h',
                  nX=1,
                  nY=1,
                  par=list(alpha1_nonVAR,
                           lambda1,
                           alpha2_nonVAR,
                           lambda2,
                           beta_nonVAR))
    # Columns names to use during prediction
    columns_indx=c(paste0('alpha1_',i),
                   'lambda1',
                   paste0('alpha2_',i),
                   'lambda2',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    data.rec.spec=data[which(data$indx==i),]

    # Define a 'prediction' object for total predictive uncertainty
    totalU= RBaM::prediction(X=data.rec.spec['hrec'], # stage values
                             spagFiles='hrec_TotalU.spag', # file where predictions are saved
                             data.dir=temp.folder.Recession, # a copy of data files will be saved here
                             doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                             doStructural=TRUE, # propagate structural uncertainty ?
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM::BaM(mod=M_nonVAR,
               data=D,
               workspace = temp.folder.Recession,
               dir.exe = file.path(find.package("RBaM"), "bin"),
               pred=totalU,
               remnant = remnant_prior,
               doCalib=FALSE,
               doPred=TRUE)

    # Total uncertainty propagation
    env_hrec_TotalU=utils::read.table(file.path(temp.folder.Recession,'hrec_TotalU.env'),header=TRUE)

    # Discharge simulation
    hrecsim <-  residus['Y1_sim'][which(data$indx==i),]

    # Residual calculation
    residuals <- residus['Y1_res'][which(data$indx==i),]

    # Residual standard deviation
    residual_sd <- env_hrec_TotalU$Stdev

    # residual data frame
    ResultsResiduals[[i]]=data.frame(time=data.rec.spec['time_rec'],
                                     H=data.rec.spec['hrec'],
                                     H_rec_sim=hrecsim,
                                     H_res=residuals,
                                     uH_rec=data.rec.spec['uHrec'],
                                     uH_sim=residual_sd)

    # Create a table to store the parameters for all recessions
    local.parameters.temp=data.frame(alpha1=summary.MCMC[indx_MAP,
                                                         which(paste0('alpha1_',i)==colnames(summary.MCMC))],
                                     u_alpha1=summary.MCMC[indx_st_dev,
                                                           which(paste0('alpha1_',i)==colnames(summary.MCMC))],
                                     alpha2=summary.MCMC[indx_MAP,
                                                         which(paste0('alpha2_',i)==colnames(summary.MCMC))],
                                     u_alpha2=summary.MCMC[indx_st_dev,
                                                           which(paste0('alpha2_',i)==colnames(summary.MCMC))],
                                     beta=summary.MCMC[indx_MAP,
                                                       which(paste0('beta_',i)==colnames(summary.MCMC))],
                                     u_beta=summary.MCMC[indx_st_dev,
                                                         which(paste0('beta_',i)==colnames(summary.MCMC))])

    if(i==1){
      local.parameters <- local.parameters.temp
    }else{
      local.parameters <- rbind(local.parameters,local.parameters.temp)
    }
  }

  parameters= data.frame(
    local.parameters['alpha1'],
    local.parameters['u_alpha1'],
    lambda1=summary.MCMC[indx_MAP,
                        which('lambda1'==colnames(summary.MCMC))],
    u_lambda1=summary.MCMC[indx_st_dev,
                          which('lambda1'==colnames(summary.MCMC))],
    local.parameters['alpha2'],
    local.parameters['u_alpha2'],
    lambda2=summary.MCMC[indx_MAP,
                        which('lambda2'==colnames(summary.MCMC))],
    u_lambda2=summary.MCMC[indx_st_dev,
                          which('lambda2'==colnames(summary.MCMC))],
    local.parameters['beta'],
    local.parameters['u_beta'],
    gamma1=summary.MCMC[indx_MAP,]$Y1_gamma1,
    u_gamma1=summary.MCMC[indx_st_dev,]$Y1_gamma1,
    gamma2=summary.MCMC[indx_MAP,]$Y1_gamma2,
    u_gamma2=summary.MCMC[indx_st_dev,]$Y1_gamma2)

  # ggplot(ResultsResiduals[[39]],aes(x=time_rec))+
  #   geom_point(aes(y=hrec),col='black')+
  #   geom_point(aes(y=H_rec_sim),col='red')

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.Recession,
                       dir.destination=file.path(temp.folder.Recession, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.Recession ,
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

#' Exponential recession model with five parameters
#'
#' Equation of recession model below:
#' \deqn{h(t) =  \alpha_1 \cdot \exp(-\lambda_1 \cdot t -\lambda_2 \cdot t^2 - \lambda_3 \cdot t^(0.5)) + \beta(k)}
#' This model includes one recession-specific parameters: \eqn{\beta(k)}
#' and three stable parameters: \eqn{\alpha_1}, \eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\lambda_3}.
#'
#' Default values for this recession are shown in ‘Details’
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param hrec real vector, stage value of the recessions
#' @param uHrec real vector, uncertainty of stage value of the recessions
#' @param indx  integer, factor used to gather the data of a same recession
#' @param beta.object object, prior knowledge on the parameter representing the asymptotic stage
#' @param alpha1.object object, prior knowledge on the parameter representing the initial stage of the first exponential function
#' @param lambda1.object object, prior knowledge on the first parameter of the exponential function
#' @param lambda2.object object, prior knowledge on the second parameter of the exponential function
#' @param lambda3.object object, prior knowledge on the third parameter of the exponential function
#' @param Dataset.object object, dataset given by user
#' @param nCyclesrec  integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @details
#' By default, prior values are the same for all recession for the parameter alpha1, alpha2 and beta
#'
#' Default values are:
#' \describe{
#'   \item{Starting in 200, we assume that \eqn{\alpha_1 \sim U(0, 1000).}}
#'   \item{Starting in \eqn{\exp(-\log(0.5) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_1 \sim \text{LogNormal}(\mu = -\log(0.5) + \log(\log(2)), \sigma = 1).}}
#'   \item{Starting in \eqn{\exp(-\log(50) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_2 \sim \text{LogNormal}(\mu = -\log(80) + \log(\log(2)), \sigma = 0.5).}}
#'   \item{Starting in \eqn{\exp(-\log(50) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_3 \sim \text{LogNormal}(\mu = -\log(80) + \log(\log(2)), \sigma = 0.5).}}
#'   \item{Starting in 1000, we assume that \eqn{\beta \sim U(-10000, 10000).}}
#' }
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : list, results after fitting the recession model to each recession observed
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage observed
#'        \item H_rec_sim: real value, stage simulated
#'        \item H_res: real value, residual between stage observed and simulated
#'        \item uH_rec: real value, uncertainty in stage observed (as a standard deviation)
#'        \item uH_sim: real value, uncertainty in stage simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the recession model
#'   \itemize{
#'        \item alpha1 : real value, parameter representing the initial stage
#'        \item beta : real value, parameter representing the asymptotic stage
#'        \item lambda1 : real value, parameter describing the first parameter of the exponential function
#'        \item lambda2 : real value, parameter describing the second parameter of the exponential function
#'        \item lambda3 : real value, parameter describing the third parameter of the exponential function
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRecession_BR1 <- function(time_rec,hrec,uHrec,indx,
                             beta.object=NULL,
                             alpha1.object=NULL,
                             lambda1.object=NULL,
                             lambda2.object=NULL,
                             lambda3.object=NULL,
                             Dataset.object=NULL,
                             nCyclesrec=100,
                             burnrec=0.5,
                             nSlimrec=max(nCyclesrec/10,1),
                             temp.folder.Recession=file.path(tempdir(),'BaM','Recession')){

  # Create data frame to used in all calculation
  data=data.frame(time_rec=time_rec,hrec=hrec,uHrec=uHrec,indx=indx)
  # Read the number of recessions
  Ncurves=max(indx)

  if(any(!is.null(beta.object)|
         !is.null(alpha1.object)|
         !is.null(lambda1.object)|
         !is.null(lambda2.object)|
         !is.null(lambda3.object))&
     is.null(Dataset.object))stop('Prior knowledge about recession models are specified, dataset must be given too using the "dataset" function of RBaM')

  # Check and assign the dataset
  if(!is.null(Dataset.object)){
    if(class(Dataset.object)!='dataset')stop('Dataset.object must be created using the "dataset" function of RBaM before it can be used')
    # Use the calibration dataset given in the input
    colnames(Dataset.object$data) = c('time_rec','hrec','uHrec','indx')
    D=Dataset.object
  }else{
    # Define the calibration dataset by specifying
    D=RBaM::dataset(X=data['time_rec'],
                    Y=data['hrec'],
                    Yu=data['uHrec'],
                    VAR.indx=data['indx'],
                    data.dir=temp.folder.Recession)
  }

  # Check prior recessions
  # Prior knowledge Fixed
  if(!is.null(alpha1.object)){
    if(class(alpha1.object)!='parameter')stop('alpha1.object must be created using the "parameter" function of RBaM before it can be used')
    if(alpha1.object$name!='alpha1')stop("name of alpha1.object must be 'alpha1'")
    alpha1=alpha1.object
  }else{
    alpha1=RBaM::parameter(name='alpha1',
                           init=200,
                           prior.dist='Uniform',
                           prior.par=c(0,1000))

  }

 if(!is.null(beta.object)){
    if(class(beta.object)!='parameter_VAR')stop('beta.object must be created using the "parameter_VAR" function of RBaM before it can be used')
    if(beta.object$name!='beta')stop("name of beta.object must be 'beta'")
    beta=beta.object
  }else{
    beta=RBaM::parameter_VAR(name='beta',
                             index='indx',
                             d=D,
                             init=rep(1000,Ncurves), # first guesses
                             prior.dist=rep('Uniform',Ncurves), # prior distributions
                             prior.par=rep(list(c(-10000 , 10000)),Ncurves)) # prior parameters
  }

  if(!is.null(lambda1.object)){
    if(class(lambda1.object)!='parameter')stop('lambda1.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda1.object$name!='lambda1')stop("name of lambda1.object must be 'lambda1'")
    lambda1=lambda1.object
  }else{
    lambda1=RBaM::parameter(name='lambda1',
                            init=exp(-log(0.5)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c(-log((0.5))+log(log(2)),1))
  }

  if(!is.null(lambda2.object)){
    if(class(lambda2.object)!='parameter')stop('lambda2.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda2.object$name!='lambda2')stop("name of lambda2.object must be 'lambda2'")
    lambda2=lambda2.object
  }else{
    lambda2=RBaM::parameter(name='lambda2',
                            init=exp(-log(50)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c( -log(80)+log(log(2)),0.5))
  }

  if(!is.null(lambda3.object)){
    if(class(lambda3.object)!='parameter')stop('lambda3.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda3.object$name!='lambda3')stop("name of lambda3.object must be 'lambda3'")
    lambda3=lambda3.object
  }else{
    lambda3=RBaM::parameter(name='lambda3',
                            init=exp(-log(50)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c( -log(80)+log(log(2)),0.5))
  }

  priors=list(alpha1,lambda1,lambda2,lambda3,beta)
  # Stitch it all together into a model object
  M=RBaM::model(ID='Recession_h',
                nX=1,nY=1, # number of input/output variables
                par=priors) # list of model parameters

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
                                                                           prior.par = c(0,1000)))))  # Run BaM executable
  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.Recession,
             mcmc=mcmc_temp,
             cook = cook_temp,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             remnant = remnant_prior)

  # Save data object and model object
  save(D,file = file.path(temp.folder.Recession,'DataObject.RData'))
  save(M,file = file.path(temp.folder.Recession,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

  # First prediction : observed data :
  # Define a 'prediction' object for total predictive uncertainty only for observed stages
  MCMC    <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Cooking.txt"),header=TRUE)
  residus   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  indx_MAP=which(c('MaxPost')==rownames(summary.MCMC))
  indx_st_dev=which(c('St.Dev.')==rownames(summary.MCMC))

  # Way to handle prediction for VAR parameters :
  ResultsResiduals = c()

  for( i in 1:Ncurves){
    beta_nonVAR=RBaM::parameter(name=beta$name,
                                init=beta$init[i],
                                prior.dist=beta$prior[[i]]$dist,
                                prior.par=beta$prior[[i]]$par)

    M_nonVAR=RBaM::model(ID='Recession_h',
                         nX=1,
                         nY=1,
                         par=list(alpha1,
                                  lambda1,
                                  lambda2,
                                  lambda3,
                                  beta_nonVAR))
    # Columns names to use during prediction
    columns_indx=c('alpha1',
                   'lambda1',
                   'lambda2',
                   'lambda3',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    data.rec.spec=data[which(data$indx==i),]

    # Define a 'prediction' object for total predictive uncertainty
    totalU= RBaM::prediction(X=data.rec.spec['hrec'], # stage values
                             spagFiles='hrec_TotalU.spag', # file where predictions are saved
                             data.dir=temp.folder.Recession, # a copy of data files will be saved here
                             doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                             doStructural=TRUE, # propagate structural uncertainty ?
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM::BaM(mod=M_nonVAR,
              data=D,
              workspace = temp.folder.Recession,
              dir.exe = file.path(find.package("RBaM"), "bin"),
              pred=totalU,
              remnant = remnant_prior,
              doCalib=FALSE,
              doPred=TRUE)

    # Total uncertainty propagation
    env_hrec_TotalU=utils::read.table(file.path(temp.folder.Recession,'hrec_TotalU.env'),header=TRUE)

    # Discharge simulation
    hrecsim <-  residus['Y1_sim'][which(data$indx==i),]

    # Residual calculation
    residuals <- residus['Y1_res'][which(data$indx==i),]

    # Residual standard deviation
    residual_sd <- env_hrec_TotalU$Stdev

    # residual data frame
    ResultsResiduals[[i]]=data.frame(time=data.rec.spec['time_rec'],
                                     H=data.rec.spec['hrec'],
                                     H_rec_sim=hrecsim,
                                     H_res=residuals,
                                     uH_rec=data.rec.spec['uHrec'],
                                     uH_sim=residual_sd)

    # Create a table to store the parameters for all recessions
    local.parameters.temp=data.frame(beta=summary.MCMC[indx_MAP,
                                                       which(paste0('beta_',i)==colnames(summary.MCMC))],
                                     u_beta=summary.MCMC[indx_st_dev,
                                                         which(paste0('beta_',i)==colnames(summary.MCMC))])

    if(i==1){
      local.parameters <- local.parameters.temp
    }else{
      local.parameters <- rbind(local.parameters,local.parameters.temp)
    }
  }

  parameters= data.frame(
    alpha1=summary.MCMC[indx_MAP,
                        which('alpha1'==colnames(summary.MCMC))],
    u_alpha1=summary.MCMC[indx_st_dev,
                          which('alpha1'==colnames(summary.MCMC))],
    lambda1=summary.MCMC[indx_MAP,
                        which('lambda1'==colnames(summary.MCMC))],
    u_lambda1=summary.MCMC[indx_st_dev,
                          which('lambda1'==colnames(summary.MCMC))],
    lambda2=summary.MCMC[indx_MAP,
                         which('lambda2'==colnames(summary.MCMC))],
    u_lambda2=summary.MCMC[indx_st_dev,
                           which('lambda2'==colnames(summary.MCMC))],
    lambda3=summary.MCMC[indx_MAP,
                         which('lambda3'==colnames(summary.MCMC))],
    u_lambda3=summary.MCMC[indx_st_dev,
                           which('lambda3'==colnames(summary.MCMC))],
    local.parameters['beta'],
    local.parameters['u_beta'],
    gamma1=summary.MCMC[indx_MAP,]$Y1_gamma1,
    u_gamma1=summary.MCMC[indx_st_dev,]$Y1_gamma1,
    gamma2=summary.MCMC[indx_MAP,]$Y1_gamma2,
    u_gamma2=summary.MCMC[indx_st_dev,]$Y1_gamma2)

  # ggplot(ResultsResiduals[[39]],aes(x=time_rec))+
  #   geom_point(aes(y=hrec),col='black')+
  #   geom_point(aes(y=H_rec_sim),col='red')

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.Recession,
                       dir.destination=file.path(temp.folder.Recession, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.Recession ,
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

#' Exponential recession model with five parameters
#'
#' Equation of recession model below:
#' \deqn{h(t) =  \alpha_1 \cdot \exp(-\lambda_1 \cdot t -\lambda_2 \cdot t^2 - \lambda_3 \cdot t^(0.5)) + \beta(k)}
#' This model includes one recession-specific parameters: \eqn{\beta(k)}
#' and three stable parameters: \eqn{\alpha_1}, \eqn{\lambda_1}, \eqn{\lambda_2}, \eqn{\lambda_3}.
#'
#' Default values for this recession are shown in ‘Details’
#'
#' @param time_rec real vector, recession duration relative to the first data detected during the recession
#' @param hrec real vector, stage value of the recessions
#' @param uHrec real vector, uncertainty of stage value of the recessions
#' @param indx  integer, factor used to gather the data of a same recession
#' @param beta.object object, prior knowledge on the parameter representing the asymptotic stage
#' @param alpha1.object object, prior knowledge on the parameter representing the initial stage of the first exponential function
#' @param lambda1.object object, prior knowledge on the first parameter of the exponential function
#' @param lambda2.object object, prior knowledge on the second parameter of the exponential function
#' @param lambda3.object object, prior knowledge on the third parameter of the exponential function
#' @param Dataset.object object, dataset given by user
#' @param nCyclesrec  integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles
#' @param burnrec real between 0 (included) and 1 (excluded), MCMC burning factor
#' @param nSlimrec integer, MCMC slim step
#' @param temp.folder.Recession directory, temporary directory to write computations
#'
#' @details
#' By default, prior values are the same for all recession for the parameter alpha1, alpha2 and beta
#'
#' Default values are:
#' \describe{
#'   \item{Starting in 200, we assume that \eqn{\alpha_1 \sim U(0, 1000).}}
#'   \item{Starting in \eqn{\exp(-\log(0.5) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_1 \sim \text{LogNormal}(\mu = -\log(0.5) + \log(\log(2)), \sigma = 1).}}
#'   \item{Starting in \eqn{\exp(-\log(50) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_2 \sim \text{LogNormal}(\mu = -\log(80) + \log(\log(2)), \sigma = 0.5).}}
#'   \item{Starting in \eqn{\exp(-\log(50) + \log(\log(2)))}, we assume that
#'     \eqn{\lambda_3 \sim \text{LogNormal}(\mu = -\log(80) + \log(\log(2)), \sigma = 0.5).}}
#'   \item{Starting in 1000, we assume that \eqn{\beta \sim U(-10000, 10000).}}
#' }
#'
#' @return List with the following components :
#' \enumerate{
#'   \item ResultsResiduals : list, results after fitting the recession model to each recession observed
#'   \itemize{
#'        \item time: real value, time
#'        \item H: real value, stage observed
#'        \item H_rec_sim: real value, stage simulated
#'        \item H_res: real value, residual between stage observed and simulated
#'        \item uH_rec: real value, uncertainty in stage observed (as a standard deviation)
#'        \item uH_sim: real value, uncertainty in stage simulated (as a standard deviation)
#'        }
#'   \item parameters : data frame, parameters of the recession model
#'   \itemize{
#'        \item alpha1 : real value, parameter representing the initial stage
#'        \item beta : real value, parameter representing the asymptotic stage
#'        \item lambda1 : real value, parameter describing the first parameter of the exponential function
#'        \item lambda2 : real value, parameter describing the second parameter of the exponential function
#'        \item lambda3 : real value, parameter describing the third parameter of the exponential function
#'        \item gamma1 ; real value, parameter describing the standard deviation of structural errors (relative to discharge)
#'        \item gamma2 : real value, parameter describing the standard deviation of structural errors (constant error)
#'        }
#' }
#' @export
fitRecession_BR2 <- function(time_rec,hrec,uHrec,indx,
                             beta.object=NULL,
                             alpha1.object=NULL,
                             lambda1.object=NULL,
                             c1.object=NULL,
                             Dataset.object=NULL,
                             nCyclesrec=100,
                             burnrec=0.5,
                             nSlimrec=max(nCyclesrec/10,1),
                             temp.folder.Recession=file.path(tempdir(),'BaM','Recession')){

  # Create data frame to used in all calculation
  data=data.frame(time_rec=time_rec,hrec=hrec,uHrec=uHrec,indx=indx)
  # Read the number of recessions
  Ncurves=max(indx)

  if(any(!is.null(beta.object)|
         !is.null(alpha1.object)|
         !is.null(lambda1.object)|
         !is.null(c1.object))&
     is.null(Dataset.object))stop('Prior knowledge about recession models are specified, dataset must be given too using the "dataset" function of RBaM')

  # Check and assign the dataset
  if(!is.null(Dataset.object)){
    if(class(Dataset.object)!='dataset')stop('Dataset.object must be created using the "dataset" function of RBaM before it can be used')
    # Use the calibration dataset given in the input
    colnames(Dataset.object$data) = c('time_rec','hrec','uHrec','indx')
    D=Dataset.object
  }else{
    # Define the calibration dataset by specifying
    D=RBaM::dataset(X=data['time_rec'],
                    Y=data['hrec'],
                    Yu=data['uHrec'],
                    VAR.indx=data['indx'],
                    data.dir=temp.folder.Recession)
  }

  # Check prior recessions
  # Prior knowledge Fixed
  if(!is.null(alpha1.object)){
    if(class(alpha1.object)!='parameter')stop('alpha1.object must be created using the "parameter" function of RBaM before it can be used')
    if(alpha1.object$name!='alpha1')stop("name of alpha1.object must be 'alpha1'")
    alpha1=alpha1.object
  }else{
    alpha1=RBaM::parameter(name='alpha1',
                           init=200,
                           prior.dist='Uniform',
                           prior.par=c(0,1000))

  }

  if(!is.null(beta.object)){
    if(class(beta.object)!='parameter_VAR')stop('beta.object must be created using the "parameter_VAR" function of RBaM before it can be used')
    if(beta.object$name!='beta')stop("name of beta.object must be 'beta'")
    beta=beta.object
  }else{
    beta=RBaM::parameter_VAR(name='beta',
                             index='indx',
                             d=D,
                             init=rep(1000,Ncurves), # first guesses
                             prior.dist=rep('Uniform',Ncurves), # prior distributions
                             prior.par=rep(list(c(-10000 , 10000)),Ncurves)) # prior parameters
  }

  if(!is.null(lambda1.object)){
    if(class(lambda1.object)!='parameter')stop('lambda1.object must be created using the "parameter" function of RBaM before it can be used')
    if(lambda1.object$name!='lambda1')stop("name of lambda1.object must be 'lambda1'")
    lambda1=lambda1.object
  }else{
    lambda1=RBaM::parameter(name='lambda1',
                            init=exp(-log(0.5)+log(log(2))),
                            prior.dist='LogNormal',
                            prior.par=c(-log((0.5))+log(log(2)),1))
  }

  if(!is.null(c1.object)){
    if(class(c1.object)!='parameter')stop('c1.object must be created using the "parameter" function of RBaM before it can be used')
    if(c1.object$name!='c1')stop("name of c1.object must be 'c1'")
    c1=c1.object
  }else{
    c1=RBaM::parameter(name='c1',
                       init=0.5,
                       prior.dist='Uniform',
                       prior.par=c(0,1))
  }

  priors=list(alpha1,lambda1,c1,beta)
  # Stitch it all together into a model object
  M=RBaM::model(ID='Recession_h',
                nX=1,nY=1, # number of input/output variables
                par=priors) # list of model parameters

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
                                                                           prior.par = c(0,1000)))))  # Run BaM executable
  RBaM:: BaM(mod=M,
             data=D,
             workspace = temp.folder.Recession,
             mcmc=mcmc_temp,
             cook = cook_temp,
             dir.exe = file.path(find.package("RBaM"), "bin"),
             remnant = remnant_prior)

  # Save data object and model object
  save(D,file = file.path(temp.folder.Recession,'DataObject.RData'))
  save(M,file = file.path(temp.folder.Recession,'ModelObject.RData'))

  # PREDICTIONS : two steps
  # First prediction : estimate total uncertainty of simulation (u_sim = u_total) at observed stages to returned as u_sim for segmentation
  # Second prediction : estimate total uncertainty of simulation discretized at Hgrid to plot (manage in PlotRCPrediction function)

  # First prediction : observed data :
  # Define a 'prediction' object for total predictive uncertainty only for observed stages
  MCMC    <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Cooking.txt"),header=TRUE)
  residus   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Residuals.txt"),header=TRUE)
  summary.MCMC   <- utils::read.table(file=file.path(temp.folder.Recession,"Results_Summary.txt"),header=TRUE)

  indx_MAP=which(c('MaxPost')==rownames(summary.MCMC))
  indx_st_dev=which(c('St.Dev.')==rownames(summary.MCMC))

  # Way to handle prediction for VAR parameters :
  ResultsResiduals = c()

  for( i in 1:Ncurves){
    beta_nonVAR=RBaM::parameter(name=beta$name,
                                init=beta$init[i],
                                prior.dist=beta$prior[[i]]$dist,
                                prior.par=beta$prior[[i]]$par)

    M_nonVAR=RBaM::model(ID='Recession_h',
                         nX=1,
                         nY=1,
                         par=list(alpha1,
                                  lambda1,
                                  c1,
                                  beta_nonVAR))
    # Columns names to use during prediction
    columns_indx=c('alpha1',
                   'lambda1',
                   'c1',
                   paste0('beta_',i),
                   'Y1_gamma1',
                   'Y1_gamma2',
                   'LogPost')

    parSamples=MCMC[columns_indx]

    data.rec.spec=data[which(data$indx==i),]

    # Define a 'prediction' object for total predictive uncertainty
    totalU= RBaM::prediction(X=data.rec.spec['hrec'], # stage values
                             spagFiles='hrec_TotalU.spag', # file where predictions are saved
                             data.dir=temp.folder.Recession, # a copy of data files will be saved here
                             doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                             doStructural=TRUE, # propagate structural uncertainty ?
                             parSamples=parSamples # pass the reduced MCMC data frame to use for this prediction
    )

    RBaM::BaM(mod=M_nonVAR,
              data=D,
              workspace = temp.folder.Recession,
              dir.exe = file.path(find.package("RBaM"), "bin"),
              pred=totalU,
              remnant = remnant_prior,
              doCalib=FALSE,
              doPred=TRUE)

    # Total uncertainty propagation
    env_hrec_TotalU=utils::read.table(file.path(temp.folder.Recession,'hrec_TotalU.env'),header=TRUE)

    # Discharge simulation
    hrecsim <-  residus['Y1_sim'][which(data$indx==i),]

    # Residual calculation
    residuals <- residus['Y1_res'][which(data$indx==i),]

    # Residual standard deviation
    residual_sd <- env_hrec_TotalU$Stdev

    # residual data frame
    ResultsResiduals[[i]]=data.frame(time=data.rec.spec['time_rec'],
                                     H=data.rec.spec['hrec'],
                                     H_rec_sim=hrecsim,
                                     H_res=residuals,
                                     uH_rec=data.rec.spec['uHrec'],
                                     uH_sim=residual_sd)

    # Create a table to store the parameters for all recessions
    local.parameters.temp=data.frame(beta=summary.MCMC[indx_MAP,
                                                       which(paste0('beta_',i)==colnames(summary.MCMC))],
                                     u_beta=summary.MCMC[indx_st_dev,
                                                         which(paste0('beta_',i)==colnames(summary.MCMC))])

    if(i==1){
      local.parameters <- local.parameters.temp
    }else{
      local.parameters <- rbind(local.parameters,local.parameters.temp)
    }
  }

  parameters= data.frame(
    alpha1=summary.MCMC[indx_MAP,
                        which('alpha1'==colnames(summary.MCMC))],
    u_alpha1=summary.MCMC[indx_st_dev,
                          which('alpha1'==colnames(summary.MCMC))],
    lambda1=summary.MCMC[indx_MAP,
                         which('lambda1'==colnames(summary.MCMC))],
    u_lambda1=summary.MCMC[indx_st_dev,
                           which('lambda1'==colnames(summary.MCMC))],
    c1=summary.MCMC[indx_MAP,
                         which('c1'==colnames(summary.MCMC))],
    u_c1=summary.MCMC[indx_st_dev,
                           which('c1'==colnames(summary.MCMC))],
    local.parameters['beta'],
    local.parameters['u_beta'],
    gamma1=summary.MCMC[indx_MAP,]$Y1_gamma1,
    u_gamma1=summary.MCMC[indx_st_dev,]$Y1_gamma1,
    gamma2=summary.MCMC[indx_MAP,]$Y1_gamma2,
    u_gamma2=summary.MCMC[indx_st_dev,]$Y1_gamma2)

  # ggplot(ResultsResiduals[[39]],aes(x=time_rec))+
  #   geom_point(aes(y=hrec),col='black')+
  #   geom_point(aes(y=H_rec_sim),col='red')

  # Save results from first prediction
  copy_files_to_folder(dir.source=temp.folder.Recession,
                       dir.destination=file.path(temp.folder.Recession, 'Residual'))

  # Second prediction : Hgrid :
  invisible(remove_files(dir.source = temp.folder.Recession ,
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


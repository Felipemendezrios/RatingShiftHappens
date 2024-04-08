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


#' #' Fit rating curve using BaRatin model with fixed parameters
#' #'
#' #' Rating curve estimated by using BaRatin approach with the parameter c fixed to 1.67 and only one hydraulic control type channel. The input data must be already ordered.
#' #'
#' #' @param time real vector, time
#' #' @param H real vector, stage
#' #' @param Q real vector, discharge
#' #' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' #' @param b.distr
#' #' @param a.distr
#' #' @param b.prior
#' #' @param st_b.prior
#' #' @param Bc.prior
#' #' @param KS.prior
#' #' @param S0.prior
#' #' @param st_Bc.prior
#' #' @param st_KS.prior
#' #' @param st_S0.prior
#' #' @param remnant.err.model
#' #' @param g1.prior
#' #' @param g2.prior
#' #' @param g1.distr.type
#' #' @param g2.distr.type
#' #' @param hmax_Xtra
#' #'
#' #' @return Not yet
# fitRC_BaRatin_simplifie=function(time,H,Q,uQ,
#                                  b.distr = 'Gaussian' ,a.distr = 'LogNormal',
#                                  b.prior, st_b.prior,
#                                  Bc.prior, KS.prior, S0.prior,
#                                  st_Bc.prior, st_KS.prior, st_S0.prior,
#                                  remnant.err.model = 'Linear',
#                                  g1.prior = c(0, 1000, 0.1) , g2.prior = c(0, 100, 0.1),
#                                  g1.distr.type='Uniform', g2.distr.type='Uniform',
#                                  hmax_Xtra=1.5*max(H)){
#   # Hardcoded variables
#   hyr_contr_matrix=matrix(c(1),nrow=1,ncol=1)
#   c = 1.67
#
#   if(is.null(check_square_matrix(hyr_contr_matrix)))stop('Matrix must be square. Try another hydraulic matrix control')
#   if(is.null(check_vector_lengths(hyr_contr_matrix,
#                                   b.prior, Bc.prior, KS.prior, S0.prior)))stop('Prior information of the parameters does not have the same length')
#   if(is.null(check_vector_lengths(hyr_contr_matrix,
#                                   st_b.prior,st_Bc.prior,st_KS.prior, st_S0.prior)))stop('Prior information of the parameters does not have the same length')
#
#   # How to handle uncertainty of the parameters? because depends on the distribution, number of parameters can vary
#   if(is.null(check_param_distribution(b.distr,st_b.prior)))stop('The number of the parameters does not match with the specified distribution for st_b.prior')
#   if(is.null(check_param_distribution(a.distr,st_Bc.prior)))stop('The number of the parameters does not match with the specified distribution for st_Bc.prior')
#   if(is.null(check_param_distribution(a.distr,st_KS.prior)))stop('The number of the parameters does not match with the specified distribution for st_KS.prior')
#   if(is.null(check_param_distribution(a.distr,st_S0.prior)))stop('The number of the parameters does not match with the specified distribution for st_S0.prior')
#   if(is.null(check_param_distribution(g1.distr.type,g1.prior)))stop('The number of the parameters does not match with the specified distribution for g1.prior')
#   if(is.null(check_param_distribution(g2.distr.type,g2.prior)))stop('The number of the parameters does not match with the specified distribution for g2.prior')
#
#   npar = 3*ncol(hyr_contr_matrix) # Three parameters by control
#   priors <- vector(mode = 'list',length = npar)
#
#   for(i in 1 :ncol(hyr_contr_matrix)){
#     b1=RBaM::parameter(name=paste0('b',i),
#                        init=b.prior[i],
#                        prior.dist = b.distr ,
#                        prior.par = st_b.prior)
#       parameter(name='k1',init=-0.5,prior.dist='Uniform',prior.par=c(-1.5,0))
#
#   }
#   # put all parameters in priors
#   data=data.frame(time=time,H=H,Q=Q)
#
#   ##ending
#   # residual data frame
#   ResultsResiduals=data.frame(time=data$time,
#                               H=data$H,
#                               Q_obs=data$Q,
#                               Q_sim=qsim,
#                               Q_res=residuals,
#                               uQ_obs=uQ,
#                               uQ_sim=residual_sd
#   )
#   parameters=data.frame(a=Q0,b=mu)
#
#   return(list(ResultsResiduals=ResultsResiduals,
#               parameters=parameters))
# }

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
  if(length(time)<=2){ # because second degree polynomial by default (loess function)
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



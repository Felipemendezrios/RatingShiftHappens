#' fit rating curve LOESS model
#'
#' LOESS model used to estimate a simple rating curve
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#'
#' @return data frame with the following components :
#' \enumerate{
#'   \item time: real value, time
#'   \item H: real value, stage
#'   \item Q_obs: real value, discharge observed
#'   \item Q_sim: real value, discharge simulated
#'   \item Q_res: real value, residual between discharge observed and simulated
#'   \item uQ_obs: real value, uncertainty in discharge observed (as a standard deviation)
#'   \item uQ_sim: real value, uncertainty in discharge simulated (as a standard deviation)
#' }
#' @export
#'
#' @examples
#' # Dataset
#' subset = RhoneRiver[1:20,]
#'
#' # LOESS regression to estimate a simple rating curve
#' fit=fitRC_loess(time=subset$Time,H=subset$H,Q=subset$Q,uQ=subset$uQ)
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
    warning('NULL was returned because it not possible to perform LOESS regression with only two points.
    A second degree polynomial requires at least three points for prediction')
    return(NA)
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
  return(ResultsResiduals)
}


#' BaRatin model
#'
#' Rating curve estimated by using BaRatin approach
#'
#' @param time real vector, time
#' @param H real vector, stage
#' @param Q real vector, discharge
#' @param uQ real vector, uncertainty in discharge (as a standard deviation)
#' @param matrix real data frame, hydraulic controls specified as bonnifait matrix
#' @param priors real list, prior information of each parameter
#' @param remnant real vector, remnant model information
#'
#' @return real data frame
fitRC_BaRatin=function(time,H,Q,uQ,matrix,priors,remnant){
  stop("fitRC_BaRatin has not yet been implemented")
}




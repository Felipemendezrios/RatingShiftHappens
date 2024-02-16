fitRC_loess<-function(time,H,Q,uQ){
  if(length(time)<=2){ # because second degree polynomial by default (loess function)
    warning('NULL was returned because it not possible to perform LOESS regression with only two points.
    A second degree polynomial requires at least three points for prediction')
    return(NA)
  }
  data=data.frame(time=time,H=H,Q=Q)

  # LOESS model for estimating rating curve
  mod <- loess(Q ~ H, data = data)

  # Prediction of discharge Q_rc from the water level (H)
  qsim <- predict(mod, newdata = data$H)

  # Residual calculation
  residuals <- data$Q - qsim

  # Residual standard deviation
  residual_sd <- sd(residuals)  # incorrect interpreted as predictive uncertainty

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


fitRC_BaRatin=function(time,H,Q,uQ,matrix,priors,remnant){
  stop("fitRC_BaRatin has not yet been implemented")
}



fitRC <- function(time,H,Q,uQ=0,funk=fitRC_loess,...){

  if(any(is.na(time)) | any(is.na(H)) | any(is.na(Q)) | any(is.na(uQ))){
    stop('Missing values not allowed in time, stage, discharge or uncertainty')
  }

  check <- check_vector_lengths(time,H,Q,uQ)
  if(is.null(check)){
    stop('time, hauteur, discharge or uncertainty do not have the same length')
  }
  return(results = funk(time,H,Q,uQ,...))
}




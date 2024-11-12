#' RatingShiftHappens catalog
#'
#' Rating curve and models with their equation available
#'
#' @param printOnly Logical, should the catalog be returned or only printed?
#' @return If \code{printOnly==FALSE}, a list with the following fields:
#' \enumerate{
#'   \item models: available models
#'   \item RC_Equations: available RC equations
#' }
#' @examples
#' catalog <- GetCatalog_RC()
#' GetCatalog_RC(printOnly=TRUE)
#' @export
GetCatalog_RC<-function(printOnly=FALSE){
  # available models
  models=c('fitRC_loess',
           'fitRC_BaRatinBAC',
           'fitRC_BaRatinKAC',
           'fitRC_exponential',
           'fitRC_LinearRegression')

  RC_Equations=c('Loess_Equation',
                 'BaRatinBAC_Equation',
                 'BaRatinKAC_Equation',
                 'Exponential_Equation',
                 'LinearRegression_Equation')

  if(printOnly){
    message('MODELS:')
    print(models)
    message('RC Equations:')
    print(RC_Equations)
  } else{
    return(list(models=models,
                RC_Equations=RC_Equations))
  }
}

#' RatingShiftHappens catalog of recessions model
#'
#' Available recession equations and their estimation
#'
#' @return list with the following fields depending on the model chosen:
#' \enumerate{
#'   \item Equation: recession equation
#'   \item funk: recession estimation
#' }
#' @examples
#' catalog <- GetCatalog_Recession()
#' # Example for the first model in the catalog
#' catalog[[1]]$Equation
#' catalog[[1]]$funk
#' @export
GetCatalog_Recession<-function(){
  # available models: update manually!
  M3=list(Equation=Recession_M3_Equation,
          funk=Estimation_Recession_M3)

  BR1=list(Equation=Recession_B1_Equation,
           funk=Estimation_Recession_BR1)

  return(list(M3=M3,
              BR1=BR1))
}

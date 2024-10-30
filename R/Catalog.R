#' RatingShiftHappens catalog
#'
#' Rating curve and recession models with their equation available
#'
#' @param printOnly Logical, should the catalog be returned or only printed?
#' @return If \code{printOnly==FALSE}, a list with the following fields:
#' \enumerate{
#'   \item models: available models
#'   \item RC Equations: available RC Equations
#' }
#' @examples
#' catalog <- GetCatalog()
#' GetCatalog(printOnly=TRUE)
#' @export
GetCatalog<-function(printOnly=FALSE){
  # available models
  models=c('fitRC_loess',
           'fitRC_BaRatinBAC',
           'fitRC_BaRatinKAC',
           'fitRC_exponential',
           'fitRC_LinearRegression',
           'fitRecession_M3',
           'fitRecession_BR1',
           'fitRecession_BR2')

  RC_Equations=c('Loess_Equation',
              'BaRatinBAC_Equation',
              'BaRatinKAC_Equation',
              'Exponential_Equation',
              'LinearRegression_Equation')

  Recession_Equations=c('Recession_M3_Equation',
                        'Recession_B2_Equation')
  if(printOnly){
    message('MODELS:')
    print(models)
    message('RC Equations:')
    print(RC_Equations)
    message('Recession Equations:')
    print(Recession_Equations)
  } else{
    return(list(models=models,
                RC_Equations=RC_Equations,
                Recession_Equations=Recession_Equations))
  }
}

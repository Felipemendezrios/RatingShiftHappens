#' RatingShiftHappens catalog
#'
#' Rating curve and recession models with their equation available
#'
#' @param printOnly Logical, should the catalog be returned or only printed?
#' @return If \code{printOnly==FALSE}, a list with the following fields:
#' \enumerate{
#'   \item models: available models
#'   \item Equations: available equations
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

  Equations=c('Loess_Equation',
              'BaRatinBAC_Equation',
              'BaRatinKAC_Equation',
              'Exponential_Equation',
              'LinearRegression_Equation',
              'Recession_M3_Equation',
              'Recession_BR1_Equation',
              'Recession_B2_Equation')
  if(printOnly){
    message('MODELS:')
    print(models)
    message('Equations:')
    print(Equations)
  } else{
    return(list(models=models,
                Equations=Equations))
  }
}

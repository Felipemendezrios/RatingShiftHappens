#' RatingShiftHappens catalog
#'
#' Rating curve models available
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
           'fitRC_BaRatin',
           'fitRC_exponential',
           'fitRC_LinearRegression')

  Equations=c('Loess_Equation',
              'BaRatin_Equation',
              'Exponential_Equation',
              'LinearRegression_Equation')
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

#' RatingShiftHappens catalog
#'
#' Rating curve models available
#'
#' @param printOnly Logical, should the catalog be returned or only printed?
#' @return If \code{printOnly==FALSE}, a list with the following fields:
#' \describe{
#'   \item{models}{available models}
#' }
#' @examples
#' catalog <- GetCatalog()
#' GetCatalog(printOnly=TRUE)
#' @export
GetCatalog<-function(printOnly=FALSE){
  # available models
  models=c('fitRC_loess',
           'fitRC_BaRatin')
  if(printOnly){
    message('MODELS:')
    print(models)
  } else{
    return(list(models=models))
  }
}

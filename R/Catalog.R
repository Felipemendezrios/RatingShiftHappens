#' Catalog of rating curves
#'
#' Rating curve and models with their equation available
#'
#' @param printOnly Logical, should the catalog be returned or only printed?
#' @return If \code{printOnly==FALSE}, a list with the following fields:
#' \enumerate{
#'   \item models: available models
#'   \item equations: available RC equations
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

  equations=c('Equation_LOESS',
                 'BaRatinBAC_Equation',
                 'BaRatinKAC_Equation',
                 'Exponential_Equation',
                 'LinearRegression_Equation')

  if(printOnly){
    message('MODELS:')
    print(models)
    message('RC Equations:')
    print(equations)
  } else{
    return(list(models=models,
                equations=equations))
  }
}

#' Catalog of recession models
#'
#' Available recession equations and their estimation
#'
#' @param printOnlyEquations Logical, should the catalog be returned or only printed?
#'
#' @return If \code{printOnlyEquations==FALSE}, a list with the following fields depending on the model chosen:
#' \enumerate{
#'   \item Equation: recession equation
#'   \item funk: recession estimation
#' }
#' @examples
#' # Print only equations
#' GetCatalog_Recession(printOnlyEquations=TRUE)
#'
#' names(GetCatalog_Recession())
#'
#' catalog <- GetCatalog_Recession()
#' # Example for the first model in the catalog
#' catalog[[1]]$Equation
#' catalog[[1]]$funk
#' @export
GetCatalog_Recession<-function(printOnlyEquations=FALSE){

  if(printOnlyEquations!=TRUE){
    # available models: update manually!
    M3=list(Equation=Recession_M3_Equation,
            funk=Estimation_Recession_M3)
    BR1=list(Equation=Recession_B1_Equation,
             funk=Estimation_Recession_BR1)
    return(list(M3=M3,
                BR1=BR1))
  }else{
    message('Available recession equations:')
    print(list(M3=Recession_M3_Equation(),
               BR1=Recession_B1_Equation()))
  }
}

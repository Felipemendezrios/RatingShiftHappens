#' Model 1:
#'
#' Recession model with two exponential and asymptotic: \deqn{\alpha_{1_k} \cdot \exp(-\lambda_1 \cdot t) + \alpha_{2_k} \cdot \exp(-\lambda_2 \cdot t) + \beta_k}
#'
#' @return string, model equation
#' @details
#' Parameter's description:
#'
#' \enumerate{
#'   \item{ t real value, time in days}
#'   \item{ alpha_1_k real value, scaling coefficient in meters}
#'   \item{ lambda_1 real value, fast runoff in 1/days}
#'   \item{ alpha_2_k real value, scaling coefficient in meters}
#'   \item{ lambda_2 real value, slow emptying of the aquifer in 1/days}
#'   \item{ beta_k real value, asymptotic stage in meters}
#' }
#' @export
Recession_M3_Equation <- function(){
  return('alpha_1_k*exp(-lambda_1*t)+alpha_2_k*exp(-lambda_2*t)+beta_k')
}

#' Model 2:
#'
#' Recession with one exponential and four parameters: \deqn{\alpha \cdot \exp(-\lambda \cdot t^c) + \beta_k}
#'
#' @return string, model equation
#' @details
#' Parameter's description:
#' \enumerate{
#'   \item{t real value, time in days}
#'   \item{lambda_1 real value, decay rate constant coefficient 1/(days^c)}
#'   \item{beta_k real value, asymptotic stage in meters}
#'   \item{alpha_1 real value, scaling coefficient in meters}
#'   \item{c real value, recession exponent in meters}
#' }
#' @export
Recession_B2_Equation <- function(){
  return('alpha*exp(-lambda*t^c)+beta_k')
}

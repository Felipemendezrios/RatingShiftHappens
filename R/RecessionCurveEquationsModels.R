#' Model 1:
#'
#' Recession model with two exponential and asymptotic: \deqn{\alpha_{1_k} \cdot \exp(-\lambda_1 \cdot t) + \alpha_{2_k} \cdot \exp(-\lambda_2 \cdot t) + \beta_k}
#'
#' @return Parameter's description:
#'
#' \enumerate{
#'    \item{t real value, time}
#'    \item{alpha_1_k real value, initial stage of the first exponential function}
#'    \item{lambda_1 real value, fast runoff}
#'    \item{alpha_2_k real value, initial stage of the second exponential function}
#'    \item{lambda_2 real value, slow emptying of the aquifer}
#'    \item{beta_k real value, asymptotic stage}
#' }
#' @export
Recession_M3_Equation <- function(){
  return('alpha_1_k * exp(-lambda_1*t) + alpha_2_k * exp(-lambda_2*t) + beta_k')
}

#' Model 2:
#'
#' Recession with one exponential and four parameters: \deqn{\alpha \cdot \exp(-\lambda \cdot t^c) + \beta_k}
#'
#' @return Parameter's description:
#' \enumerate{
#'   \item{t real value, time}
#'   \item{lambda_1 real value, first parameter of the exponential function}
#'   \item{beta_k real value, asymptotic stage}
#'   \item{alpha_1 real value, initial stage of the exponential function}
#'   \item{c real value, recession exponent}
#' }
#' @export
Recession_B2_Equation <- function(){
  return('alpha*exp(-lambda*t^c)+beta_k')
}

#' Recession model with two exponential and asymptotic
#'
#' @param t real value, time
#' @param alpha_1_k real value, initial stage of the first exponential function
#' @param lambda_1 real value, fast runoff
#' @param alpha_2_k real value, initial stage of the second exponential function
#' @param lambda_2 real value, slow emptying of the aquifer
#' @param beta_k real value, asymptotic stage
#'
#' @export
Recession_M3_Equation <- function(t,alpha_1_k,lambda_1,alpha_2_k,lambda_2,beta_k){
  alpha_1_k * exp(-lambda_1*t) + alpha_2_k * exp(-lambda_2*t) + beta_k
}

#' Recession with one exponential and five parameters
#'
#' @param t real value, time
#' @param lambda_1 real value, first parameter of the exponential function
#' @param lambda_2 real value, second parameter of the exponential function
#' @param lambda_3 real value, third parameter of the exponential function
#' @param beta_k real value, asymptotic stage
#' @param alpha_1 real value, initial stage of the exponential function
#'
#' @export
Recession_BR1_Equation <- function(t,alpha_1,lambda_1,lambda_2,lambda_3,beta_k){
  alpha_1*exp(-lambda_1*t-lambda_2*t^2-lambda_3*t^0.5)+beta_k
}

#' Recession with one exponential and four parameters
#'
#' @param t real value, time
#' @param lambda_1 real value, first parameter of the exponential function
#' @param beta_k real value, asymptotic stage
#' @param alpha_1 real value, initial stage of the exponential function
#' @param c real value, recession exponent
#'
#' @export
Recession_B2_Equation <- function(t,alpha_1,lambda_1,c,beta_k){
  alpha_1*exp(-lambda_1*t^c)+beta_k
}

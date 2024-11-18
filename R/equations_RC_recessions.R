#' Equation of LOESS regression
#'
#' @return NULL
#' @export
Loess_Equation<-function(){
  return(NULL)
}

#' Equation of Linear interpolation \deqn{Q(h)= a \cdot H + b}
#'
#' @param H real value, stage
#' @param a real value, slope
#' @param b real value, intercept
#'
#' @return real value, discharge
#' @export
LinearRegression_Equation <- function(H,a,b){
  a*H+b
}

#' Equation of exponential regression
#'
#' \deqn{Q(h)=a \cdot \exp(b \cdot H)}
#'
#' @param H real value, stage
#' @param a real value, initial state
#' @param b real value, growth constant
#'
#' @return real value, discharge
#' @export
Exponential_Equation <- function(H,a,b){
  a*exp(b*H)
}

#' Equation of BaRatin model (b-a-c)
#'
#' \deqn{Q(h) = a \cdot (h-b)^{c} \quad \text{for } (h>k) \quad (\text{and } Q=0 \quad \text{if } h \leq b)}
#'
#' @param H real value, stage
#' @param a real value, coefficient describing geometry characteristics
#' @param b real value, offset
#' @param c real value, exponent
#'
#' @return real value, discharge
#' @export
BaRatinBAC_Equation <- function(H,a,b,c){
  a*(H-b)^c
}

#' Equation of BaRatin model (k-a-c)
#'
#' \deqn{Q(h) = a \cdot (h-b)^{c} \quad \text{for } (h>k) \quad (\text{and } Q=0 \quad \text{if } h \leq b)}
#'
#' @param H real value, stage
#' @param a real value, coefficient describing geometry characteristics
#' @param b real value, offset
#' @param c real value, exponent
#'
#' @return real value, discharge
#' @export
BaRatinKAC_Equation <- function(H,a,b,c){
  a*(H-b)^c
}

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
#' Recession with one exponential and four parameters: \deqn{\alpha_k \cdot \exp(-\lambda \cdot t^c) + \beta_k}
#'
#' @return string, model equation
#' @details
#' Parameter's description:
#' \enumerate{
#'   \item{t real value, time in days}
#'   \item{lambda real value, decay rate constant coefficient 1/(days^c)}
#'   \item{beta_k real value, asymptotic stage in meters}
#'   \item{alpha_k real value, scaling coefficient in meters}
#'   \item{c real value, recession exponent in meters}
#' }
#' @export
Recession_B1_Equation <- function(){
  return('alpha_k*exp(-lambda*t^c)+beta_k')
}


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
#' @return equation, \deqn{Q(h)= a \cdot H + b}
#' @export
LinearRegression_Equation <- function(H,a,b){
  a*H+b
}

#' Equation of exponential regression \deqn{Q(h)=a \cdot \exp(b \cdot H)}
#'
#' @param H real value, stage
#' @param a real value, initial state
#' @param b real value, growth constant
#'
#' @return equation, \deqn{Q(h)=a \cdot \exp(b \cdot H)}
#' @export
Exponential_Equation <- function(H,a,b){
  a*exp(b*H)
}

#' Equation of BaRatin model (b-a-c)
#'
#' @param H real value, stage
#' @param a real value, coefficient describing geometry characteristics
#' @param b real value, offset
#' @param c real value, exponent
#'
#' @return equation, \deqn{Q(h)=a \cdot (H - b)^ c}
#' @export
BaRatinBAC_Equation <- function(H,a,b,c){
  a*(H-b)^c
}

#' Equation of BaRatin model (k-a-c)
#'
#' @param H real value, stage
#' @param a real value, coefficient describing geometry characteristics
#' @param b real value, offset
#' @param c real value, exponent
#'
#' @return equation, \deqn{Q(h)=a \cdot (H - b)^ c}
#' @export
BaRatinKAC_Equation <- function(H,a,b,c){
  a*(H-b)^c
}


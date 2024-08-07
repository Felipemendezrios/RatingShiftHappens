% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RatingCurveEquationsModels.R
\name{BaRatinKAC_Equation}
\alias{BaRatinKAC_Equation}
\title{Equation of BaRatin model (k-a-c)}
\usage{
BaRatinKAC_Equation(H, a, b, c)
}
\arguments{
\item{H}{real value, stage}

\item{a}{real value, coefficient describing geometry characteristics}

\item{b}{real value, offset}

\item{c}{real value, exponent}
}
\value{
equation, \deqn{Q(h)=a \cdot (H - b)^ c}
}
\description{
Equation of BaRatin model (k-a-c)
}

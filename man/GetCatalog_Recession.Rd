% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catalog.R
\name{GetCatalog_Recession}
\alias{GetCatalog_Recession}
\title{Catalog of recession models}
\usage{
GetCatalog_Recession(printOnlyEquations = FALSE)
}
\arguments{
\item{printOnlyEquations}{Logical, should the catalog be returned or only printed?}
}
\value{
If \code{printOnlyEquations==FALSE}, a list with the following fields depending on the model chosen:
\enumerate{
\item Equation: recession equation
\item funk: recession estimation
}
}
\description{
Available recession equations and their estimation
}
\examples{
# Print only equations
GetCatalog_Recession(printOnlyEquations=TRUE)

names(GetCatalog_Recession())

catalog <- GetCatalog_Recession()
# Example for the first model in the catalog
catalog[[1]]$Equation
catalog[[1]]$funk
}

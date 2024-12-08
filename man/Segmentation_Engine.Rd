% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Segmentation.R
\name{Segmentation_Engine}
\alias{Segmentation_Engine}
\title{Segmentation engine}
\source{
\url{https://theses.hal.science/tel-03211343}

\url{https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2020WR028607}
}
\usage{
Segmentation_Engine(
  obs,
  time = 1:length(obs),
  u = 0 * obs,
  nS = 2,
  nMin = 1,
  nCycles = 100,
  burn = 0.5,
  nSlim = max(nCycles/10, 1),
  temp.folder = file.path(tempdir(), "BaM"),
  mu_prior = list(NULL)
)
}
\arguments{
\item{obs}{real vector, observations}

\item{time}{vector, time in POSIXct, string or numeric format}

\item{u}{real vector, uncertainty in observations (as a standard deviation)}

\item{nS}{integer, number of segments}

\item{nMin}{integer, minimum number of observations by segment}

\item{nCycles}{integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles}

\item{burn}{real between 0 (included) and 1 (excluded), MCMC burning factor}

\item{nSlim}{integer, MCMC slim step}

\item{temp.folder}{directory, temporary directory to write computations}

\item{mu_prior}{list, object describing prior knowledge about residual between the rating curve and observation if user-defined (see details)}
}
\value{
list with the following components:
\enumerate{
\item summary: list, summarize the information to present to the user
\itemize{
\item data: data frame, all data with their respective periods after segmentation
\item shift: data frame, all detected shift time in numeric or POSIXct format in UTC
}
\item plot : list, data formatted to use as input for some plot functions
\itemize{
\item density.tau: data frame, a table with three columns. The first column indicates the specific shift being analyzed.
The second column contains the values assessed during the MCMC exploration. The last columns shows the probability density associated
with each tested value
\item density.inc.tau: data frame, all information about the 95\% credibility interval and the Maximum a posterior (MAP) estimation
for shift times with their associated probability densities
}
\item tau: real vector, estimated shift times in numeric or POSIXct format in UTC
\item segments: list, segment maximum a posterior (MAP) value indexed by the list number
\item mcmc: data frame, MCMC simulation
\item data.p: list, separate and assign information by identified stable period indexed by the list number
\item DIC: real, DIC estimation
\item origin.date.p: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
}
}
\description{
Segmentation procedure for a \strong{known} given number of segments
}
\details{
The residuals between the gaugings and the RC are defined as follows:
\deqn{r_i = \tilde{Q_i} - \hat{Q_i} \quad i = 1, \dots, N}
where \eqn{\tilde{Q_i}}  is the gauged discharge, \eqn{\hat{Q_i}} is the RC-estimated discharge and N is the number of gaugings.
Non-informative prior knowledge has been provided by default. If the user wants to modify this, it is advised to keep
in mind that residual must tend towards zero to obtain satisfactory results. Besides, a same prior will be assigned for all calculation.
Additionally, please ensure that the prior has been created using the \code{parameter} function from \code{RBaM} package;
otherwise, an error message will appear as shown below:
\itemize{
\item Error in mu_list$name : the $ operator is invalid for atomic vectors
}

More information about prior knowledge on mu parameter available in \code{sources}.
}
\examples{
# Run segmentation engine function at two segments
# for data set : RhoneRiverAMAX (further details on ?RhoneRiverAMAX)
res=Segmentation_Engine(obs=RhoneRiverAMAX$H,
                        time=RhoneRiverAMAX$Year,
                        u=RhoneRiverAMAX$uH,nS=2,
                        mu_prior = list(RBaM::parameter(name=paste0('mu'),
                                        init=c(0),
                                        prior.dist = 'Gaussian',
                                        prior.par = c(0,50))))
# Data information
knitr::kable(head(res$summary$data),
             align = 'c',row.names = FALSE)
# Shift information
knitr::kable(head(res$summary$shift),
             align = 'c',row.names = FALSE)
# Plot segmentation
PlotSegmentation(res$summary,
                 res$plot)
}

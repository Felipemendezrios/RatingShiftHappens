% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelAndSegment.R
\name{recursive.ModelAndSegmentation}
\alias{recursive.ModelAndSegmentation}
\title{Recursive modelling  and segmentation}
\usage{
recursive.ModelAndSegmentation(
  H,
  Q,
  time = 1:length(H),
  uQ = 0 * Q,
  uH = 0 * H,
  nSmax = 2,
  nMin = 1,
  nCycles = 100,
  burn = 0.5,
  nSlim = max(nCycles/10, 1),
  temp.folder = file.path(tempdir(), "BaM"),
  funk = fitRC_exponential,
  ...
)
}
\arguments{
\item{H}{real vector, stage}

\item{Q}{real vector, discharge}

\item{time}{vector, time in POSIXct, string or numeric format}

\item{uQ}{real vector, uncertainty in discharge (as a standard deviation)}

\item{uH}{real vector, uncertainty in stage record in meters (as a standard deviation)}

\item{nSmax}{integer, maximum number of segments to assess}

\item{nMin}{integer, minimum number of observations by segment}

\item{nCycles}{integer, number of MCMC adaptation cycles. Total number of simulations equal to 100*nCycles.}

\item{burn}{real between 0 (included) and 1 (excluded), MCMC burning factor}

\item{nSlim}{integer, MCMC slim step}

\item{temp.folder}{directory, temporary directory to write computations}

\item{funk}{the function for estimating the rating curve to be applied: see ‘Details’}

\item{...}{optional arguments to funk}
}
\value{
List with the following components :
\enumerate{
\item summary: list, summarize the information to present to the user
\itemize{
\item data: data frame, all data of (H,Q and uQ) with their respective periods after segmentation
\item shift: data frame, all detected shift time
\item param.equation: data frame, parameters estimation
}
\item res: list, provide all the information of the periods from tree structure
\itemize{
\item tau: real vector, estimated shift times
\item segments: list, segment maximum a posterior (MAP) value indexed by the list number
\item mcmc: data frame, MCMC simulation
\item data.p: list, separate and assign information by identified stable period indexed by the list number
\item DIC: real, DIC estimation
\item nS: integer, optimal number of segments following DIC criterion
}
\item tree: data frame, provide tree structure
}
}
\description{
Recursive procedure for an \strong{unknown} number of segments, updating the rating curve at each iteration before segmentation
}
\details{
Some functions for estimating the rating curve are available in this package.
Use \code{GetCatalog()$models} to discover the supported models. More information in \code{?GetCatalog()}.
}
\examples{
# Define fit model for rating curve
fit=fitRC_exponential
equation=Exponential_Equation

# Apply recursive model and segmentation
results=recursive.ModelAndSegmentation(H=ArdecheRiverMeyrasGaugings$H,
                                       Q=ArdecheRiverMeyrasGaugings$Q,
                                       time=ArdecheRiverMeyrasGaugings$Date,
                                       uQ=ArdecheRiverMeyrasGaugings$uQ,
                                       nSmax=3,nMin=2,funk=fit)

# Data information
knitr::kable(head(results$summary$data),
             align = 'c',row.names = FALSE)

# Shift information
knitr::kable(head(results$summary$shift),
             align = 'c',row.names = FALSE)

# Parameters estimation of the rating curve
results$summary$param.equation

# Have a look at recursion tree
results$tree

# Terminal nodes
terminal = results$tree$indx[which(results$tree$nS==1)]
terminal

# Visualize tree structure
plotTree(results$tree)

# See the arguments of the specified fit model for the rating curve
args(equation)

# See parameters estimates for each rating curve
results$summary$param.equation

# It must be pass as input data a and b as parameters of the specified fit model in this case
a=results$summary$param.equation$a
b=results$summary$param.equation$b

# Plot the rating curve after segmentation following a regression exponential
plotRC_ModelAndSegmentation(summary=results$summary,
                            equation=equation,
                            a=a,
                            b=b)

# Plot the rating curves after segmentation with zoom user-defined
plotRC_ModelAndSegmentation(summary=results$summary,
                            equation = equation,
                            autoscale = FALSE,
                            Hmin_user = 1,
                            Hmax_user = 2,
                            H_step_discretization = 0.01,
                            a=a,b=b)

# Plot the rating curves after segmentation in log scale
plotRC_ModelAndSegmentation(summary=results$summary,
                            logscale=TRUE,
                            equation = equation,
                            a=a,
                            b=b)

# Plot the rating curves after segmentation in log scale with zoom
plotRC_ModelAndSegmentation(summary=results$summary,
                            a=a,
                            b=b,
                            logscale=TRUE,
                            equation = equation,
                            autoscale = FALSE,
                            Hmin_user = 0.5,
                            Hmax_user = 2,
                            H_step_discretization = 0.01)

# Plot shift times in stage record
plot_H_ModelAndSegmentation(summary=results$summary,
                            plot_summary=results$plot)

# Plot shift times in discharge observations
plot_Q_ModelAndSegmentation(summary=results$summary,
                            plot_summary=results$plot)

# Plot residual
plotResidual_ModelAndSegmentation(summary=results$summary,
                                  plot_summary=results$plot)


# example with BaRatin method (k-a-c)
fit=fitRC_BaRatinKAC

# Hydraulic matrix control is also required.
# The `control_matrix_builder` was developed to help the user to create this matrix.
# In this case, 3 hydraulic control has been set

controlMatrix=matrix(c(1,0,0,0,1,1,0,0,1),ncol=3,nrow=3)

# Prior information about input data is required.
# The `prior_infor_param_builder` function was developed to help the user to create these objects.
#' # Prior information for Ardeche River at Meyras

a1=RBaM::parameter(name='a1',init=14.17,prior.dist='Gaussian',prior.par=c(14.17,3.65))
k1=RBaM::parameter(name='k1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.6,0.5))
c1=RBaM::parameter(name='c1',init=1.5,prior.dist='Gaussian',prior.par=c(1.5,0.025))
a2=RBaM::parameter(name='a2',init=26.5,prior.dist='Gaussian',prior.par=c(26.5,8.4))
k2=RBaM::parameter(name='k2',init=0,prior.dist='Gaussian',prior.par=c(0,0.5))
c2=RBaM::parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
a3=RBaM::parameter(name='a3',init=31.82,prior.dist='Gaussian',prior.par=c(31.8,10.9))
k3=RBaM::parameter(name='k3',init=1.2,prior.dist='Gaussian',prior.par=c(1.2,0.2))
c3=RBaM::parameter(name='c3',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))

a.object=list(a1,a2,a3)
k.object=list(k1,k2,k3)
c.object=list(c1,c2,c3)

resultsBaRatin=recursive.ModelAndSegmentation(H=ArdecheRiverMeyrasGaugings$H,
                                              Q=ArdecheRiverMeyrasGaugings$Q,
                                              time=ArdecheRiverMeyrasGaugings$Date,
                                              uQ=ArdecheRiverMeyrasGaugings$uQ,
                                              nSmax=3,
                                              nMin=2,
                                              funk=fit,
                                              a.object=a.object,
                                              k.object=k.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
                                              )

# Visualize tree structure
plotTree(resultsBaRatin$tree)

# Terminal nodes
terminal = resultsBaRatin$tree$indx[which(resultsBaRatin$tree$nS==1)]
terminal

# Plot the rating curves after using BaRatin method. Function specially created to this method
plotRCPrediction(Hgrid=data.frame(seq(-1,2,by=0.01)),
                 autoscale=FALSE,
                 temp.folder=file.path(tempdir(),'BaM'),
                 CalibrationData='CalibrationData.txt',
                 allnodes=FALSE,
                 nodes=terminal)

# Plot shift times in stage record
plot_H_ModelAndSegmentation(summary=resultsBaRatin$summary,
                            plot_summary=resultsBaRatin$plot)

# Plot shift times in discharge observations
plot_Q_ModelAndSegmentation(summary=resultsBaRatin$summary,
                            plot_summary=resultsBaRatin$plot)

# Plot residual
plotResidual_ModelAndSegmentation(summary=resultsBaRatin$summary,
                                  plot_summary=resultsBaRatin$plot)

# Plot all gaugings after segmentation in a plot H-Q
plotGaugingsSegmented(summary=resultsBaRatin$summary)
}

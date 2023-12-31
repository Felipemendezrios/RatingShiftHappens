% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plots.R
\name{plotTree}
\alias{plotTree}
\title{Plot recursive segmentation tree}
\usage{
plotTree(tree)
}
\arguments{
\item{tree}{data frame, tree resulting from the call of recursive.segmentation}
}
\value{
a ggplot
}
\description{
Plot the tree resulting from the recursive segmentation procedure
}
\examples{

# Apply recursive segmentation
results=recursive.segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH)

# plot recursion tree
plotTree(results$tree)
}

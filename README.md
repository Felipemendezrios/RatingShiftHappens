Rating Shift Happens
================
Felipe MENDEZ and Benjamin RENARD (INRAE, RiverLy and RECOVER).
December 2023

## Introduction

The goal of `RatingShiftHappens` package is to create a tools package
for detecting, visualizing and estimating rating shifts. This package
was derived from [BayDERS](https://github.com/MatteoDarienzo/BayDERS)
developed by Darienzo in 2021.

This documentation provides the description of several functions
available to the segmentation process.

Three fundamental functions are available so far :

1.  segmentation.engine
2.  segmentation
3.  recursive.segmentation

## Installation

You can install the development version of
[RatingShiftHappens](https://github.com/Felipemendezrios/RatingShiftHappens)
from GitHub using the following command. Please note that the
[RBaM](https://github.com/BaM-tools/RBaM) package is also required to
run this package.

``` r
# Before first use, install Rating Shift Happens and RBaM packages ones and for all, following these commands: 

# devtools::install_github("Felipemendezrios/RatingShiftHappens")
# devtools::install_github('BaM-tools/RBaM') 

library(RatingShiftHappens)
```

Functions will be explained more precisely below along with a example.

## Segmentation procedure for a *known* given number of segments

This basic example demonstrates the segmentation of annual maximum
stages (H, m) for the Rhone River at Beaucaire, France, along with the
associated uncertainties expressed as standard deviations (uH), divided
into two groups. More information about the data set, please refer to
the documentation available in `?RhoneRiver`.

``` r
 # Run segmentation engine function at two segments
 res=segmentation.engine(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH,nS=2)

 # Data information
 knitr::kable(head(res$summary$data),
              align = 'c',row.names = F)
```

| time | obs  |   u    | I95\_lower | I95\_upper | period |
|:----:|:----:|:------:|:----------:|:----------:|:------:|
| 1816 | 5.69 | 0.4500 |  4.808016  |  6.571984  |   1    |
| 1817 | 4.56 | 0.4075 |  3.761315  |  5.358685  |   1    |
| 1818 | 4.72 | 0.4300 |  3.877216  |  5.562785  |   1    |
| 1819 | 5.08 | 0.4125 |  4.271515  |  5.888485  |   1    |
| 1820 | 5.12 | 0.4925 |  4.154718  |  6.085282  |   1    |
| 1821 | 5.40 | 0.3900 |  4.635614  |  6.164386  |   1    |

``` r
 # Shift information
 knitr::kable(head(res$summary$shift),
              align = 'c',row.names = F)
```

|   tau   | I95\_lower | I95\_upper |
|:-------:|:----------:|:----------:|
| 1969.12 |  1967.09   |  1971.35   |

``` r
 # Plot segmentation
 plotSegmentation(res$summary)
```

<img src="man/readme/README-segmentation.engine-1.png" width="100%" />

## For more advanced details :

MCMC sampling demonstrate all combinations of parameters estimated.

``` r
knitr::kable(head(res$mcmc),align = 'c')
```

|   mu1   |   mu2   |  tau1   | structural\_sd | LogPost  |
|:-------:|:-------:|:-------:|:--------------:|:--------:|
| 5.20155 | 7.83626 | 1967.03 |    1.07363     | -309.191 |
| 5.33280 | 7.87847 | 1970.18 |    1.10196     | -307.813 |
| 5.33280 | 7.82841 | 1970.42 |    1.00044     | -309.262 |
| 5.34227 | 8.03573 | 1967.69 |    1.09610     | -309.200 |
| 5.52388 | 7.24338 | 1967.76 |    1.14877     | -314.065 |
| 5.24600 | 7.32321 | 1969.19 |    1.14107     | -311.731 |

A few functions are provided with the `RBaM` package to explore MCMC
samples.

``` r
  # Trace plot for each parameter, useful to assess convergence.
  plots=RBaM::tracePlot(res$mcmc)
  gridExtra::grid.arrange(grobs=plots,ncol=3)
```

<img src="man/readme/README-unnamed-chunk-4-1.png" width="100%" />

``` r
  # Density plot for each parameter
  plots=RBaM::densityPlot(res$mcmc)
  gridExtra::grid.arrange(grobs=plots,ncol=3)
```

<img src="man/readme/README-unnamed-chunk-4-2.png" width="100%" />

In this example, the focus will be on exploring the uncertainty
associated with the first shift time.

``` r
 Shift=data.frame(time=res$mcmc$tau1)
 ggplot2::ggplot(Shift,ggplot2::aes(x=time))+
   ggplot2::geom_histogram(ggplot2::aes(y=..density..),col=1,fill='white',bins=80)+
   ggplot2::labs(title='Histogram with density of first shift')+
   ggplot2::theme_bw()+
   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
   ggplot2::geom_density(col=4,lwd=1,fill=4,alpha=0.25)
```

<img src="man/readme/README-unnamed-chunk-5-1.png" width="100%" />

## Segmentation procedure for an *unknown* number of segments

This is a basic example which shows you how to segment the same dataset
with an **unknown** number of segments :

``` r
 # Run segmentation engine function at two segments
 res=segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH,nSmax=3)

 # Get lower DIC value and optimal number of segments (to define optimal solution)
 DIC.df = data.frame(nS=c(1:3),DIC=c(res$results[[1]]$DIC,res$results[[2]]$DIC,res$results[[3]]$DIC))
 nSopt=res$nS
 
 ggplot2::ggplot(DIC.df,ggplot2::aes(x=nS,y=DIC,col=factor(nS)))+
   ggplot2::geom_point(size=3,show.legend = F)+
   ggplot2::geom_segment(ggplot2::aes(x=nSopt,y=min(DIC)*1.03,xend=nSopt,yend=min(DIC)*1.005),
                         arrow=ggplot2::arrow(length=ggplot2::unit(0.5,'cm')),
                         color='BLACK',lwd=1, show.legend = F)+
   ggplot2::theme_bw()
```

<img src="man/readme/README-segmentation-1.png" width="100%" />

``` r
 # Data information
 knitr::kable(head(res$results[[nSopt]]$summary$data),
              align = 'c',row.names = F)
```

| time | obs  |   u    | I95\_lower | I95\_upper | period |
|:----:|:----:|:------:|:----------:|:----------:|:------:|
| 1816 | 5.69 | 0.4500 |  4.808016  |  6.571984  |   1    |
| 1817 | 4.56 | 0.4075 |  3.761315  |  5.358685  |   1    |
| 1818 | 4.72 | 0.4300 |  3.877216  |  5.562785  |   1    |
| 1819 | 5.08 | 0.4125 |  4.271515  |  5.888485  |   1    |
| 1820 | 5.12 | 0.4925 |  4.154718  |  6.085282  |   1    |
| 1821 | 5.40 | 0.3900 |  4.635614  |  6.164386  |   1    |

``` r
 # Shift information
 knitr::kable(head(res$results[[nSopt]]$summary$shift),
              align = 'c',row.names = F)
```

|   tau   | I95\_lower | I95\_upper |
|:-------:|:----------:|:----------:|
| 1969.12 |  1967.09   |  1971.35   |

``` r
 # Plot segmentation
 plotSegmentation(res$summary)
```

<img src="man/readme/README-segmentation-2.png" width="100%" />

## Recursive segmentation procedure for an *unknown* number of segments

This is a basic example which shows you how to segment the data set with
an **unknown** number of segments using a recursive process:

``` r
 # Apply recursive segmentation
 results=recursive.segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH,nSmax=3)
 
 # Data information
 knitr::kable(head(results$summary$data),
              align = 'c',row.names = F)
```

| time | obs  |   u    | I95\_lower | I95\_upper | period |
|:----:|:----:|:------:|:----------:|:----------:|:------:|
| 1816 | 5.69 | 0.4500 |  4.808016  |  6.571984  |   1    |
| 1817 | 4.56 | 0.4075 |  3.761315  |  5.358685  |   1    |
| 1818 | 4.72 | 0.4300 |  3.877216  |  5.562785  |   1    |
| 1819 | 5.08 | 0.4125 |  4.271515  |  5.888485  |   1    |
| 1820 | 5.12 | 0.4925 |  4.154718  |  6.085282  |   1    |
| 1821 | 5.40 | 0.3900 |  4.635614  |  6.164386  |   1    |

``` r
 # Shift information
 knitr::kable(head(results$summary$shift),
              align = 'c',row.names = F)
```

|   tau   | I95\_lower | I95\_upper |
|:-------:|:----------:|:----------:|
| 1969.12 |  1967.09   |  1971.35   |

``` r
 # Have a look at recursion tree
 results$tree
#>   indx level parent nS
#> 1    1     1      0  2
#> 2    2     2      1  1
#> 3    3     2      1  1

 # Visualize tree structure
 plotTree(results$tree)
```

<img src="man/readme/README-unnamed-chunk-6-1.png" width="100%" />

``` r
 # Plot segmentation
 plotSegmentation(summary=results$summary)
```

<img src="man/readme/README-unnamed-chunk-6-2.png" width="100%" />

#' Plot recursive segmentation tree
#'
#' Plot the tree resulting from the recursive segmentation procedure
#'
#' @param tree data frame, tree resulting from the call of recursive.segmentation
#'
#' @return a ggplot
#'
#' @examples
#'
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH)
#'
#' # plot recursion tree
#' plotTree(results$tree)
#' @export
#' @import ggplot2
plotTree <- function(tree){
  DF=tree
  n=NROW(DF)
  # Add columns for x/y's of nodes and arrows
  DF$y=-1*DF$level+stats::rnorm(n,0,0.2)
  DF$x=0*DF$y
  DF$xstart=DF$xend=DF$ystart=DF$yend=NA*DF$y
  DF$isTerminal=FALSE

  if(n>1){
    for(i in 2:n){
      # count nodes having the same parent as current
      mask=DF$parent==DF$parent[i]
      K=sum(mask)
      # deduce x's
      shifts=(1:K)-mean(1:K)
      DF$x[mask]=DF$x[DF$parent[i]]+shifts
      # arrows
      DF$xstart[i]=DF$x[DF$parent[i]]
      DF$xend[i]=DF$x[i]
      DF$ystart[i]=DF$y[DF$parent[i]]-0.1
      DF$yend[i]=DF$y[i]+0.1
      # is node terminal?
      DF$isTerminal[i]=sum(DF$parent==DF$indx[i])==0
    }
  } else {
    DF$isTerminal[1]=TRUE
  }

  DF$fontface='plain';DF$fontface[DF$isTerminal]='bold'
  # plot
  g=ggplot(DF)+
    geom_segment(aes(x=xstart,y=ystart,xend=xend,yend=yend),
                 arrow=arrow(type='closed',length=unit(0.03, "npc")))+
    geom_label(aes(x,y,label=indx,fill=as.factor(parent),fontface=fontface,size=isTerminal),color='black',alpha=0.8)+
    scale_fill_brewer('Parent node',palette='Set1')+
    scale_size_manual('Terminal node',values=c(4,6))+
    theme_void()
  return(g)
}

#' Plots segmentation
#'
#' Plots  focusing on the observed data and shift time and their estimation
#'
#' @param summary data frame, summary data resulting from any segmentation function without estimation of rating curve
#' @param plot_summary list, plot data resulting from any segmentation function
#'
#' @return  List with the following components :
#' \enumerate{
#'   \item final_plot: ggplot, observed data presenting the shift time and their estimation
#'   \item observation_and_shift: ggplot, observed data indexed by period, and the estimated shift time is indicated vertically
#'   \item shift_time_density: ggplot, density plot of the shift times following MCMC results, with
#'         vertical lines indicating the 95% credibility interval and a red cross representing shift time assignment
#'   }
#' @export
#' @importFrom ggnewscale new_scale_color
#' @import patchwork
plotSegmentation <- function(summary,
                             plot_summary) {

  # Plot Segmentation showing shift time and observed data indexed by period
  data=summary$data
  shift=summary$shift

  # Add some colors to the palette for observations
  colourCount_obs = length(unique(data$period))
  getPalette_obs =  scales::viridis_pal(option='D')

  # Add some colors to the palette for shift
  colourCount_tau = length(unique(shift$tau))
  getPalette_tau_MAP =  scales::viridis_pal(option='H')

  # Plot shift times
  obs_shift_plot=ggplot(data)+
    geom_vline(data=shift,
               aes(xintercept = tau,
                   col=factor(tau)),
               alpha=0.8)+
    geom_point(data=shift,
               aes(x=tau,
                   y=min(data$I95_lower)),
               col='red',
               shape=4,
               size=3)+
    coord_cartesian(ylim = c(min(data$I95_lower),max(data$I95_upper)))+
    labs(col='Global shift time')

  if(is.numeric(shift$tau)){
    obs_shift_plot=obs_shift_plot+
      scale_color_manual(values=getPalette_tau_MAP(colourCount_tau),
                         labels=round(shift$tau,2))
  }else{
    obs_shift_plot=obs_shift_plot+
      scale_color_manual(values=getPalette_tau_MAP(colourCount_tau),
                         labels=round(shift$tau,units='days'))
  }

  # Plot observations by period
  obs_shift_plot=
    obs_shift_plot+
    ggnewscale::new_scale_color()+
    geom_point(aes(x=time,
                   y=obs,
                   col=factor(period)))+
    geom_errorbar(aes(x=time,
                      y=obs,
                      ymin=I95_lower,
                      ymax=I95_upper,
                      col=factor(period)),
                  width=3)+
    labs(y='Observation',
         x='Time',
         col='Period',
         title = 'Segmentation')+
    scale_color_manual(values = getPalette_obs(colourCount_obs))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Plot probability density function of the shift time
  colourCount_shift = length(unique(plot_summary$density.inc.tau$Shift))

  getPalette_tau_density = scales::viridis_pal(option = "H")
  getPalette_tau_trait =  scales::viridis_pal(option='E')

  # Customize labels
  label_shift <- paste0('Shift ',seq(1,colourCount_shift))
  if(length(which(colnames(plot_summary$density.inc.tau)=='id_iteration'))!=0){
    label_wrap <- stats::setNames(paste('Shift detected in this node : ',
                                        unique(plot_summary$density.inc.tau$id_iteration)),
                                  unique(plot_summary$density.inc.tau$id_iteration))
  }

  pdf_shift_plot= ggplot(plot_summary$density.tau,
                         aes(x=Value,
                             y = Density,
                             fill=Shift))+
    geom_area(alpha=0.4,
              position = "identity")+
    geom_segment(data=plot_summary$density.inc.tau,
                 aes(x=tau_lower_inc,
                     y=0 ,
                     xend=tau_lower_inc,
                     yend=density_tau_lower_inc,
                     col=factor(Shift)),
                 alpha=0.7)+
    geom_segment(data=plot_summary$density.inc.tau,
                 aes(x=tau_upper_inc,
                     y=0 ,
                     xend=tau_upper_inc,
                     yend=density_tau_upper_inc,
                     col=factor(Shift)),
                 alpha=0.7)+
    geom_point(data=plot_summary$density.inc.tau,
                 aes(x=taU_MAP,
                     y=density_taU_MAP),
                 col='red',
               shape=4,
               size=2)+
    scale_fill_manual(values=getPalette_tau_density(colourCount_shift),
                      labels=label_shift)+
    scale_color_manual(values=getPalette_tau_trait(colourCount_shift),
                       labels=label_shift)+
    theme_bw()+
    coord_cartesian(xlim = c(min(data$time),max(data$time)))+
    labs(x='Time',
         y='Probability density',
         fill=NULL,
         col=NULL)+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  if(length(which(colnames(plot_summary$density.inc.tau)=='id_iteration'))!=0){
    pdf_shift_plot=pdf_shift_plot+
      facet_wrap(~id_iteration,
                 ncol=1,
                 scales='free_y',
                 labeller = labeller(id_iteration=label_wrap))
  }

  return(list(final_plot=obs_shift_plot/pdf_shift_plot,
              observation_and_shift=obs_shift_plot,
              shift_time_density=pdf_shift_plot))
}

#' Plot Rating curve after segmentation
#'
#' Plot of all the rating curve after using the function `recursive.ModelAndSegmentation`, along with associated uncertainties
#'
#' @param summary data.frame, summary data resulting from the function `recursive.ModelAndSegmentation` saved in `$summary`
#' @param equation the equation to be considered for plotting : see ‘Details’
#' @param ... extra arguments for the specified equation : see 'Details'. Respect the name of the arguments (a, b, ..., c, control matrix) depending of the equation used.
#' @param Hmin_user real value, minimal stage record set by user for plotting in meters
#' @param Hmax_user real value, minimal stage record set by user for plotting in meters
#' @param H_step_discretization real positive non-zero value, time step of the discretization in meters
#' @param autoscale logical, auto scale following data for plotting
#' @param logscale logical, `TRUE`= log scale in discharge axis, `FALSE` = normal scale in discharge axis
#'
#'
#' @return ggplot, rating curve after segmentation
#' @details
#' Some equations for estimating the rating curve are available in this package.
#' Use `GetCatalog()$Equations` to discover the equations supported. More information in `?GetCatalog()`.
#' For the extra information, you can view the arguments through `args()`  by enclosing the equation specified in the input arguments in brackets as shown in the example of `?recursive.ModelAndSegmentation``.
#' Information about parameters values are available from the function `recursive.ModelAndSegmentation` saved in `$summary$param.equation`. Save each parameter estimates separately
#' Please ensure that you enter the same number of rating curve parameters in the extra information as shown in the example `?recursive.ModelAndSegmentation`.
#' @export
plotRC_ModelAndSegmentation=function(summary,
                                     equation=Exponential_Equation,
                                     ...,
                                     Hmin_user = 0,
                                     Hmax_user = 2,
                                     H_step_discretization=0.01,
                                     autoscale=TRUE,
                                     logscale = FALSE){
  # Check discretization
  if(H_step_discretization<=0)stop('The discretization must be a positive non-zero value')

  # Summary data
  summary_data=summary$data
  # summary equation parameters
  param=summary$param.equation
  # summary shift times adding first and last discharge measurement
  summary_shift=rbind(min(summary_data$time),
                      summary$shift,
                      max(summary_data$time))

  # Check negative and or zero values in discharge
  columns_to_check = c('Q_I95_lower','Q_I95_upper','Qsim','Qsim_I95_lower','Qsim_I95_upper')
  if(any(summary_data<=0)){
    if(logscale==TRUE){
      warning(paste0('Rating curve presented in log scale does not accept negative values and zero. Some values have been replace by 1e-6 '))
      DF_remplacement = replace_negatives_or_zero_values(data_frame = summary_data,
                                                         columns = columns_to_check,
                                                         consider_zero = TRUE,
                                                         replace = 1e-6  )
    }else{
      warning(paste0('Rating curve does not accept negative values zero. Some values have been replace by zero'))
      DF_remplacement = replace_negatives_or_zero_values(data_frame = summary_data,
                                                         columns = columns_to_check,
                                                         consider_zero = FALSE,
                                                         replace = 0 )
    }
    columns.match = match(columns_to_check,colnames(summary_data))
    columns <- columns.match[!is.na(columns.match)]
    summary_data[,columns]=DF_remplacement
  }

  # Add some colors to the palette for observations
  colourCount_obs = length(unique(summary_data$period))
  getPalette_obs =  scales::viridis_pal(option='D')

  ### limit of plotting
  if(autoscale==TRUE){
    H_min_plot_limit = min(summary_data$H)
    H_max_plot_limit = max(summary_data$H)
  }else{
    H_min_plot_limit = Hmin_user
    H_max_plot_limit = Hmax_user
  }

  # Check H_step_discretization
  if(H_step_discretization>(H_max_plot_limit-H_min_plot_limit))stop('Discretization is not possible, verify H_step_discretization')
  if(((H_max_plot_limit-H_min_plot_limit)/H_step_discretization<100))warning('Discretization is performed with less than 100 points. Consider to reduce H_step_discretization')

  # Check H limit defined by user
  if(H_max_plot_limit<=H_min_plot_limit)stop('Plot view need to be verify, Hmax_user is lower than Hmin_user. Try autoscale=TRUE')

  # Structural uncertainty associated to the simulation of discharge by rating curve. Parametric uncertainty has not been propagated
  uq_lower=unique(stats::qnorm(0.025)*summary_data$uQ_sim)
  uq_upper=unique(stats::qnorm(0.975)*summary_data$uQ_sim)

  # Estimate the grid each H_step_discretization in meters
  H_grid= seq(H_min_plot_limit, H_max_plot_limit, by = H_step_discretization)

  # Summary table of observations, extract data from summary_data and add period in time format (numeric or date depending on input format)
  summary_obs=data.frame(H_obs=summary_data$H,
                         Q_obs=summary_data$Q,
                         Q_I95_lower=summary_data$Q_I95_lower,
                         Q_I95_upper=summary_data$Q_I95_upper,
                         period=summary_data$period,
                         period_date_begin=0,
                         period_date_end=0)


  Q_sim_grid <- data.frame(matrix(NA, nrow = length(H_grid), ncol = length(param)))

  for (i in seq_along(H_grid)) {
    # Discharge simulation following equation specified in input arguments, discretized in H_grid
    Q_sim_grid [i,] <- equation(H_grid [i],...)
  }
  # Lower uncertainty of discharge simulations
  Q_sim_lower <-  data.frame(sapply(1:ncol(Q_sim_grid), function(col_index) {
    Q_sim_grid[, col_index] + uq_lower[col_index]}))
  # Upper uncertainty of discharge simulations
  Q_sim_upper <-  data.frame(sapply(1:ncol(Q_sim_grid), function(col_index) {
    Q_sim_grid[, col_index] + uq_upper[col_index]}))


  # Summarize results index by period discretized in H_grid if possible if not at stage observed
  period_date_begin_obs <- c()
  period_date_end_obs <- c()
  summary_sim_grid <- c()

  for(i in 1:length(param)){ # scan all rating curve estimates
    if(all(!is.na(Q_sim_grid[,i]))){ # Verify if results can be discretized in H_grid

      summary_sim_grid_p=data.frame(H_sim=H_grid,
                                    Qsim=Q_sim_grid[,i],
                                    Qsim_I95_lower=Q_sim_lower[,i],
                                    Qsim_I95_upper=Q_sim_upper[,i],
                                    period=i,
                                    period_date_begin=summary_shift$tau[i],
                                    period_date_end=summary_shift$tau[i+1])

    }else{

      warning(paste0('Fit rating curve for the period ', i,
                     ' without any parameters. Discharge simulation estimated at observed stage'))
      summary_data_p=summary_data[which(summary_data$period==i),]

      summary_sim_grid_p =  data.frame(H_sim=summary_data_p$H,
                                       Qsim=summary_data_p$Qsim,
                                       Qsim_I95_lower=summary_data_p$Qsim_I95_lower,
                                       Qsim_I95_upper=summary_data_p$Qsim_I95_upper,
                                       period=summary_data_p$period,
                                       period_date_begin=summary_shift$tau[i],
                                       period_date_end=summary_shift$tau[i+1])

    }
    summary_sim_grid=rbind(summary_sim_grid,summary_sim_grid_p)

    # Period in time format for observations
    period_date_begin_obs_p=rep(summary_shift$tau[i],length(which(summary_data$period==i)))
    period_date_begin_obs=rbind(period_date_begin_obs,data.frame(begin_validity=period_date_begin_obs_p))

    period_date_end_obs_p=period_date_end=rep(summary_shift$tau[i+1],length(which(summary_data$period==i)))
    period_date_end_obs=rbind(period_date_end_obs,data.frame(end_validity=period_date_end_obs_p))
  }

  # Check negative and or zero values in discharge simulations
  if(any(summary_sim_grid<=0)){
    if(logscale==TRUE){
      DF_remplacement = replace_negatives_or_zero_values(data_frame = summary_sim_grid,
                                                         columns = columns_to_check,
                                                         consider_zero = TRUE,
                                                         replace = 1e-6   )
    }else{
      DF_remplacement = replace_negatives_or_zero_values(data_frame = summary_sim_grid,
                                                         columns = columns_to_check,
                                                         consider_zero = FALSE,
                                                         replace = 0 )
    }
    columns.match = match(columns_to_check,colnames(summary_sim_grid))
    columns <- columns.match[!is.na(columns.match)]
    summary_sim_grid[,columns]=DF_remplacement
  }

  # Assign period date of observations
  summary_obs$period_date_begin=period_date_begin_obs$begin_validity
  summary_obs$period_date_end=period_date_end_obs$end_validity

  # Plot RC steps :
  min_x <- min(summary_sim_grid$H_sim,summary_obs$H_obs)
  max_x <- max(summary_sim_grid$H_sim,summary_obs$H_obs)
  min_y <- min(summary_sim_grid$Qsim_I95_lower,summary_obs$Q_I95_lower)
  max_y <- max(summary_sim_grid$Qsim_I95_upper,summary_obs$Q_I95_upper)

  # Create an empty ggplot to paste graphs over it
  RC_plot <- ggplot()

  RC_plot = RC_plot+
    # Plot Simulation (MAP and uncertainty)
    geom_ribbon(data=summary_sim_grid,
                aes(x=H_sim,
                    ymin=Qsim_I95_lower,
                    ymax=Qsim_I95_upper ,
                    fill=factor(period)),
                alpha=0.3,
                na.rm = FALSE)+
    geom_line(data=summary_sim_grid,
              aes(x=H_sim,
                  y=Qsim,
                  col=factor(period_date_end)))+
    # Plot observations with uncertainty :
    geom_errorbar(data=summary_obs,
                  aes(x=H_obs,
                      y=Q_obs,
                      ymin=Q_I95_lower,
                      ymax=Q_I95_upper,
                      col=factor(period_date_end)),
                  width=0.1)+
    geom_point(data=summary_obs,
               aes(x=H_obs,
                   y=Q_obs,
                   col=factor(period_date_end)))
  # Customize plot
  RC_plot = RC_plot +
    coord_cartesian(xlim=c(H_min_plot_limit,H_max_plot_limit))+
    scale_x_continuous(breaks = seq(floor(min_x / 0.5) * 0.5,
                                    ceiling(max_x / 0.5) * 0.5,
                                    by=0.1))+
    labs(x='Stage (m)',
       y='Discharge (m3/s)',
       title = 'Rating curves after segmentation')+
    scale_color_manual(values = getPalette_obs(colourCount_obs),
                       name = 'Final validity date')+
    scale_fill_manual(values = getPalette_obs(colourCount_obs),
                      name = 'Period')+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Log scale
  if(logscale==TRUE){

    breaks <- 10^seq(floor(log10(min_y)), ceiling(log10(max_y)), by = 1)

    # Custom label formatting function
    custom_format <- function(x) {
      ifelse(x >= 1,
             as.character(round(x)),
             sapply(x, function(val) {
               formatted <- format(val, scientific = FALSE)
               if (grepl("\\.0+$", formatted)) {
                 formatted <- sub("\\.0+$", "", formatted)
               }
               formatted
             })
      )
    }
  RC_plot= RC_plot+ scale_y_continuous(trans = "log10", breaks = breaks, labels=custom_format)
  }

  return(RC_plot)
}

#' Plot stage time series after segmentation
#'
#' Plots focusing on the observed data and shift times after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#' @param plot_summary list, plot data resulting from any segmentation function
#' @param uH vector or value, uncertainty concerning the observed stage in meters
#'
#'
#' @return  List with the following components :
#' \enumerate{
#'   \item final_plot: ggplot, observed data presenting the shift time and their estimation
#'   \item observation_and_shift: ggplot, observed data indexed by period, and the estimated shift time is indicated vertically,
#'   \item shift_time_density: ggplot, density plot of the shift times following MCMC results, with
#'         vertical lines indicating the 95% credibility interval and a red cross representing shift time assignment
#'   }
#' @export
plot_H_ModelAndSegmentation <- function(summary,
                                        plot_summary,
                                        uH=NA){
  # Adapt summary to use PlotSegmentation function to plot segmentation of stage time series
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$H,
                             I95_lower=summary$data$H-uH,
                             I95_upper=summary$data$H+uH,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                        shift=summary$shift)

  plot_segmentation=plotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary)

  plot_segmentation[[1]][[1]]$labels$y='Stage record (m)'
  plot_segmentation[[2]]$labels$y='Stage record (m)'

  return(plot_segmentation)
}
#' Plot discharge time series after segmentation
#'
#' Plots focusing on the observed data and shift times after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#' @param plot_summary list, plot data resulting from any segmentation function
#' @param uH vector or value, uncertainty concerning the observed stage in meters
#'
#'
#' @return  List with the following components :
#' \enumerate{
#'   \item final_plot: ggplot, observed data presenting the shift time and their estimation
#'   \item observation_and_shift: ggplot, observed data indexed by period, and the estimated shift time is indicated vertically,
#'   \item shift_time_density: ggplot, density plot of the shift times following MCMC results, with
#'         vertical lines indicating the 95% credibility interval and a red cross representing shift time assignment
#'   }
#' @export
plot_Q_ModelAndSegmentation <- function(summary,
                                        plot_summary,
                                        uH=NA){
  # Adapt summary to use PlotSegmentation function to plot segmentation of discharge measurements
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$Q,
                             I95_lower=summary$data$Q_I95_lower,
                             I95_upper=summary$data$Q_I95_upper,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plot_segmentation=plotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary)

  plot_segmentation[[1]][[1]]$labels$y='Discharge (m3/s)'
  plot_segmentation[[2]]$labels$y='Discharge (m3/s)'

  return(plot_segmentation)
}

#' Plot residuals with updating rating curve
#'
#' Plot residuals after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#' @param plot_summary list, plot data resulting from any segmentation function
#'
#' @return ggplot, residual segmentation
#' @export
plotResidual_ModelAndSegmentation <- function(summary,
                                              plot_summary){

  if(is.null(summary$shift))stop('Any shift time detected')

  # Adapt summary to use PlotSegmentation function to plot segmentation of discharge measurements
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$Qres,
                             I95_lower=summary$data$Q_I95_lower-summary$data$Qsim,
                             I95_upper=summary$data$Q_I95_upper-summary$data$Qsim,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plotresidual=plotSegmentation(summary = summary_plot,
                                plot_summary = plot_summary)

  plotresidual[[1]][[1]]$labels$y='Residual (m3/s) : Observed discharge minus simulated discharge'
  plotresidual[[2]]$labels$y='Residual (m3/s) : Observed discharge minus simulated discharge'

  return(plotresidual)
}

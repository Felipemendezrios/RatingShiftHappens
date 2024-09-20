#' Plot recursive segmentation tree
#'
#' Plot the tree resulting from any recursive segmentation procedures
#'
#' @param tree data frame, tree resulting from the call of recursive.segmentation and recursive.ModelAndSegmentation
#'
#' @return a ggplot
#'
#' @examples
#'
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs=CongoRiverBrazzavilleMINAN$Q,
#'                                time=CongoRiverBrazzavilleMINAN$year,
#'                                u=CongoRiverBrazzavilleMINAN$uQ,
#'                                nSmax=2)
#'
#' # plot recursion tree
#' plotTree(results$tree)
#' @export
#' @import ggplot2
#' @importFrom stats rnorm setNames
#' @importFrom scales viridis_pal
plotTree <- function(tree){
  DF=tree
  n=NROW(DF)

  if(n==1)stop('Any shift detected, non tree structure to plot')
  # Add columns for x/y's of nodes and arrows
  DF$y=-1*DF$level
  DF$x=0*DF$y
  DF$xstart=DF$xend=DF$ystart=DF$yend=NA*DF$y
  DF$isTerminal=FALSE

  # Assign x coordinates for all nodes in a same level
  x_ref_plot=0 # central axis fixed at 0 for the first level
  for(i in 2:max(DF$level)){# first level skipped, always x = 0
    indx_table=which(DF$level==i)
    # Get the number of nodes in this level and assign coordinates with middle point at x=0
    half_n <- (length(indx_table) - 1) / 2
    coord_local_node <- seq(-half_n, half_n)

    id_parents=unique(DF$parent[indx_table])

    # if same parent -> reinitialize x offset
    if(length(id_parents)==1){
      x_ref_plot = DF$x[id_parents]
    }
    # Allocate coordinates of nodes at this level adding x_ref_plot estimated of the last level (offset) to move the reference axis
    DF$x[indx_table] = coord_local_node+x_ref_plot

    # Give a weight to each node in function of the number of children. This helps to offset the reference axis to plot new nodes
    # There are more children on the left than on the right, so moving to the left makes more sense for the plot of the next level.
    weight_x_plot =  ifelse(coord_local_node+x_ref_plot < 0, -1, ifelse(coord_local_node+x_ref_plot > 0, 1, 0)) *  DF$nS[indx_table]/max(DF$nS[indx_table])

    x_ref_plot=x_ref_plot+
        (sum(weight_x_plot)) # move the xref_plot to plot nodes in this level

}
  # Assign segments to link nodes
  if(n>1){
    for(i in 2:n){
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
  # Define the colors and remove the level 0
  tree_levels <- unique(DF$parent)[-1]
  tree_colors <-scales::viridis_pal(option='D')(length(tree_levels))

  # Create a named vector of colors for the legend
  color_legend <- stats::setNames(tree_colors, as.character(tree_levels))
  # Get the first node
  DF_initial_node=DF[1,]
  # Get tree structure without first node
  DF_tree_plot=DF[-1,]

  # plot
  g=ggplot(DF_tree_plot)+
    # plot initial node
    geom_point(data=DF_initial_node,
               aes(x,y),
               size=11, col='gray', shape=16, alpha=0.5)+
    # Text for the initial node
    geom_text(data=DF_initial_node,
              aes(x,y,label=indx),
              size=4)+
    # Arrows to link the nodes
    geom_segment(aes(x=xstart,y=ystart,xend=xend,yend=yend),
                 arrow=arrow(type='closed',length=unit(0.015, "npc")))+
    # Plot nodes of the tree structure
    geom_point(aes(x,y,col=factor(parent),shape=factor(isTerminal)),
              size=8,alpha=0.8)+
    geom_text(aes(x,y,label=indx,fontface=fontface),
              size=4)+
    scale_color_manual(name='Parent node',
                       values=color_legend,
                       breaks=tree_levels,
                       labels = as.character(tree_levels))+
    scale_shape_manual('Terminal node',values=c(15,17))+
    theme_void()+
    theme(plot.margin = margin(20, 20, 20, 20),
          legend.key.size = unit(1, "cm"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=10))

  return(g)
}

#' Plots segmentation
#'
#' Plots  focusing on the observed data and shift time and their estimation
#'
#' @param summary data frame, summary data resulting from any segmentation function without estimation of rating curve
#' @param plot_summary list, plot data resulting from any segmentation function
#' @param shift_data_reported vector, shift declared and stored by the hydrometric unit
#' @param show_unc_interval logical, plot of posterior shift estimation as : if `TRUE` uncertainty interval at 95% and if `FALSE` density distribution
#'
#' @return  List with the following components :
#' \enumerate{
#'   \item final_plot: ggplot, observed data presenting the shift time and their estimation
#'   \item observation_and_shift: ggplot, observed data indexed by period, and the estimated shift time is indicated vertically
#'   \item shift_time_density: ggplot, density plot of the shift times following MCMC results, with
#'         vertical lines indicating the 95% credibility interval and a red cross representing shift time assignment
#'   }
#' @export
#' @import patchwork
#' @importFrom ggnewscale new_scale_color
#' @importFrom scales viridis_pal
#' @importFrom stats setNames quantile
plotSegmentation <- function(summary,
                             plot_summary,
                             shift_data_reported=NULL,
                             show_unc_interval=FALSE) {
  # Verify if shift has been declared
  if(!is.null(shift_data_reported)){
    # Is a data frame?
    if(is.data.frame(shift_data_reported))stop('shift_data_reported should not be a data frame')
    # Have they the same format and are they different to numeric format ?
   if(any(class(shift_data_reported)!=class(summary$data$time)) &
      is.numeric(shift_data_reported)!=is.numeric(summary$data$time))stop('shift_data_reported should be in the same format as the time stored in summary')
  }

  if(nrow(summary$shift)==0){

   data=summary$data

   obs_shift_plot=ggplot(data=data)+
     geom_point(aes(x=time,
                    y=obs,
                    col=factor(period)))+
     geom_errorbar(aes(x=time,
                       y=obs,
                       ymin=I95_lower,
                       ymax=I95_upper,
                       col=factor(period)),
                   width=3)

   #Add shift declared
   if(!is.null(shift_data_reported)){

     shift_reported <- data.frame(
       time = shift_data_reported,
       period = 'Shift(s) declared'
     )

     obs_shift_plot=obs_shift_plot +
       geom_point(data=shift_reported,
                  aes(x=time,
                      y=min(data$I95_lower)*0.9,
                      col=factor(period)),
                  shape=4,
                  size=3)
   }

   obs_shift_plot=obs_shift_plot+
     labs(y='Observation',
          x='Time',
          col='Period',
          title = 'Segmentation')+
     scale_color_manual(values = c(getPalette_obs(1),'black'))+
     coord_cartesian(xlim = c(min(data$time),max(data$time)))+
     theme_bw()+
     theme(plot.title = element_text(hjust=0.5,
                                     face='bold',
                                     size=15),
           legend.title.align=0.5)

   pdf_shift_plot = ggplot(data=data,
                           aes(x=time,
                               y=obs))+
     coord_cartesian(xlim = c(min(data$time),max(data$time)))+
     labs(x='Time',
          fill=NULL,
          col=NULL)+
     theme_bw()+
     theme(axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())

   # Check type of plot defined by user
   if(show_unc_interval==TRUE){
     pdf_shift_plot =
       pdf_shift_plot +
       labs(y=NULL)
   }else{
     pdf_shift_plot =
       pdf_shift_plot +
       labs(y='Probability density')
   }
   #Add shift declared
   if(!is.null(shift_data_reported)){

     pdf_shift_plot=pdf_shift_plot +
       geom_point(data=shift_reported,
                  aes(x=time,
                      y=0,
                      col=factor(period)),
                  shape=4,
                  size=3)+
       scale_color_manual(values = c('black'))
   }

   return(list(final_plot=obs_shift_plot/pdf_shift_plot,
               observation_and_shift=obs_shift_plot,
               shift_time_density=pdf_shift_plot))
  }

  # Plot Segmentation showing shift time and observed data indexed by period
  data=summary$data
  shift=summary$shift
  shift$period='Shift(s) detected'

  if(!is.null(shift_data_reported)){
    shift_reported <- data.frame(
      time = shift_data_reported,
      period = 'Shift(s) declared'
    )
  }
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
                   y=min(data$I95_lower),
                   col=factor(period)),
               shape=4,
               size=3)+
    labs(col='Global shift(s) time')

  # Add shift declared
  if(!is.null(shift_data_reported)){
    obs_shift_plot=obs_shift_plot +
      geom_point(data=shift_reported,
                 aes(x=time,
                     y=min(data$I95_lower)*0.9,
                     col=factor(period)),
                 shape=4,
                 size=3)

    if(is.numeric(shift$tau)){
      obs_shift_plot=
        obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red',
                                    'black'),
                           labels=c(round(shift$tau,2),
                                    'Shift(s) detected',
                                    'Shift(s) declared'))

    }else if(lubridate::is.POSIXct(shift$tau)){
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red',
                                    'black'),
                           labels=c(as.character(round(shift$tau,units='days')),
                                    'Shift(s) detected',
                                    'Shift(s) declared'))
    }else{
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red',
                                    'black'),
                           labels=c(as.character(shift$tau),
                                    'Shift(s) detected',
                                    'Shift(s) declared'))
    }

  }else{
    # Any shift declared
    if(is.numeric(shift$tau)){
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red'),
                           labels=c(round(shift$tau,2),
                                    'Shift(s) detected'))

    }else if(lubridate::is.POSIXct(shift$tau)){
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                           'red'),
                           labels=c(as.character(round(shift$tau,units='days')),
                                    'Shift(s) detected'))
    }else{
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                           'red'),
                           labels=c(as.character(shift$tau),
                                    'Shift(s) detected'))
    }
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

  # Customize labels
  label_shift <- paste0('Shift ',seq(1,colourCount_shift))
  if(length(which(colnames(plot_summary$density.inc.tau)=='id_iteration'))!=0){
    label_wrap <- stats::setNames(paste('Shift(s) detected in node',
                                        unique(plot_summary$density.inc.tau$id_iteration),':'),
                                  unique(plot_summary$density.inc.tau$id_iteration))
  }

  # Plot density
  plot_summary$density.inc.tau$period='Shift(s) detected'

  # Check type of plot defined by user
  if(show_unc_interval==TRUE){

    if(length(which(colnames(plot_summary$density.inc.tau)=='id_iteration'))!=0){

      # Assign a "yplot" increasing depending if more than one shift time is detected in a same node
      DF_yplot_assigned=plot_summary$density.tau %>%
        group_by(Shift,id_iteration)%>%
        summarize(.groups = 'drop_last')%>%arrange(id_iteration, Shift) %>%
        mutate(yplot = cur_group_id() )

      IC_data=merge(data.frame(plot_summary$density.tau %>%
                                 group_by(Shift,id_iteration)%>%
                                 summarize(
                                   Lower_inc = stats::quantile(Value, probs=0.025),
                                   .groups = 'drop_last'
                                   )),
                    data.frame(plot_summary$density.tau %>%
                                 group_by(Shift,id_iteration)%>%
                                 summarize(
                                   Upper_inc = stats::quantile(Value,probs= 0.975),
                                   .groups = 'drop_last',
                                 )))

      IC_to_merge=full_join(IC_data,
                            DF_yplot_assigned,
                            by=c('Shift','id_iteration'))

      plot_summary$density.tau =
        full_join(plot_summary$density.tau,
                  IC_to_merge,
                  by=c('Shift','id_iteration'))
      # Assign a "yplot" increasing depending if more than one shift time is detected in a same node
      DF_unc_yplot_assigned=
        plot_summary$density.inc.tau%>%
        group_by(Shift,id_iteration)%>%
        summarize(.groups = 'drop_last')%>%arrange(id_iteration, Shift) %>%
        mutate(yplot = cur_group_id() )

      plot_summary$density.inc.tau =
        full_join(plot_summary$density.inc.tau ,
                  DF_unc_yplot_assigned,
                  by=c('Shift','id_iteration'))
    }else{
      IC_data=merge(data.frame(plot_summary$density.tau %>%
                                 group_by(Shift)%>%
                                 summarize(
                                   Lower_inc = stats::quantile(Value, probs=0.025)
                                 )),
                    data.frame(plot_summary$density.tau %>%
                                 group_by(Shift)%>%
                                 summarize(
                                   Upper_inc = stats::quantile(Value,probs= 0.975)
                                 )
                    ))
      # Add yplot to add for plotting
      IC_data$yplot = seq(nrow(IC_data),1)
      plot_summary$density.inc.tau$yplot= seq(nrow(IC_data),1)

      plot_summary$density.tau =
        full_join(plot_summary$density.tau ,
                  IC_data,
                  by='Shift')
    }


    # Plot limits
    minxplot=min(min(data$time),min(plot_summary$density.tau$Lower_inc))
    maxxplot=max(max(data$time),max(plot_summary$density.tau$Upper_inc))

    pdf_shift_plot =
      ggplot(plot_summary$density.tau,
             aes(x=Value,
                 y = yplot))+
      # Interval
      geom_errorbar(aes(xmin=Lower_inc,
                        xmax=Upper_inc,
                        col=factor(Shift)),
                width=0.08)+
      # MaxPost
      geom_point(data=plot_summary$density.inc.tau,
                 aes(x=taU_MAP,
                     y=yplot,
                     col=factor(period)),
                 shape=4,
                 size=4)
  }else{

    # Plot limits
    minxplot=min(min(data$time),min(plot_summary$density.tau$Value))
    maxxplot=max(max(data$time),max(plot_summary$density.tau$Value))

    pdf_shift_plot = ggplot(plot_summary$density.tau,
                            aes(x=Value,
                                y = Density))+
      geom_area(aes(fill=factor(Shift)),
                alpha=0.4,
                position = "identity")+
      geom_point(data=plot_summary$density.inc.tau,
                 aes(x=taU_MAP,
                     y=density_taU_MAP,
                     col=factor(period)),
                 shape=4,
                 size=3)
  }

  # Add shift declared
  if(!is.null(shift_data_reported)){

    # Check type of plot defined by user
    if(show_unc_interval==TRUE){
      label_unc_interval <- paste0('Uncertainty interval Shift ',seq(1,colourCount_shift))

      pdf_shift_plot=
        pdf_shift_plot +
        geom_point(data=shift_reported,
                   aes(x=time,
                       y=0.8,
                       col=factor(period)),
                   shape=4,
                   size=3)+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_shift),
                           'red',
                           'black'),
                          labels=c(label_unc_interval,
                          'Shift(s) detected',
                          'Shift(s) declared'))+
        labs(y=NULL,
             fill=NULL,
             col=NULL)

    }else{
      pdf_shift_plot =
        pdf_shift_plot +
        geom_point(data=shift_reported,
                   aes(x=time,
                       y=min(plot_summary$density.tau$Density)*0.9,
                       col=factor(period)),
                   shape=4,
                   size=3)+
        scale_color_manual(values=c('red',
                                    'black'),
                           labels=c('Shift(s) detected',
                                    'Shift(s) declared'))+
        scale_fill_manual(values=getPalette_tau_MAP(colourCount_shift),
                          labels=label_shift)+
        labs(y='Probability density',
             fill='Posterior distribution',
             col=NULL)
    }


  }else if(show_unc_interval==TRUE){
    pdf_shift_plot=
      pdf_shift_plot+
      scale_color_manual(values=c(getPalette_tau_MAP(colourCount_shift),
                                  'red'),
                         labels=c(label_unc_interval,
                                  'Shift(s) detected'))+
      labs(y=NULL,
           fill=NULL,
           col=NULL)
  }else{
    pdf_shift_plot=
      pdf_shift_plot+
      scale_color_manual(values='red',
                         labels='Shift(s) detected')+
      scale_fill_manual(values=getPalette_tau_MAP(colourCount_shift),
                        labels=label_shift)+
      labs(y='Probability density',
           fill='Posterior distribution',
           col=NULL)
  }

  pdf_shift_plot=
    pdf_shift_plot+
    labs(x='Time')+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  # Adjust X axis plot
  obs_shift_plot =
    obs_shift_plot +
    coord_cartesian(xlim=c(minxplot,maxxplot),
                    ylim = c(min(data$I95_lower),max(data$I95_upper)))

  pdf_shift_plot =
    pdf_shift_plot+
    coord_cartesian(xlim=c(minxplot,maxxplot))

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
#' @importFrom scales viridis_pal
#' @importFrom stats qnorm
plotRC_ModelAndSegmentation=function(summary,
                                     equation=Exponential_Equation,
                                     ...,
                                     Hmin_user = 0,
                                     Hmax_user = 2,
                                     H_step_discretization=0.01,
                                     autoscale=TRUE,
                                     logscale = FALSE){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with recursive.ModelAndSegmentation function.
                                           If not please use plotSegmentation() function')

  if(identical(equation,BaRatinKAC_Equation)||identical(equation,BaRatinBAC_Equation))stop('To plot the rating curve using Baratin method, you must to use the function plotRCPrediction')

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


  Q_sim_grid <- data.frame(matrix(NA, nrow = length(H_grid), ncol = nrow(param)))

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

  for(i in 1:nrow(param)){ # scan all rating curve estimates
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
  ### limit of plotting
  if(autoscale==TRUE){
    min_x <- min(summary_sim_grid$H_sim,summary_obs$H_obs)
    max_x <- max(summary_sim_grid$H_sim,summary_obs$H_obs)
    min_y <- min(summary_sim_grid$Qsim_I95_lower,summary_obs$Q_I95_lower)
    max_y <- max(summary_sim_grid$Qsim_I95_upper,summary_obs$Q_I95_upper)
  }else{
    H_min_plot_limit = Hmin_user
    H_max_plot_limit = Hmax_user

    subsetzoomdata=summary_obs[dplyr::between(summary_obs$H_obs,H_min_plot_limit,H_max_plot_limit),]

    min_x <- min(summary_sim_grid$H_sim,subsetzoomdata$H_obs)
    max_x <- max(summary_sim_grid$H_sim,subsetzoomdata$H_obs)
    min_y <- min(summary_sim_grid$Qsim_I95_lower,subsetzoomdata$Q_I95_lower)
    max_y <- max(summary_sim_grid$Qsim_I95_upper,subsetzoomdata$Q_I95_upper)
  }

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
    coord_cartesian(xlim=c(H_min_plot_limit,H_max_plot_limit),
                    ylim=c(min_y,max_y))+
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
          legend.title = element_text(hjust=0.5))

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
#' @param ... optional arguments, as vector, to consider shift declared and stored by the hydrometric unit
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
                                        uH=0,
                                        ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with recursive.ModelAndSegmentation function.
                                           If not please use plotSegmentation() function')
  # Adapt summary to use PlotSegmentation function to plot segmentation of stage time series
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$H,
                             I95_lower=summary$data$H-uH,
                             I95_upper=summary$data$H+uH,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plot_segmentation=plotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary,
                                     ...)

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
#' @param ... optional arguments, as vector, to consider shift declared and stored by the hydrometric unit
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
                                        uH=NA,
                                        ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with recursive.ModelAndSegmentation function.
                                           If not please use plotSegmentation() function')
  # Adapt summary to use PlotSegmentation function to plot segmentation of discharge measurements
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$Q,
                             I95_lower=summary$data$Q_I95_lower,
                             I95_upper=summary$data$Q_I95_upper,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plot_segmentation=plotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary,
                                     ...)

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
#' @param ... optional arguments, as vector, to consider shift declared and stored by the hydrometric unit
#'
#' @return ggplot, residual segmentation
#' @export
plotResidual_ModelAndSegmentation <- function(summary,
                                              plot_summary,
                                              ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with recursive.ModelAndSegmentation function.
                                           If not please use plotSegmentation() function')

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
                                plot_summary = plot_summary,
                                ...)

  plotresidual[[1]][[1]]$labels$y='Residual (m3/s) : Observed - simulated '
  plotresidual[[2]]$labels$y='Residual (m3/s) : Observed - simulated'

  return(plotresidual)
}

#' Plot all gaugings after segmentation
#'
#' @param summary list, summary data resulting from the `recursive.ModelAndSegmentation` function
#' @param show_gauging logical, if `TRUE` : gaugings are plotted
#' @param show_RC logical, if `TRUE` : rating curves estimation are plotted
#' @param logscale logical, if `TRUE` : log scale is applied
#'
#' @return ggplot, gaugings after segmentation
#' @export
#' @importFrom scales viridis_pal

plotGaugingsSegmented <- function(summary,
                                  show_gauging=TRUE,
                                  show_RC=FALSE,
                                  logscale=FALSE){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with recursive.ModelAndSegmentation function.
                                           If not please use plotSegmentation() function')
  if(is.null(summary$data))stop('Input data does not match the format of the "recursive.ModelAndSegmentation" results')
  plot_RC_customized = ggplot(summary$data, aes(x=H,
                                                y=Q,
                                                ymin=Q_I95_lower,
                                                ymax=Q_I95_upper,
                                                col=factor(period),
                                                fill=factor(period)))
  if(show_RC){
    plot_RC_customized = plot_RC_customized +
      geom_ribbon(aes(ymin=Qsim_I95_lower,
                      ymax=Qsim_I95_upper,
                      col=NULL),
                  alpha=0.3)
  }

  if(show_gauging){
    plot_RC_customized = plot_RC_customized +
      geom_errorbar()+
      geom_point()
  }

  if(logscale){

    breaks <- 10^seq(floor(log10(min(summary$data$Q))), ceiling(log10(max(summary$data$Q))), by = 1)

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
    plot_RC_customized= plot_RC_customized +
      scale_y_continuous(trans = "log10", breaks = breaks, labels=custom_format)

  }

  # Plot customized
  colourCount_obs = length(unique(summary$data$period))
  getPalette_obs =  scales::viridis_pal(option='D')

  plot_RC_customized=plot_RC_customized+
    scale_color_manual(values = getPalette_obs(colourCount_obs),
                       name = 'Period')+
    scale_fill_manual(values = getPalette_obs(colourCount_obs),
                      name = 'Period')+
    labs(x='Stage (m)',
         y='Discharge (m3/s)')+
    theme_bw()+
    guides(color = guide_legend(override.aes = list(fill = NA)))

  return(plot_RC_customized)
}


#' Plot rating curve from a tree structure
#'
#' Plot all rating curve specified in the nodes after using `recursive.ModelAndSegmentation` function.
#'
#' @param Hgrid data frame, grid user-defined for plotting rating curve
#' @param autoscale logical, auto scale following data for plotting Hgrid
#' @param temp.folder directory, temporary directory to write computations
#' @param CalibrationData character, name of the calibration data used in the `recursive.ModelAndSegmentation` function. It must to match or an error message will be appear
#' @param allnodes logical, if TRUE all node will be plotted and the parameter nodes in the input data function will not be considered. If FALSE, it must to specify the nodes to plot.
#' @param nodes integer vector, the nodes from tree structure for plotting rating curve.
#'
#' @return ggplot, rating curve with MAP, uncertainties and gauging data. More specification in 'Details'
#' @details
#' The function allows plotting any node from the tree structure, as shown in the example from `recursive.ModelAndSegmentation`.
#' The rating curves consist of the optimal rating curve after estimation, along with the parametric and total uncertainty.
#' Gauging data used for calibration during segmentation have been plotted. Hence, the smaller the number of the node, the more gauging data have been used for calibration.
#' @export
#' @importFrom RBaM prediction
plotRCPrediction <- function(Hgrid=data.frame(grid=seq(-1,2,by=0.01)),
                             autoscale=FALSE,
                             temp.folder=file.path(tempdir(),'BaM'),
                             CalibrationData='CalibrationData.txt',
                             allnodes=FALSE,
                             nodes=1){

  if(!is.data.frame(Hgrid))(stop('Hgrid must be a data frame'))
  if(!dir.exists(file.path(temp.folder,'it_1')))(stop('Segmentation using the BaRatin method is required before using this function.
                                                      Please, ensure that you have specified the correct path in input data'))

  if(allnodes==FALSE & any(nodes>length(list.files(temp.folder,pattern = 'it_'))))stop('At least one specified node does not exist. Check node input data')
  if(allnodes==FALSE & !is.vector(nodes))(stop('Nodes must be a vector'))
  if(any(CalibrationData==list.files(file.path(temp.folder,'it_1')))==FALSE)stop('CalibrationData given in input data does not exist in the directory specified in temp.folder. Please, check the name of calibration data used in the recursive.ModelAndSegmentation function')

  # Extract all data required from specified nodes
  allData=list()

  if(allnodes){
    nodes_to_plot=seq(1,length(grep(list.dirs(temp.folder,recursive  = FALSE),pattern='it_')==T))
  }else{
    nodes_to_plot=nodes
  }

  for(i in 1:length(nodes_to_plot)){

    temp.folder.RCPlot=file.path(temp.folder,paste0('it_',nodes_to_plot[[i]]))

    CalData=read.table(file.path(temp.folder.RCPlot,
                                 CalibrationData),
                       header = T)

    # Upload Data object
    load(file.path(temp.folder.RCPlot,'DataObject.RData'))

    # Upload Model object
    load(file.path(temp.folder.RCPlot,'ModelObject.RData'))

    if(autoscale==TRUE){
      Hgrid=data.frame(grid=seq(min(CalData$H),max(CalData$H),by=0.01))
    }
    totalU=RBaM::prediction(X=Hgrid, # stage values
                            spagFiles='QRC_TotalU.spag', # file where predictions are saved
                            data.dir=temp.folder.RCPlot, # a copy of data files will be saved here
                            doParametric=TRUE, # propagate parametric uncertainty, i.e. MCMC samples?
                            doStructural=TRUE) # propagate structural uncertainty ?

    # Define a 'prediction' object for parametric uncertainty only - not the doStructural=FALSE
    paramU=RBaM::prediction(X=Hgrid,
                            spagFiles='QRC_ParamU.spag',
                            data.dir=temp.folder.RCPlot,
                            doParametric=TRUE,
                            doStructural=FALSE)

    # Define a 'prediction' object with no uncertainty - this corresponds to the 'maxpost' RC maximizing the posterior density
    maxpost=RBaM::prediction(X=Hgrid,
                             spagFiles='QRC_Maxpost.spag',
                             data.dir=temp.folder.RCPlot,
                             doParametric=FALSE,
                             doStructural=FALSE)

    # # Define a 'prediction' object with no uncertainty - this corresponds to the rating curve computed with prior information
    # priorU=RBaM::prediction(X=Hgrid,
    #                         spagFiles='QRC_Prior.spag',
    #                         data.dir=temp.folder.RCPlot,
    #                         priorNsim = nrow(mcmc.segm),
    #                         doParametric=TRUE,
    #                         doStructural=FALSE)

    RBaM:: BaM(mod=M,
               data=D,
               workspace = temp.folder.RCPlot,
               dir.exe = file.path(find.package("RBaM"), "bin"),
               pred=list(totalU,paramU,maxpost), # list of predictions
               # pred=priorU, # list of predictions
               doCalib=FALSE,
               doPred=TRUE)

    # Get total uncertainty from the rating curve
    TotalUenv=read.table(file.path(temp.folder.RCPlot,'QRC_TotalU.env'),header = TRUE)
    TotalUenv=replace_negatives_or_zero_values(data_frame = TotalUenv,
                                               columns = 'all',
                                               consider_zero = TRUE,
                                               replace = 0 )

    # Add parametric uncertainty
    ParametricUenv=read.table(file.path(temp.folder.RCPlot,'QRC_ParamU.env'),header = TRUE)
    ParametricUenv=replace_negatives_or_zero_values(data_frame = ParametricUenv,
                                                    columns = 'all',
                                                    consider_zero = TRUE,
                                                    replace = 0 )


    # Add maxpost rating curve
    MAPRC=read.table(file.path(temp.folder.RCPlot,'QRC_Maxpost.spag'))

    allData[[i]]=list(CalData=CalData,
                      TotalUenv=TotalUenv,
                      ParametricUenv=ParametricUenv,
                      MAPRC=MAPRC,
                      HgridPlot=Hgrid)
  }

  # Plots
  PlotRCPred=list()
  for (i in 1:length(allData)){
    PlotRCPred[[i]]<-local({
      inner_list <- allData[[i]]
      i<-i

      ggplot()+
        # Total uncertainty
        geom_ribbon(data = inner_list$TotalUenv,
                    aes(x=inner_list$HgridPlot[,1],
                        ymin=q2.5,
                        ymax=q97.5,
                        fill='Total'),
                    alpha=0.65)+
        # Parametric uncertainty
        geom_ribbon(data = inner_list$ParametricUenv,
                    aes(x=inner_list$HgridPlot[,1],
                        ymin=q2.5,
                        ymax=q97.5,
                        fill='Parametric'),
                    alpha=0.65)+
        # Map
        geom_line(data=inner_list$MAPRC,
                  aes(x=inner_list$HgridPlot[,1],
                      y=inner_list$MAPRC[,1],
                      col='MAP'))+
        # Guagings
        geom_errorbar(data = inner_list$CalData,
                      aes(x=H,
                          ymin=Q-uQ,
                          ymax=Q+uQ,
                          col='Gaugings'),
                      width=0.05)+
        geom_point(data = inner_list$CalData,
                   aes(x=H,
                       y=Q,
                       col='Gaugings'))+
        labs(title=paste0('Rating curve estimation for node ',nodes_to_plot[[i]]),
             x='H[m]',
             y='Q[m3/s]',
             fill='Uncertainty',
             col=NULL)+
        scale_fill_manual(values=c('pink', 'red'))+
        scale_color_manual(values=c('blue','black'))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))
    })
  }
  return(PlotRCPred)
}


#' Plot recession extracted
#'
#' @param Rec_extracted data frame, recession extracted after using `Extraction_recession`
#' @param error_bar_plot logical, `TRUE` = plot error bar.
#'
#' @return list, plots of recession extracted. First plot indicates recession following stage record, second plot shows recession time and stage and third plot indicates the minimum recession stage.
#' @export
#' @importFrom RColorBrewer brewer.pal
plot_rec_extracted <- function(Rec_extracted,
                               error_bar_plot = FALSE){

  colors=RColorBrewer::brewer.pal(10,'Paired')

  # Plot of the extracted stage-recession
  rec.plot=ggplot(Rec_extracted,
                  aes(x=date,
                      y=hrec,
                      ymin=hrec-uHrec,
                      ymax=hrec+uHrec,
                      col=indx))+
    geom_point()

  if(error_bar_plot){
    rec.plot=rec.plot+
      geom_errorbar()
  }

  rec.plot=rec.plot+
    theme_bw()+
    labs(x = 'Time [date]',
         y = 'Stage')+
    guides(color='none')+
    scale_color_gradientn(colors=colors)

  # Plot of the recession
  rec.plot2=ggplot(Rec_extracted,
                   aes(x=time_rec,
                       y=hrec,
                       ymin=hrec-uHrec,
                       ymax=hrec+uHrec,
                       col=indx))+
    geom_point()

  if(error_bar_plot){
    rec.plot2=rec.plot2+
      geom_errorbar()
  }

  position_x=quantile(Rec_extracted$time_rec,probs = 0.95)
  position_y=quantile(Rec_extracted$hrec,probs = 0.99)

  test_long_serie=as.data.frame(Rec_extracted%>%
                                  group_by(indx)%>%
                                  summarize(max(time_rec)))

  rec.plot2=
    rec.plot2 +
    annotate("text",
             x = position_x,
             y = position_y,
             label = paste0("Total points = ",nrow(Rec_extracted)))+
    annotate("text",
             x = position_x,
             y = position_y*0.9,
             label = paste0("Total recessions = ",max(unique(Rec_extracted$indx))))+
    annotate("text",
             x = position_x,
             y = position_y*0.8,
             label = paste0("Total recessions (t > 80 days) = ", length(which(test_long_serie>80))))+
    theme_bw()+
    labs(x = 'Recession time [day]',
         y = 'Stage',
         col = '# Recession')+
    scale_color_gradientn(colors=colors)


  DF_h_min_rec <- data.frame(Rec_extracted %>%
                                 group_by(indx) %>%
                                 arrange(hrec, desc(time_rec)) %>% # Trier par hrec croissant puis par time_rec décroissant
                                 slice_head(n = 1) %>%
                                 ungroup())

  # Plot minimum h recession values
  plot_min_h_rec=
    ggplot(DF_h_min_rec,aes(x=date,
                            y=hrec,
                            ymin=hrec-uHrec,
                            ymax=hrec+uHrec))+
    geom_point()

  if(error_bar_plot){
    plot_min_h_rec=plot_min_h_rec+
      geom_errorbar()
  }

  plot_min_h_rec=plot_min_h_rec +
    theme_bw()+
    labs(x = 'Time [date]',
         y = 'Minimum Recession Stage H[m]')

  return(list(rec.plot,
              rec.plot2,
              plot_min_h_rec))
}

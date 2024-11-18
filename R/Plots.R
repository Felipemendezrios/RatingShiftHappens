#' Plot recursive segmentation tree
#'
#' Plot the tree resulting from any recursive segmentation procedures
#'
#' @param tree data frame, tree resulting from the call of Recursive_Segmentation and ModelAndSegmentation_Gaugings
#'
#' @return a ggplot
#'
#' @examples
#'
#' # Apply recursive segmentation
#' results=Recursive_Segmentation(obs=CongoRiverMINAN$Q,
#'                                time=CongoRiverMINAN$year,
#'                                u=CongoRiverMINAN$uQ,
#'                                nSmax=2)
#'
#' # plot recursion tree
#' PlotTree(results$tree)
#' @export
#' @import ggplot2
#' @importFrom stats rnorm setNames
#' @importFrom scales viridis_pal
PlotTree <- function(tree){
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
#'    \item final_plot: ggplot, observed data presenting the shift time and their estimation
#'    \item observation_and_shift: ggplot, observed data indexed by period, and the estimated shift time is indicated vertically,
#'    \item shift_time_density: ggplot, density plot of the shift times following MCMC results, with
#'          vertical lines indicating the 95% credibility interval and a red cross representing shift time assignment
#' }
#' @export
#' @import patchwork
#' @importFrom ggnewscale new_scale_color
#' @importFrom scales viridis_pal
#' @importFrom stats setNames quantile
PlotSegmentation <- function(summary,
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

  if(is.null(summary$shift)){

   data=summary$data
   getPalette_obs =  scales::viridis_pal(option='D')

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
               linewidth = 0.5,
               alpha=0.8)+
    geom_point(data=shift,
               aes(x=tau,
                   y=min(data$I95_lower),
                   col=factor(period)),
               shape=4,
               size=3)+
    labs(col='Global shift(s) \ntime')

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

    obs_shift_plot=obs_shift_plot+# Use guides to separate point and line legend representation
      guides(color = guide_legend(override.aes = list(
        linetype = c(rep("solid",colourCount_tau),'blank','blank'),   # Line for 'Activation Height', no line for 'Gaugings'
        shape = c(rep(NA,colourCount_tau), 4,4)      # Point shape for 'Gaugings', no shape for 'Activation Height'
      )))

  }else{
    # Any shift declared
    if(is.numeric(shift$tau)){
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red'),
                           labels=c(round(shift$tau,2),
                                    'Shift(s) \ndetected'))

    }else if(lubridate::is.POSIXct(shift$tau)){
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red'),
                           labels=c(as.character(round(shift$tau,units='days')),
                                    'Shift(s) \ndetected'))
    }else{
      obs_shift_plot=obs_shift_plot+
        scale_color_manual(values=c(getPalette_tau_MAP(colourCount_tau),
                                    'red'),
                           labels=c(as.character(shift$tau),
                                    'Shift(s) \ndetected'))
    }
    obs_shift_plot=obs_shift_plot+# Use guides to separate point and line legend representation
      guides(color = guide_legend(override.aes = list(
        linetype = c(rep("solid",colourCount_tau),'blank'),   # Line for 'Activation Height', no line for 'Gaugings'
        shape = c(rep(NA,colourCount_tau), 4)      # Point shape for 'Gaugings', no shape for 'Activation Height'
      )))
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
                                    'Shift(s) \ndetected',
                                    'Shift(s) \ndeclared'))+
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
                           labels=c('Shift(s) \ndetected',
                                    'Shift(s) \ndeclared'))+
        scale_fill_manual(values=getPalette_tau_MAP(colourCount_shift),
                          labels=label_shift)+
        labs(y='Probability density',
             fill='Posterior \ndistribution',
             col=NULL)
    }


  }else if(show_unc_interval==TRUE){
    pdf_shift_plot=
      pdf_shift_plot+
      scale_color_manual(values=c(getPalette_tau_MAP(colourCount_shift),
                                  'red'),
                         labels=c(label_unc_interval,
                                  'Shift(s) \ndetected'))+
      labs(y=NULL,
           fill=NULL,
           col=NULL)
  }else{
    pdf_shift_plot=
      pdf_shift_plot+
      scale_color_manual(values='red',
                         labels='Shift(s) \ndetected')+
      scale_fill_manual(values=getPalette_tau_MAP(colourCount_shift),
                        labels=label_shift)+
      labs(y='Probability density',
           fill='Posterior \ndistribution',
           col=NULL)
  }

  pdf_shift_plot=
    pdf_shift_plot+
    labs(x='Time',
         title='Shift(s) detected in node:')+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))

  # Adjust X axis plot
  obs_shift_plot =
    obs_shift_plot +
    coord_cartesian(xlim=c(minxplot,maxxplot),
                    ylim = c(min(data$I95_lower),max(data$I95_upper)))

  pdf_shift_plot =
    pdf_shift_plot+
    coord_cartesian(xlim=c(minxplot,maxxplot))

  if(length(which(colnames(plot_summary$density.inc.tau)=='id_iteration'))!=0){
    pdf_shift_plot=
      pdf_shift_plot+
      facet_grid(rows=vars(id_iteration),
                 scales='free_y',
                 space = 'free_y')
  }

  # Customize legend position
  obs_shift_plot = obs_shift_plot +
    theme(
      legend.position = "right",
      legend.box = "horizontal"
    ) +
    guides(
      color = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1)
    )

  pdf_shift_plot = pdf_shift_plot +
    theme(
      legend.position = "right",
      legend.box = "horizontal"
    ) +
    guides(
      color = guide_legend(ncol = 1),
      shape = guide_legend(ncol = 1)
    )

  final_plot = obs_shift_plot/pdf_shift_plot +  guides(col = "none")

  final_plot = final_plot + plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      legend.box = "horizontal"
    )

  return(list(final_plot=final_plot,
              observation_and_shift=obs_shift_plot,
              shift_time_density=pdf_shift_plot))
}

#' Plot Rating curve after segmentation
#'
#' Plot of all the rating curve after using the function `ModelAndSegmentation_Gaugings`, along with associated uncertainties
#'
#' @param summary data.frame, summary data resulting from the function `ModelAndSegmentation_Gaugings` saved in `$summary`
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
#' For the extra information, you can view the arguments through `args()`  by enclosing the equation specified in the input arguments in brackets as shown in the example of `?ModelAndSegmentation_Gaugings``.
#' Information about parameters values are available from the function `ModelAndSegmentation_Gaugings` saved in `$summary$param.equation`. Save each parameter estimates separately
#' Please ensure that you enter the same number of rating curve parameters in the extra information as shown in the example `?ModelAndSegmentation_Gaugings`.
#' @export
#' @importFrom scales viridis_pal
#' @importFrom stats qnorm
PlotRCSegmentation_Gaugings=function(summary,
                                     equation=EquationRC_Exponential,
                                     ...,
                                     Hmin_user = 0,
                                     Hmax_user = 2,
                                     H_step_discretization=0.01,
                                     autoscale=TRUE,
                                     logscale = FALSE){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with ModelAndSegmentation_Gaugings function.
                                           If not please use PlotSegmentation() function')

  if(identical(equation,EquationRC_BaRatinKAC)||identical(equation,EquationRC_BaRatinBAC))stop('To plot the rating curve using Baratin method, you must to use the function PlotRCSegmentation_Gaugings_Tree')

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

    # Custom label formatting function (custom_format)
    RC_plot= RC_plot+ scale_y_continuous(trans = "log10", breaks = breaks, labels=custom_format)
  }

  return(RC_plot)
}

#' Plot stage time series after segmentation
#'
#' Plots focusing on the observed data and shift times after using `ModelAndSegmentation_Gaugings` function, along with associated uncertainties
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
PlotSegmentation_Gaugings_Hdt <- function(summary,
                                          plot_summary,
                                          uH=0,
                                          ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with ModelAndSegmentation_Gaugings function.
                                           If not please use PlotSegmentation() function')
  # Adapt summary to use PlotSegmentation function to plot segmentation of stage time series
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$H,
                             I95_lower=summary$data$H-uH,
                             I95_upper=summary$data$H+uH,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plot_segmentation=PlotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary,
                                     ...)

  plot_segmentation[[1]][[1]]$labels$y='Stage record (m)'
  plot_segmentation[[2]]$labels$y='Stage record (m)'

  return(plot_segmentation)
}
#' Plot discharge time series after segmentation
#'
#' Plots focusing on the observed data and shift times after using `ModelAndSegmentation_Gaugings` function, along with associated uncertainties
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
PlotSegmentation_Gaugings_Qdt <- function(summary,
                                          plot_summary,
                                          uH=NA,
                                          ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with ModelAndSegmentation_Gaugings function.
                                           If not please use PlotSegmentation() function')
  # Adapt summary to use PlotSegmentation function to plot segmentation of discharge measurements
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$Q,
                             I95_lower=summary$data$Q_I95_lower,
                             I95_upper=summary$data$Q_I95_upper,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plot_segmentation=PlotSegmentation(summary = summary_plot,
                                     plot_summary = plot_summary,
                                     ...)

  plot_segmentation[[1]][[1]]$labels$y='Discharge (m3/s)'
  plot_segmentation[[2]]$labels$y='Discharge (m3/s)'

  return(plot_segmentation)
}

#' Plot residuals with updating rating curve
#'
#' Plot residuals after using `ModelAndSegmentation_Gaugings` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#' @param plot_summary list, plot data resulting from any segmentation function
#' @param ... optional arguments, as vector, to consider shift declared and stored by the hydrometric unit
#'
#' @return ggplot, residual segmentation
#' @export
PlotSegmentation_Gaugings_Residual <- function(summary,
                                               plot_summary,
                                               ...){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with ModelAndSegmentation_Gaugings function.
                                                         If not please use PlotSegmentation() function')

  if(is.null(summary$shift))stop('Any shift time detected')

  # Estimate residual uncertainty :
  u_res=sqrt(summary$data$uQ_sim^2+summary$data$uQ^2)

  # Adapt summary to use PlotSegmentation function to plot segmentation of discharge measurements
  data_adapted <- data.frame(time=summary$data$time,
                             obs=summary$data$Qres,
                             I95_lower=summary$data$Qres+stats::qnorm(0.025)*u_res,
                             I95_upper=summary$data$Qres+stats::qnorm(0.975)*u_res,
                             period=summary$data$period)

  summary_plot = list(data=data_adapted,
                      shift=summary$shift)

  plotresidual=PlotSegmentation(summary = summary_plot,
                                plot_summary = plot_summary,
                                ...)

  plotresidual[[1]][[1]]$labels$y=expression(Residuals ~ m^3/s)
  plotresidual[[2]]$labels$y=expression(Residuals ~ m^3/s)

  return(plotresidual)
}

#' Plot all gaugings after segmentation
#'
#' @param summary list, summary data resulting from the `ModelAndSegmentation_Gaugings` function
#' @param show_gauging logical, if `TRUE` : gaugings are plotted
#' @param show_RC logical, if `TRUE` : rating curves estimation are plotted
#' @param logscale logical, if `TRUE` : log scale is applied
#'
#' @return ggplot, gaugings after segmentation
#' @export
#' @importFrom scales viridis_pal

PlotSegmentation_Gaugings_QdH <- function(summary,
                                          show_gauging=TRUE,
                                          show_RC=FALSE,
                                          logscale=FALSE){
  if(length(which(colnames(summary$data)=='H'))==0)stop('Be sure that segmentation has been computed with ModelAndSegmentation_Gaugings function.
                                           If not please use PlotSegmentation() function')
  if(is.null(summary$data))stop('Input data does not match the format of the "ModelAndSegmentation_Gaugings" results')
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

    # Custom label formatting function (custom_format)
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
#' Plot all rating curve specified in the nodes after using `ModelAndSegmentation_Gaugings` function.
#'
#' @param Hgrid data frame, grid user-defined for plotting rating curve
#' @param autoscale logical, auto scale following data for plotting Hgrid
#' @param temp.folder directory, temporary directory to write computations
#' @param CalibrationData character, name of the calibration data used in the `ModelAndSegmentation_Gaugings` function. It must to match or an error message will be appear
#' @param allnodes logical, if TRUE all node will be plotted and the parameter nodes in the input data function will not be considered. If FALSE, it must to specify the nodes to plot.
#' @param nodes integer vector, the nodes from tree structure for plotting rating curve.
#'
#' @return List with the following components :
#' \enumerate{
#'     \item ggplot, rating curve with MAP, uncertainties and gauging data. More specification in 'Details'
#'     \item list, data frame with numeric data of each rating curve plotted before
#' }
#' @details
#' The function allows plotting any node from the tree structure, as shown in the example from `ModelAndSegmentation_Gaugings`.
#' The rating curves consist of the optimal rating curve after estimation, along with the parametric and total uncertainty.
#' Gauging data used for calibration during segmentation have been plotted. Hence, the smaller the number of the node, the more gauging data have been used for calibration.
#' @export
#' @importFrom RBaM prediction
#' @importFrom stringr str_detect
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
PlotRCSegmentation_Gaugings_Tree <- function(Hgrid=data.frame(grid=seq(-1,2,by=0.01)),
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
  if(any(CalibrationData==list.files(file.path(temp.folder,'it_1')))==FALSE)stop('CalibrationData given in input data does not exist in the directory specified in temp.folder. Please, check the name of calibration data used in the ModelAndSegmentation_Gaugings function')

  # Extract all data required from specified nodes
  allData=list()

  if(allnodes){
    nodes_to_plot=seq(1,length(grep(list.dirs(temp.folder,recursive  = FALSE),pattern='it_')==T))
  }else{
    nodes_to_plot=nodes
  }

  # Error model fixed for all functions in the package
  remnant_prior <- list(RBaM::remnantErrorModel(funk = "Linear",
                                                par = list(RBaM::parameter(name="gamma1",
                                                                           init=1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)),
                                                           RBaM::parameter(name="gamma2",
                                                                           init=0.1,
                                                                           prior.dist = "Uniform",
                                                                           prior.par = c(0,1000)))))

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
               remnant = remnant_prior,
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

    #Read stage activation to plot

    summary_estimation=read.table(file.path(temp.folder.RCPlot,
                                            "Results_Summary.txt"),
                                  header = T)
    k_MAP=unlist(summary_estimation[which(c('MaxPost')==rownames(summary_estimation)),
                                    which(stringr::str_detect(colnames(summary_estimation),'k'))])

    MCMC_cooked=read.table(file.path(temp.folder.RCPlot,
                                     "Results_Cooking.txt"),
                           header = T)

    k_MCMC=MCMC_cooked[,which(stringr::str_detect(colnames(MCMC_cooked),'k'))]

    k_uncertainty=data.frame(apply(k_MCMC, 2, function(x) quantile(x, probs = c(0.025, 0.975))))

    # Adapt data frame to plot

    # Convert the data frame to long format
    k_long <- k_uncertainty %>%
      tibble::rownames_to_column("percentile") %>%
      tidyr::pivot_longer(-percentile, names_to = "k", values_to = "value")

    # Create a data frame for credibility intervals
    Credibility_interval_data <- data.frame(
      k = unique(k_long$k),  # Repeat each k for 2 rows (2.5% and 97.5%)
      xmin = k_long$value[k_long$percentile == "2.5%"],
      xmax = k_long$value[k_long$percentile == "97.5%"],
      ymin = 0,
      ymax = max(TotalUenv$q97.5)
    )

    allData[[i]]=list(CalData=CalData,
                      TotalUenv=TotalUenv,
                      ParametricUenv=ParametricUenv,
                      MAPRC=MAPRC,
                      HgridPlot=Hgrid,
                      k_MAP=k_MAP,
                      k_uncertainty=Credibility_interval_data)
  }

  # Plots
  PlotRCPred=list()
  for (i in 1:length(allData)){
    PlotRCPred[[i]]<-local({
      inner_list <- allData[[i]]
      i<-i

      ggplot()+
        # Stage activation estimated (first control)
        geom_vline(aes(xintercept=inner_list$k_MAP),
                   col='green',
                   linewidth=0.5,
                   show.legend = FALSE)+
        geom_rect(data = inner_list$k_uncertainty,
                  aes(xmin = xmin,
                      xmax = xmax,
                      ymin = ymin,
                      ymax = ymax,
                      fill = "Posterior \nactivation stage"),
                  alpha = 0.15)+
        # Total uncertainty
        geom_ribbon(data = inner_list$TotalUenv,
                    aes(x=inner_list$HgridPlot[,1],
                        ymin=q2.5,
                        ymax=q97.5,
                        fill='Posterior \nparametric and \nstructural \nuncertainty'),
                    alpha=0.65)+
        # Parametric uncertainty
        geom_ribbon(data = inner_list$ParametricUenv,
                    aes(x=inner_list$HgridPlot[,1],
                        ymin=q2.5,
                        ymax=q97.5,
                        fill='Posterior \nparametric \nuncertainty'),
                    alpha=0.65)+
        # Map
        geom_line(data=inner_list$MAPRC,
                  aes(x=inner_list$HgridPlot[,1],
                      y=inner_list$MAPRC[,1],
                      col='Posterior \nrating curve'))+
        # Gaugings
        geom_errorbar(data = inner_list$CalData,
                      aes(x=H,
                          ymin=Q+stats::qnorm(0.025)*uQ,
                          ymax=Q+stats::qnorm(0.975)*uQ,
                          col='Gaugings'),
                      width=0.05)+
        geom_point(data = inner_list$CalData,
                   aes(x=H,
                       y=Q,
                       col='Gaugings'))+
        coord_cartesian(xlim = c(min(inner_list$k_uncertainty$xmin),
                                 max(inner_list$HgridPlot)))+
        labs(title=paste0('Rating curve estimation for node ',nodes_to_plot[[i]]),
             x='H[m]',
             y='Q[m3/s]',
             fill='Uncertainty',
             col=NULL)+
        scale_fill_manual(values=c('green','pink', 'red'))+
        scale_color_manual(values=c('blue','black'))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))+
        guides(color = guide_legend(override.aes = list(
          shape = c(19,NA)
        )))
    })
  }
  return(list(PlotRCPred=PlotRCPred,
              allData=allData))
}


#' Plot recession extracted
#'
#' @param Rec_extracted data frame, recession extracted after using `Extract_Recessions`
#' @param error_bar_plot logical, `TRUE` = plot error bar.
#'
#' @return list, plots of recession extracted. First plot indicates recession following stage record, second plot shows recession time and stage and third plot indicates the minimum recession stage.
#' @export
#' @importFrom RColorBrewer brewer.pal
PlotExtract_Recessions <- function(Rec_extracted,
                                   error_bar_plot = FALSE){

  colors=RColorBrewer::brewer.pal(10,'Paired')

  # Plot of the extracted stage-recession
  rec.plot=ggplot(Rec_extracted,
                  aes(x=date,
                      y=hrec,
                      ymin=hrec+stats::qnorm(0.025)*uHrec,
                      ymax=hrec+stats::qnorm(0.975)*uHrec,
                      col=indx))+
    geom_point()

  if(error_bar_plot){
    rec.plot=rec.plot+
      geom_errorbar()
  }

  rec.plot=rec.plot+
    theme_bw()+
    labs(x = 'Time [date]',
         y = 'Stage [m]')+
    guides(color='none')+
    scale_color_gradientn(colors=colors)

  # Plot of the recession
  rec.plot2=ggplot(Rec_extracted,
                   aes(x=time_rec,
                       y=hrec,
                       ymin=hrec+stats::qnorm(0.025)*uHrec,
                       ymax=hrec+stats::qnorm(0.975)*uHrec,
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
                               arrange(hrec, desc(time_rec)) %>%
                               slice_head(n = 1) %>%
                               ungroup())

  # Plot minimum h recession values
  plot_min_h_rec=
    ggplot(DF_h_min_rec,aes(x=date,
                            y=hrec))+
    geom_point()

  if(error_bar_plot){
    plot_min_h_rec=plot_min_h_rec+
      geom_errorbar()
  }

  plot_min_h_rec=plot_min_h_rec +
    theme_bw()+
    labs(x = 'Time [date]',
         y = 'Minimum recession stage H[m]')

  return(list(rec.plot,
              rec.plot2,
              plot_min_h_rec))
}


#' Comparison between simulated and observed data-recession
#'
#' @param model_rec list, results obtained by using `ModelAndSegmentation_Recessions` function
#' @param spec_recession integer vector, number of recession to plot the observed and simulated recession data separately
#' @param recession_rejected logical, `TRUE` = includes plot of all recession rejected
#' @param all_recession logical, `TRUE` = plot all recessions with observed and simulated recession data
#'
#' @return ggplot, comparison between recession observed and simulated and the status if this is accepted or rejected
#'
#' @details
#' If all_recession = `TRUE`, spec_recession will not be plotted, because they have already been plotted
#' @export
#' @importFrom gridExtra arrangeGrob grid.arrange
PlotSimObs_Recessions <- function(model_rec,
                                  spec_recession=NULL,
                                  recession_rejected=TRUE,
                                  all_recession=FALSE){

  if(!is.null(spec_recession)){
    if(any(spec_recession<=0))stop('spec_recession must be positive')
    if(any(spec_recession %% 1 != 0))stop('spec_recession must be a integer')
    if(any(spec_recession > model_rec$n.rec.max))stop('spec_recession must be between the number of curves specified in the data')
  }
  if(all(is.null(spec_recession) & isFALSE(all_recession) &
         isFALSE(recession_rejected)))stop('A recession must be selected in spec_recession or all recesssion or the plot recession rejected')

  if(all_recession == TRUE){
    DF.obs.vs.sim = model_rec$summary.residual
  }else if(!is.null(spec_recession)){
    if(recession_rejected == TRUE){
      # Get recession rejected
      indx.rec.rejected=unique(model_rec$summary.residual$indx[which(model_rec$summary.residual$status=='Rejected')])
      indx.plot=c(spec_recession,indx.rec.rejected)
      DF.obs.vs.sim = model_rec$summary.residual[model_rec$summary.residual$indx %in% indx.plot,]
    }else{
      DF.obs.vs.sim = model_rec$summary.residual[model_rec$summary.residual$indx %in% spec_recession,]
    }
  }else if(recession_rejected == TRUE){
    # Get recession rejected
    indx.rec.rejected=unique(model_rec$summary.residual$indx[which(model_rec$summary.residual$status=='Rejected')])
    if(length(indx.rec.rejected)==0)stop('No simulated recessions have been rejected. Please specify either all recessions or only some')
    DF.obs.vs.sim = model_rec$summary.residual[model_rec$summary.residual$indx %in% indx.rec.rejected,]
  }

  # Plot MAP simulation and observed data
  PlotREC.obs.sim=list()
  indx_to_plot=unique(DF.obs.vs.sim$indx)

  # Plot the first plot to get legend
  data_sample <- DF.obs.vs.sim[which(DF.obs.vs.sim$indx==indx_to_plot[1]),]

  plot_sample =
    ggplot(data=data_sample,
           aes(x=X1_obs,
               y=Y1_obs,
               ymin=Y1_obs+stats::qnorm(0.025)*uH_obs,
               ymax=Y1_obs+stats::qnorm(0.975)*uH_obs,
               col=factor('Observations')))+
    geom_pointrange()+
    geom_line(aes(x=X1_obs,
                  y=Y1_sim,
                  col=factor('Simulation (MAP)')))+
    labs(title=paste0('Recession number ',indx_to_plot[1]),
         subtitle = paste0('Status: ', unique(data_sample$status)),
         x='Time [days]',
         y='H [m]',
         col= NULL)+
    scale_color_manual(values=c('Simulation (MAP)'='orange',
                                'Observations'='black'))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(color = ifelse(unique(data_sample$status)=='Accepted',
                                                     'green4',
                                                     'red'),
                                       size=12,
                                       hjust = 0.5),
          legend.key.height = unit(0.8,'cm'),
          legend.text = element_text(size=12))

  common_legend=get_legend(plot_sample)

  for (i in 1:length(indx_to_plot)){

    PlotREC.obs.sim[[i]]<-local({
      inner_list <- DF.obs.vs.sim[which(DF.obs.vs.sim$indx==indx_to_plot[i]),]

      ggplot(data=inner_list,
             aes(x=X1_obs,
                 y=Y1_obs,
                 ymin=Y1_obs+stats::qnorm(0.025)*uH_obs,
                 ymax=Y1_obs+stats::qnorm(0.975)*uH_obs,
                 col=factor('Observations')))+
        geom_pointrange()+
        geom_line(aes(x=X1_obs,
                      y=Y1_sim,
                      col=factor('Simulation (MAP)')))+
        labs(title=paste0('Recession number ',indx_to_plot[i]),
             subtitle = paste0('Status: ', unique(inner_list$status)),
             x='Time [days]',
             y='H [m]',
             col= NULL)+
        scale_color_manual(values=c('Simulation (MAP)'='orange',
                                    'Observations'='black'))+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(color= ifelse(unique(inner_list$status)=='Accepted',
                                                        'green4',
                                                        'red'),
                                          size=12,
                                          hjust = 0.5),
              legend.position = "none")
    })
  }
  # Define number of columns
  ncol_plot = ifelse(length(PlotREC.obs.sim)==1,1,ifelse(length(PlotREC.obs.sim)<=6,2,4))

  grobs <- lapply(PlotREC.obs.sim, ggplotGrob)

  grid_PlotREC.obs.sim <- gridExtra::grid.arrange(
    do.call(gridExtra::arrangeGrob, c(grobs, ncol = ncol_plot)),
    common_legend,
    widths = c(5,1))

  return(grid_PlotREC.obs.sim)
}


#' Plot recession modeled and segmented
#'
#' @param model_rec list, results obtained by using `ModelAndSegmentation_Recessions` function
#'
#' @return list of plots with the following components:
#'  \enumerate{
#'   \item plot.rec.segmented: ggplot, recession data indexed by period, and  estimated shift time is indicated vertically
#'   \item plot.b.segmented: ggplot, asymptotic height indexed by period, and the estimated shift time is indicated vertically
#'  }
#' @export
PlotSegmentation_Recessions <- function(model_rec){

  # Plot: recession extracted showing the segmentation
  h.t.summary = model_rec$summary.rec.extracted

  # Adapt the data frame to PlotSegmentation function
  h.t.summary$data$I95_upper=h.t.summary$data$hrec+(stats::qnorm(0.975)*h.t.summary$data$uHrec) # 95% interval uncertainty
  h.t.summary$data$I95_lower=h.t.summary$data$hrec+(stats::qnorm(0.025)*h.t.summary$data$uHrec) # 95% interval uncertainty

  colnames(h.t.summary$data) <- c('indx','time_rec','time','obs','u','status','period','I95_lower','I95_upper')

  plot.rec.segmented = PlotSegmentation(summary=h.t.summary,
                                        plot_summary=model_rec$plots)$final_plot

  plot.rec.segmented[[1]] = plot.rec.segmented[[1]]+
    ylab('Recession-stage (m)')

  # Plot segmentation asymptotic height
  plot.b.segmented = PlotSegmentation(summary=model_rec$summary.results.segm,
                                      plot_summary=model_rec$plots)$final_plot

  plot.b.segmented[[1]] = plot.b.segmented[[1]]+
    ylab('River bed estimation (m)')

  return(list(plot.rec.segmented=plot.rec.segmented,
              plot.b.segmented=plot.b.segmented))
}

#' Plot segmentation using all stage record
#'
#' @param time vector, time in POSIXct format is mandatory
#' @param obs real vector , stage (m)
#' @param u real vector, uncertainty of the stage
#' @param plot_summary list, plot data resulting from any segmentation function
#'
#' @return ggplot, plot all stage record after segmentation showing homogeneous periods
#' @export
#'
#' @details
#' For an example, please see `?ModelAndSegmentation_Recessions`.
#'
PlotSegmentation_Recessions_Hdt <- function(time,
                                            obs,
                                            u,
                                            plot_summary){
  if(any(is.na(obs) | is.na(u) | is.na(time)))stop('NA value in input data. Be sure to remove them before running function')
  if(is.null(check_vector_lengths(obs,time)))stop('The input data have not the same length')
  if(!lubridate::is.POSIXct(time))stop('time must be POSIXct format')

  data = data.frame(time=time,
                    obs=obs,
                    u=u,
                    I95_lower=obs+stats::qnorm(0.025)*u,
                    I95_upper=obs+stats::qnorm(0.975)*u)

  if(!is.null(plot_summary)){
    shift.info=plot_summary$density.inc.tau


    shift.plot=data.frame(tau=shift.info$taU_MAP,
                          I95_lower=shift.info$tau_lower_inc,
                          I95_upper=shift.info$tau_upper_inc,
                          id_iteration=shift.info$id_iteration)

    shift.plot=shift.plot[order(shift.plot$tau),]

    # Get origin date of the segmentation
    origin.date <- min(data$time)

    # Date transformation function to passe to numeric format if necessary
    DateTransformed <- DateFormatTransform(date=data$time)
    data$time <- DateTransformed$time
    origin.date <- DateTransformed$origin

    shifts.numeric=lubridate::time_length(lubridate::interval(origin.date, shift.plot$tau), "day")

    intervals.time.shift=c(-Inf,shifts.numeric,Inf) # intervals defined by time shifts

    ## cut data depending on the intervals time shift
    data$period <- cut(
      data$time,
      breaks = intervals.time.shift,
      labels = 1:(length(intervals.time.shift)-1),
      right = FALSE
    )

    data$time=NumericFormatTransform(numeric.date = data$time,
                                     origin.date = origin.date)

    plot.data=list(data=data,
                   shift=shift.plot)
  }else{
    # Any shift detected
    plot.data=list(data=data.frame(data,period=1),
                   shift=shift.plot)
  }
  PlotSegmentation(summary=plot.data,
                   plot_summary=plot_summary)$final_plot
}

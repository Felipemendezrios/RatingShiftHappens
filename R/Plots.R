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

#' Plot segmentation
#'
#' Plot the final segmentation displaying shift time along with uncertainties
#'
#' @param summary list, summary data resulting from any segmentation function
#'
#' @return a ggplot
#'
#' @examples
#' # Apply recursive segmentation
#' results=recursive.segmentation(obs=RhoneRiver$H,time=RhoneRiver$Year,u=RhoneRiver$uH)
#'
#' # plot recursive segmentation
#' plotSegmentation(summary=results$summary)
#' @export
plotSegmentation <- function(summary) {

  data=summary$data
  shift=summary$shift

  # Add some colors to the palette for observations
  colourCount_obs = length(unique(data$period))
  getPalette_obs =  scales::viridis_pal(option='D')

  # Add some colors to the palette for shift
  colourCount_tau = length(unique(shift$tau))
  getPalette_tau = scales::viridis_pal(option = "C")

  # Plot shift times
  g=ggplot(data)+
    geom_vline(xintercept = shift$tau,alpha=0.8)+
    coord_cartesian(ylim = c(min(data$I95_lower),max(data$I95_upper)))+
    geom_rect(data = shift,
              aes(xmin = I95_lower,
                  xmax = I95_upper,
                  ymin = -Inf,
                  ymax = Inf,
                  fill=(factor(tau))),
              alpha=0.4)+
    labs(fill='Shift time')

  # Plot observations by period
  g=g+
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
         col='Period',
         title = 'Segmentation')+
    scale_color_manual(values = getPalette_obs(colourCount_obs))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  if(is.numeric(shift$tau)){
    g=g+
      scale_fill_manual(values=getPalette_tau(colourCount_tau),
                        labels=round(shift$tau,2))
  }else{
    g=g+
      scale_fill_manual(values=getPalette_tau(colourCount_tau),
                        labels=round(shift$tau,units='days'))
  }

  return(g)
}

#' Plot Rating curve after segmentation
#'
#' Plot of all the rating curve after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary data.frame, summary data resulting from `recursive.ModelAndSegmentation` function saved in `summary$data`
#' @param equation the equation to be considered for plotting : see ‘Details’
#' @param Hmin_user real value, minimal stage record set by user for plotting in meters
#' @param Hmax_user real value, minimal stage record set by user for plotting in meters
#' @param H_step_discretization real positive non-zero value, time step of the discretization in meters
#' @param autoscale logical, auto scale following data for plotting
#' @param logscale logical, `TRUE`= log scale in discharge axis, `FALSE` = normal scale in discharge axis
#'
#' @return ggplot, rating curve after segmentation
#' @details
#' Some equations for estimating the rating curve are available in this package.
#' Use `GetCatalog()$Equations` to discover the equations supported. More information in `?GetCatalog()`.
#' @export
plotRC_ModelAndSegmentation=function(summary,
                                     equation = LinearRegression_Equation,
                                     Hmin_user = 0,
                                     Hmax_user = 2,
                                     H_step_discretization=0.01,
                                     autoscale=TRUE,
                                     logscale = FALSE){
  # Check discretization
  if(H_step_discretization<=0)stop('The discretization must be a positive non-zero value')

  # Summary data
  summary_data=summary$data
  # summary equation paramaters
  summary_param=summary$param.equation
  # summary shift times adding first and last discharge measurement
  summary_shift=rbind(min(summary_data$time),
                      summary$shift,
                      max(summary_data$time))

  # Check negative and or zero values in discharge
  columns_to_check = c('Q_I95_lower','Q_I95_upper','Qsim','Qsim_I95_lower','Qsim_I95_upper')
  if(any(summary_data<=0)){
    if(logscale==TRUE){
      warning(paste0('Rating curve presented in log scale does not accept negative values and zero. Some values have been replace by NA'))
      DF_remplacement = replace_negatives_or_zero_values(data_frame = summary_data,
                                                         columns = columns_to_check,
                                                         consider_zero = TRUE,
                                                         replace = NA )
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
  colourCount_period = length(unique(summary_data$period))
  getPalette_period =  scales::viridis_pal(option='D')

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

  # Convert parameters from equation in a list of each component is related to rating curve
  param <- split(summary_param,seq_len(nrow(summary_param)))

  # Uncertainty associated to the simulation of discharge by rating curve
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

  # Summary table of all simulation results discretized at H_grid if possible
  summary_sim_grid <- c()
  period_date_begin_obs <- c()
  period_date_end_obs <- c()

  for(i in 1:length(param)){

    # Discharge simulation:
    # Condition for estimating simulation using equation, at least two parameters (a,b) but possibility to integrate in the future equation with only one parameter (a)
    if(!is.na(param[[i]]$a)){

      # BaRatin Case : three parameters (a,b,c)
      if(!is.na(param[[i]]$c)){
        Q_sim_grid_p <- equation(H_grid, a = param[[i]]$a, b = param[[i]]$b, c = param[[i]]$c)

        # Linear and exponential regression: two parameters (a,b)
      }else if(!is.na(param[[i]]$b)){
        Q_sim_grid_p <- equation(H_grid, a = param[[i]]$a, b = param[[i]]$b)

        # Regression with only a parameters (a)
      }else{
        Q_sim_grid_p <- equation(H_grid, a = param[[i]]$a)
      }

      Q_sim_lower_p <- Q_sim_grid_p + uq_lower[i]
      Q_sim_upper_p <- Q_sim_grid_p + uq_upper[i]

      summary_sim_grid_p = data.frame(H_sim=H_grid,
                                      Qsim=Q_sim_grid_p,
                                      Qsim_I95_lower=Q_sim_lower_p,
                                      Qsim_I95_upper=Q_sim_upper_p,
                                      period=i,
                                      period_date_begin=summary_shift$tau[i],
                                      period_date_end=summary_shift$tau[i+1])

      summary_sim_grid=rbind(summary_sim_grid,summary_sim_grid_p)

    }else{
      # When any paramaters has been estimated. Discharge simulation at stage observed
      summary_data_p=summary_data[which(summary_data$period==i),]

      summary_sim_grid_p =  data.frame(H_sim=summary_data_p$H,
                                       Qsim=summary_data_p$Qsim,
                                       Qsim_I95_lower=summary_data_p$Qsim_I95_lower,
                                       Qsim_I95_upper=summary_data_p$Qsim_I95_upper,
                                       period=summary_data_p$period,
                                       period_date_begin=summary_shift$tau[i],
                                       period_date_end=summary_shift$tau[i+1])

      summary_sim_grid=rbind(summary_sim_grid,summary_sim_grid_p)
    }

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
                                                         replace = NA )
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
  labs(x='Stage (m)',
       y='Discharge (m3/s)',
       title = 'Rating curves after segmentation')+
    scale_color_manual(values = getPalette_period(colourCount_period),
                       name = 'Final validity date')+
    scale_fill_manual(values = getPalette_period(colourCount_period),
                      name = 'Period')+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Log scale
  if(logscale==TRUE){
    RC_plot=RC_plot+coord_trans(y='log10',xlim = c(H_min_plot_limit,H_max_plot_limit))
  }

  return(RC_plot)
}

#' Plot shift times
#'
#' Plot shift times after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#'
#' @return ggplot, stage record and shift times
#' @export
plotStage_ModelAndSegmentation <- function(summary){

  data=summary$data
  shift=summary$shift
  # Add some colors to the palette for observations
  colourCount_period = length(unique(data$period))
  getPalette_period =  scales::viridis_pal(option='D')

  # Add some colors to the palette for shift
  colourCount_tau = length(unique(shift$tau))
  getPalette_tau = scales::viridis_pal(option = "C")

  plotStageSegmentation=ggplot(data)
  if(!is.null(shift)){
    plotStageSegmentation=plotStageSegmentation+geom_vline(xintercept = shift$tau,alpha=0.8)+
      geom_rect(data = shift,
                aes(xmin = I95_lower,
                    xmax = I95_upper,
                    ymin = -Inf,
                    ymax = Inf,
                    fill=(factor(tau))),
                alpha=0.4)+
      labs(fill='Shift time')

    if(is.numeric(shift$tau)){
      plotStageSegmentation=plotStageSegmentation+
        scale_fill_manual(values=getPalette_tau(colourCount_tau),
                          labels=round(shift$tau,2))
    }else{
      plotStageSegmentation=plotStageSegmentation+
        scale_fill_manual(values=getPalette_tau(colourCount_tau),
                          labels=round(shift$tau,units='days'))
    }
  }

  plotStageSegmentation=plotStageSegmentation+
    coord_cartesian(ylim = c(min(data$H),max(data$H)))+
    geom_point(aes(x=time,
                   y=H,
                   col=factor(period)),
               show.legend = FALSE)+
    labs(x='Time',
         y='Stage m',
         col='Period',
         title = 'Stage record and segmentation')+
    scale_color_manual(values = getPalette_period(colourCount_period))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  return(plotStageSegmentation)
}

#' Plot residuals with updating rating curve
#'
#' Plot residuals after using `recursive.ModelAndSegmentation` function, along with associated uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#'
#' @return ggplot, residual segmentation
#' @export
plotResidual_ModelAndSegmentation <- function(summary){

  data=summary$data
  shift=summary$shift

  if(is.null(shift))stop('Any shift time detected')

  plotresidual=plotSegmentation(summary = list(data=data.frame(time=data$time,
                                                               obs=data$Qres,
                                                               u=NA,
                                                               I95_lower=(data$Q_I95_lower-data$Qsim),
                                                               I95_upper=(data$Q_I95_upper-data$Qsim),
                                                               period=data$period),
                                               shift=shift))

  plotresidual=plotresidual+
    guides(col='none')+
    ylab('Residual (m3/s)')

  return(plotresidual)
}

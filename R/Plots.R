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
#' @param summary list, summary data resulting from `recursive.ModelAndSegmentation` function
#' @param logscale logical, `TRUE`= log scale, `FALSE` = normal scale
#' @return ggplot, rating curve after segmentation
#' @export
plotRC_ModelAndSegmentation <- function(summary,logscale=FALSE) {

  data=summary$data

  # Check negative and or zero values in discharge
  columns_to_check = c('Q_I95_lower','Q_I95_upper','Qsim','Qsim_I95_lower','Qsim_I95_upper')
  if(any(data<=0)){
    if(logscale==TRUE){
      warning(paste0('Rating curve presented in log scale does not accept negative values and zero. Some values have been replace by NA'))
      DF_remplacement = replace_negatives_or_zero_values(data_frame = data,
                                                         columns = columns_to_check,
                                                         consider_zero = TRUE,
                                                         replace = NA )
    }else{
      warning(paste0('Rating curve does not accept negative values zero. Some values have been replace by zero'))
      DF_remplacement = replace_negatives_or_zero_values(data_frame = data,
                                                         columns = columns_to_check,
                                                         consider_zero = FALSE,
                                                         replace = 0 )
    }
    data[,match(columns_to_check,colnames(data))]=DF_remplacement
  }

  # Add some colors to the palette for observations
  colourCount_period = length(unique(data$period))
  getPalette_period =  scales::viridis_pal(option='D')

  # Plot RC by period
  plotRC=ggplot(data)+
    geom_ribbon(aes(x=H,
                    ymin=Qsim_I95_lower,
                    ymax=Qsim_I95_upper,
                    group=factor(period),
                    fill=factor(period)),
                alpha=0.5,
                na.rm = FALSE,
                show.legend = FALSE)+
    geom_line(aes(x=H,
                  y=Qsim,
                  group=factor(period),
                  col=factor(period)),
              alpha=0.7)+
    geom_errorbar(aes(x=H,
                      y=Q,
                      ymin=Q_I95_lower,
                      ymax=Q_I95_upper,
                      group=factor(period),
                      col=factor(period)),
                  width=0.1)+
    geom_point(aes(x=H,
                   y=Q,
                   group=factor(period),
                   col=factor(period)))+
    labs(x='Stage (m)',
         y='Discharge m3/s',
         col='Period',
         title = 'Rating curve and segmentation')+
    scale_color_manual(values = getPalette_period(colourCount_period))+
    scale_fill_manual(values = getPalette_period(colourCount_period))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  if(logscale==TRUE){
    plotRC=plotRC+coord_trans(y='log10')
  }
  return(plotRC)
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

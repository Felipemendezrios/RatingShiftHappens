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
#' @import  ggplot2
plotSegmentation <- function(summary) {

  data=summary$data
  shift=summary$shift

  # Add some colors to the palette for observations
  colourCount_obs = length(unique(data$period))
  getPalette_obs =  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

  # Plot observations by period
  g=ggplot(data)+
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

  # Add some colors to the palette for observations
  colourCount_tau = length(unique(shift$tau))
  getPalette_tau = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))
  # Plot shift times
  g=g+
    geom_vline(xintercept = shift$tau,alpha=0.8)+
    coord_cartesian(ylim = c(min(data$I95_lower),max(data$I95_upper)))+
    geom_rect(data = shift,
              aes(xmin = I95_lower,
                  xmax = I95_upper,
                  ymin = -Inf,
                  ymax = Inf,
                  fill=(factor(tau))),
              alpha=0.6)+
    labs(fill='Shift time')

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

#' Plot segmentation with the updated rating curve
#'
#' Plot the final segmentation with the updated rating curve,
#' displaying shift time along with uncertainties
#'
#' @param summary list, summary data resulting from model and segmentation function
#'
#' @return a ggplot
#'
#' @examples
#' # Apply recursive segmentation
#' results=recursive.ModelAndSegmentation(H=RhoneRiver$H,
#'                                        Q=RhoneRiver$Q,
#'                                        time=RhoneRiver$Year,
#'                                        uQ=RhoneRiver$uH)
#'
#' # plot recursive segmentation
#' plotModelAndSegmentation(summary=results$summary)
#' @export
#' @import  ggplot2
plotModelAndSegmentation <- function(summary) {

  data=summary$data
  shift=summary$shift

  # Add some colors to the palette for observations
  colourCount_obs = length(unique(data$period))
  getPalette_obs =  grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

  # Plot observations by period
  plotRC=ggplot(data)+
    geom_point(aes(x=H,
                   y=Q,
                   col=factor(period)))+
    geom_errorbar(aes(x=H,
                      y=Q,
                      ymin=I95_lower,
                      ymax=I95_upper,
                      col=factor(period)),
                  width=0.1)+
    labs(x='Stage (m)',
         y='Discharge m3/s',
         col='Period',
         title = 'Rating curve and segmentation')+
    scale_color_manual(values = getPalette_obs(colourCount_obs))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Add some colors to the palette for observations
  colourCount_tau = length(unique(shift$tau))
  getPalette_tau = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Pastel1"))

  # Plot observations by period
  PlotStageSegmentation=ggplot(data)+
    geom_point(aes(x=time,
                   y=H,
                   col=factor(period)))+
    labs(x='Time',
         y='Stage m',
         col='Period',
         title = 'Stage record and segmentation')+
    scale_color_manual(values = getPalette_obs(colourCount_obs))+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Plot shift times
  PlotStageSegmentation=PlotStageSegmentation+
    geom_vline(xintercept = shift$tau,alpha=0.8)+
    coord_cartesian(ylim = c(min(data$H),max(data$H)))+
    geom_rect(data = shift,
              aes(xmin = I95_lower,
                  xmax = I95_upper,
                  ymin = -Inf,
                  ymax = Inf,
                  fill=(factor(tau))),
              alpha=0.6)+
    labs(fill='Shift time')

  if(is.numeric(shift$tau)){
    PlotStageSegmentation=PlotStageSegmentation+
      scale_fill_manual(values=getPalette_tau(colourCount_tau),
                        labels=round(shift$tau,2))
  }else{
    PlotStageSegmentation=PlotStageSegmentation+
      scale_fill_manual(values=getPalette_tau(colourCount_tau),
                        labels=round(shift$tau,units='days'))
  }

  return(list(plotRC=plotRC,
              PlotStageSegmentation=PlotStageSegmentation))
}

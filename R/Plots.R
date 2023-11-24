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
  DF$y=-1*DF$level
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

  # plot
  g=ggplot(DF)+
    geom_segment(aes(x=xstart,y=ystart,xend=xend,yend=yend),
                 arrow=arrow(type='closed',length=unit(0.03, "npc")))+
    geom_label(aes(x,y,label=indx,fill=isTerminal),size=5,color='white')+
    scale_fill_manual(values=c('royalblue3','red3'))+
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
    scale_color_brewer(palette = 'Dark2')+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,
                                    face='bold',
                                    size=15),
          legend.title.align=0.5)

  # Plot shift times
  g=g+
    geom_vline(xintercept = shift$tau,alpha=0.8)+
    coord_cartesian(ylim = c(min(data$obs),max(data$obs)))+
    geom_rect(data = shift,
              aes(xmin = I95_lower,
                  xmax = I95_upper,
                  ymin = -Inf,
                  ymax = Inf,
                  fill=(factor(tau))),
              alpha=0.4)+
    labs(fill='Shift time')+
    scale_fill_brewer(palette='Pastel1',
                      labels=round(shift$tau,2))

  return(g)
}

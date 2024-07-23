

#### Choose the best node to transfer posterior information to the rest!!!
library(RatingShiftHappens)

atej=ArdecheRiverMeyrasGaugings[,c('Date','H','Q')]

# Three segments identified
atej$id=c(rep(1,40),rep(2,50),rep(3,61))

set.seed(123)

# Dataset following  theoretical CDF
sample1 <- seq(min(atej$H),max(atej$H),0.01)
# Dataset following all gaugings
sample2 <- atej$H

ecdf1 <- ecdf(sample1)
ecdf2 <- ecdf(sample2)

plot(ecdf1, col = "blue", main = "Comparaison des fonctions de répartition")
plot(ecdf2, col = "red", add = TRUE)

legend("bottomright",
       legend = c("All gaugings",
                  "Theorical CDF"),
       col = c("red",
               "blue"), lty = 1)

plot(ecdf1, col = "blue", main = "Comparaison des fonctions de répartition")
plot(ecdf2, col = "red", add = TRUE)

# groups
group1= atej$H[which(atej$id==1)]
ecdfgroup1 <- ecdf(group1)
plot(ecdfgroup1, col = "yellow", add = TRUE)

group2= atej$H[which(atej$id==2)]
ecdfgroup2 <- ecdf(group2)
plot(ecdfgroup2, col = "brown", add = TRUE)

group3= atej$H[which(atej$id==3)]
ecdfgroup3 <- ecdf(group3)
plot(ecdfgroup3, col = "green", add = TRUE)

legend("bottomright",
       legend = c('All gaugings',
                  "Theorical CDF",
                  "Group 1",
                  "Group 2",
                  "Group 3"),
       col = c("red",
               "blue",
               "yellow",
               "brown",
               "green"), lty = 1)


# Tests crash:

plot(ecdf1, col = "blue", main = "Comparaison des fonctions de répartition")
plot(ecdf2, col = "red", add = TRUE)

# One value
crash_test1=mean(atej$H)*1.5

ecdfct1 <- ecdf(crash_test1)
plot(ecdfct1, col = "black", add = TRUE)
# Several extreme values
crash_test3=c(min(atej$H)*1.01,
              min(atej$H)*1.02,
              min(atej$H)*1.03,
              min(atej$H)*1.04,
              min(atej$H)*1.05,
              max(atej$H)*1.01,
              max(atej$H)*1.02,
              max(atej$H)*1.03,
              max(atej$H)*1.04,
              max(atej$H)*1.05)

ecdfct3 <- ecdf(crash_test3)
plot(ecdfct3, col = "pink", add = TRUE)

# Few extreme values
crash_test2=c(min(atej$H),max(atej$H))

ecdfct2 <- ecdf(crash_test2)
plot(ecdfct2, col = "gray", add = TRUE)

# Several mean values
crash_test5=c(mean(atej$H)*1.01,
              mean(atej$H)*1.02,
              mean(atej$H)*1.03,
              mean(atej$H)*1.04,
              mean(atej$H)*1.05,
              median(atej$H)*1.01,
              median(atej$H)*1.02,
              median(atej$H)*1.03,
              median(atej$H)*1.04,
              median(atej$H)*1.05)

ecdfct5 <- ecdf(crash_test5)

plot(ecdfct5, col = "purple", add = TRUE)

# Few mean values
crash_test4=c(mean(atej$H),median(atej$H))

ecdfct4 <- ecdf(crash_test4)

plot(ecdfct4, col = "orange", add = TRUE)




legend("bottomright",
       legend = c("All gaugings",
                  "Theorical CDF",
                  "Crash Test1",
                  "Crash Test2",
                  "Crash Test3",
                  "Crash Test4"),
       col = c("red",
               "blue",
               'black',
               'gray',
               "pink",
               'orange'), lty = 1)


# ALL
plot(ecdf1, col = "blue", main = "Comparaison des fonctions de répartition")
plot(ecdf2, col = "red", add = TRUE)

plot(ecdfgroup1, col = "yellow", add = TRUE)
plot(ecdfgroup2, col = "brown", add = TRUE)
plot(ecdfgroup3, col = "green", add = TRUE)

plot(ecdfct1, col = "black", add = TRUE)
plot(ecdfct3, col = "pink", add = TRUE)
plot(ecdfct2, col = "gray", add = TRUE)
plot(ecdfct5, col = "purple", add = TRUE)
plot(ecdfct4, col = "orange", add = TRUE)

legend("bottomright",
       legend = c("All gaugings",
                  "Theorical CDF",
                  "Group 1",
                  "Group 2",
                  "Group 3",
                  "Crash Test1",
                  "Crash Test2",
                  "Crash Test3",
                  "Crash Test4"),
       col = c("red",
               "blue",
               "yellow",
               "brown",
               "green",
               'black',
               'gray',
               "pink",
               'orange'), lty = 1)



ks.test.T=data.frame(G1=ks.test(sample1, group1)$statistic,
           G2=ks.test(sample1, group2)$statistic,
           G3=ks.test(sample1, group3)$statistic,
           CT1=ks.test(sample1, crash_test1)$statistic,
           CT2=ks.test(sample1, crash_test2)$statistic,
           CT3=ks.test(sample1, crash_test3)$statistic,
           CT4=ks.test(sample1, crash_test4)$statistic,
           Theoricals=ks.test(sample1, sample2)$statistic)

ks.test.all=data.frame(G1=ks.test(sample2, group1)$statistic,
           G2=ks.test(sample2, group2)$statistic,
           G3=ks.test(sample2, group3)$statistic,
           CT1=ks.test(sample2, crash_test1)$statistic,
           CT2=ks.test(sample2, crash_test2)$statistic,
           CT3=ks.test(sample2, crash_test3)$statistic,
           CT4=ks.test(sample2, crash_test4)$statistic,
           Theoricals=ks.test(sample2, sample1)$statistic)

comparison=rbind(ks.test.T,ks.test.all)

rownames(comparison) <- c('Theorical CDF','All gauging CDF')

comparison

library(ggplot2)
ggplot()+geom_boxplot(data=atej,aes(x=factor(id),y=H))+
  geom_boxplot(data=data.frame(theorical=sample1,id=rep("Theorical",length(sample1))),aes(x=id,y=theorical))+
  geom_boxplot(data=data.frame(All=atej$H,id=rep("All gaugings",nrow(atej))),aes(x=id,y=All))+
  geom_boxplot(data=data.frame(All=crash_test1,id=rep('crash_test1',length(crash_test1))),aes(x=id,y=All))+
  geom_boxplot(data=data.frame(All=crash_test2,id=rep("crash_test2",length(crash_test2))),aes(x=id,y=All))+
  geom_boxplot(data=data.frame(All=crash_test3,id=rep("crash_test3",length(crash_test3))),aes(x=id,y=All))+
  geom_boxplot(data=data.frame(All=crash_test4,id=rep("crash_test4",length(crash_test4))),aes(x=id,y=All))



#### FUNCTIONS COULD BE INTERESTING

#' Function to split the data frame at specified positions
#'
#' @param df data frame that will be split
#' @param positions values, positions indicating the positions for splitting the data frame
#'
#' @return list, all recession indexed by the number of the recession
#' @export
#'
#' @examples
#' df <- data.frame(ID = 1:10,
#'                  Value = c(5, -3, 4, -2, 6, -1, 3, -4, 2, 7)
#'                  )
#'
#' positions <- c(3, 6, 9)
#'
#' split_data_frame(df, positions)
split_data_frame <- function(df, positions) {
  # Add the start and end positions
  positions <- c(0, positions, nrow(df))

  # Split the data frame
  splits <- lapply(seq(length(positions) - 1), function(i) {
    df[(positions[i] + 1):positions[i + 1], ]
  })

  return(splits)
}



# cas d'étude sur le Rhône à Beaucaire
a1=RBaM::parameter(name='a1',init=152.735,prior.dist='LogNormal',prior.par=c(5.029,3.387))
b1=RBaM::parameter(name='b1',init=-2,prior.dist='Gaussian',prior.par=c(-2,0.5))
c1=RBaM::parameter(name='c1',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
a2=RBaM::parameter(name='a2',init=290.241,prior.dist='LogNormal',prior.par=c(5.671,4.019))
b2=RBaM::parameter(name='b2',init=2,prior.dist='Gaussian',prior.par=c(2,0.3))
c2=RBaM::parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))

# Hydraulic control matrix
controlMatrix=matrix(c(1,1,0,1),ncol=2,nrow=2)

a.object=list(a1,a2)
b.object=list(b1,b2)
c.object=list(c1,c2)

resultsBaRatin=recursive.ModelAndSegmentation(H=RhoneRiver$H,
                                              Q=RhoneRiver$Q,
                                              time=RhoneRiver$Date,
                                              uQ=RhoneRiver$uQ,
                                              nSmax=3,
                                              nMin=2,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(RhoneRiver$H),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)

plotTree(resultsBaRatin$tree)

terminal = resultsBaRatin$tree$indx[which(resultsBaRatin$tree$nS==1)]
terminal

PlotRCPrediction(Hgrid=data.frame(seq(0,12,by=0.01)),
                 autoscale=FALSE,
                 temp.folder=file.path(tempdir(),'BaM'),
                 CalibrationData='CalibrationData.txt',
                 allnodes=FALSE,
                 nodes=terminal)

plot_Q_ModelAndSegmentation(summary=resultsBaRatin$summary,
                            plot_summary=resultsBaRatin$plot)$final_plot

plotResidual_ModelAndSegmentation(summary=resultsBaRatin$summary,
                                  plot_summary=resultsBaRatin$plot)$final_plot


Recessions=Extraction_recession(H=,
                                uH=0,
                                time=,
                                chi=200,
                                tgood=30,
                                delta.t.min = 0,
                                delta.t.max = 50)

plot_rec_extracted(Rec_extracted = Recessions)

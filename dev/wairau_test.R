

test=recursive.ModelAndSegmentation(H=WairauRiverGaugings$H,
                               Q=WairauRiverGaugings$Q,
                               time=WairauRiverGaugings$Date,
                               uQ=WairauRiverGaugings$uQ)

# Plot shift times in stage record
plot_H_ModelAndSegmentation(summary=test$summary,
                            plot_summary=test$plot)

# Plot shift times in discharge observations
plot_Q_ModelAndSegmentation(summary=test$summary,
                            plot_summary=test$plot)

# Plot residual
plotResidual_ModelAndSegmentation(summary=test$summary,
                                  plot_summary=test$plot)

a=test$summary$param.equation$a
b=test$summary$param.equation$b

plotRC_ModelAndSegmentation(summary=test$summary,
                            equation = Exponential_Equation,
                            a=a,
                            b=b)




plotRC_ModelAndSegmentation(summary=test$summary,
                            equation = Exponential_Equation,
                            a=a,
                            b=b,
                            autoscale = FALSE,
                            Hmin_user = 2,
                            Hmax_user = 3.5,
                            H_step_discretization = 0.01)

fit=fitRC_BaRatin

# Hydraulic matrix control is also required.
# The `control_matrix_builder` was developed to help the user to create this matrix.
# In this case, 3 hydraulic control has been set

controlMatrix=matrix(c(1,0,0,1),ncol=2,nrow=2)

# Prior information about input data is required.
# The `prior_infor_param_builder` function was developed to help the user to create these objects.
#' # Prior information for Ardeche River at Meyras

a1=RBaM::parameter(name='a1',init=47.27,prior.dist='LogNormal',prior.par=c(3.86,0.2))
b1=RBaM::parameter(name='b1',init=2,prior.dist='Gaussian',prior.par=c(2,0.025))
c1=RBaM::parameter(name='c1',init=5/3,prior.dist='Gaussian',prior.par=c(5/3,0.025))

a2=RBaM::parameter(name='a2',init=189.08,prior.dist='LogNormal',prior.par=c(5.24217,0.3))
b2=RBaM::parameter(name='b2',init=2.5,prior.dist='Gaussian',prior.par=c(2.5,0.5))
c2=RBaM::parameter(name='c2',init=5/3,prior.dist='Gaussian',prior.par=c(5/3,0.025))

a.object=list(a1,a2)
b.object=list(b1,b2)
c.object=list(c1,c2)

resultsBaRatin=recursive.ModelAndSegmentation(H=WairauRiverGaugings$H,
                                              Q=WairauRiverGaugings$Q,
                                              time=WairauRiverGaugings$Date,
                                              uQ=WairauRiverGaugings$uQ,
                                              nSmax=3,
                                              nMin=2,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(WairauRiverGaugings$H),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
                                              )

plotTree(resultsBaRatin$tree)

# Terminal nodes
terminal = resultsBaRatin$tree$indx[which(resultsBaRatin$tree$nS==1)]
terminal

# Plot the rating curves after using BaRatin method. Function specially created to this method
PlotRCPrediction(Hgrid=data.frame(seq(1.5,6,by=0.01)),
                 autoscale=FALSE,
                 temp.folder=file.path(tempdir(),'BaM'),
                 CalibrationData='CalibrationData.txt',
                 allnodes=FALSE,
                 nodes=terminal)

# Plot shift times in stage record
plot_H_ModelAndSegmentation(summary=resultsBaRatin$summary,
                            plot_summary=resultsBaRatin$plot)

# Plot shift times in discharge observations
plot_Q_ModelAndSegmentation(summary=resultsBaRatin$summary,
                            plot_summary=resultsBaRatin$plot)

# Plot residual
plotResidual_ModelAndSegmentation(summary=resultsBaRatin$summary,
                                  plot_summary=resultsBaRatin$plot)

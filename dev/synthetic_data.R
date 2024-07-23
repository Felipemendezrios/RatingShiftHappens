# CLASS 1 :
Class1=synthetic_gauging_datasets[[1]][[1]]
Class1_shifttimes=synthetic_gauging_datasets[[1]][[2]]

controlMatrix=matrix(c(1),ncol=1,nrow=1)

a1=RBaM::parameter(name='a1',init=25,prior.dist='Gaussian',prior.par=c(25,5))
b1=RBaM::parameter(name='b1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.6,2))
c1=RBaM::parameter(name='c1',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.005))

a.object=list(a1)
b.object=list(b1)
c.object=list(c1)


results_class1=recursive.ModelAndSegmentation(H=Class1$h,
                                              Q=Class1$Q,
                                              time=Class1$time,
                                              uQ=Class1$uQ,
                                              nSmax=3,
                                              nMin=1,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(Class1$h),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)


# CLASS 3 :
Class3=synthetic_gauging_datasets[[2]][[1]]
Class3_shifttimes=synthetic_gauging_datasets[[2]][[2]]


results_class3=recursive.ModelAndSegmentation(H=Class3$h,
                                              Q=Class3$Q,
                                              time=Class3$time,
                                              uQ=Class3$uQ,
                                              nSmax=3,
                                              nMin=1,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(Class1$h),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)

# CLASS 4 :
Class4=synthetic_gauging_datasets[[3]][[1]]
Class4_shifttimes=synthetic_gauging_datasets[[3]][[2]]


results_Class4=recursive.ModelAndSegmentation(H=Class4$h,
                                              Q=Class4$Q,
                                              time=Class4$time,
                                              uQ=Class4$uQ,
                                              nSmax=3,
                                              nMin=1,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(Class1$h),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)

# CLASS 6 :
Class6=synthetic_gauging_datasets[[4]][[1]]
Class6_shifttimes=synthetic_gauging_datasets[[4]][[2]]


results_Class6=recursive.ModelAndSegmentation(H=Class6$h,
                                              Q=Class6$Q,
                                              time=Class6$time,
                                              uQ=Class6$uQ,
                                              nSmax=3,
                                              nMin=1,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(Class1$h),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)

# CLASS 9 :
Class9=synthetic_gauging_datasets[[5]][[1]]
Class9_shifttimes=synthetic_gauging_datasets[[5]][[2]]


results_Class9=recursive.ModelAndSegmentation(H=Class9$h,
                                              Q=Class9$Q,
                                              time=Class9$time,
                                              uQ=Class9$uQ,
                                              nSmax=3,
                                              nMin=1,
                                              funk=fitRC_BaRatin,
                                              HmaxGrid=max(Class1$h),
                                              a.object=a.object,
                                              b.object=b.object,
                                              c.object=c.object,
                                              controlMatrix=controlMatrix
)

# CLASS 10 :
Class10=synthetic_gauging_datasets[[6]][[1]]
Class10_shifttimes=synthetic_gauging_datasets[[6]][[2]]

a1_10=RBaM::parameter(name='a1',init=14,prior.dist='Gaussian',prior.par=c(14,5))
b1_10=RBaM::parameter(name='b1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.6,2))
c1_10=RBaM::parameter(name='c1',init=1.5,prior.dist='Gaussian',prior.par=c(1.5,0.005))

a2_10=RBaM::parameter(name='a2',init=26,prior.dist='Gaussian',prior.par=c(26,5))
b2_10=RBaM::parameter(name='b2',init=-0.4,prior.dist='Gaussian',prior.par=c(-0.4,2))
c2_10=RBaM::parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.005))

a3_10=RBaM::parameter(name='a3',init=40,prior.dist='Gaussian',prior.par=c(40,10))
b3_10=RBaM::parameter(name='b3',init=0.2,prior.dist='Gaussian',prior.par=c(0.2,3))
c3_10=RBaM::parameter(name='c3',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.005))

a10.object=list(a1_10,a2_10,a3_10)
b10.object=list(b1_10,b2_10,b3_10)
c10.object=list(c1_10,c2_10,c3_10)

controlMatrix10=matrix(c(1,0,0,0,1,1,0,0,1),ncol=3,nrow=3)

results_Class10=recursive.ModelAndSegmentation(H=Class10$h,
                                               Q=Class10$Q,
                                               time=Class10$time,
                                               uQ=Class10$uQ,
                                               nSmax=3,
                                               nMin=1,
                                               funk=fitRC_BaRatin,
                                               HmaxGrid=max(Class1$h),
                                               a.object=a10.object,
                                               b.object=b10.object,
                                               c.object=c10.object,
                                               controlMatrix=controlMatrix10
)

results_class1$summary$shift
Class1_shifttimes

results_class3$summary$shift
Class3_shifttimes

results_Class4$summary$shift
Class4_shifttimes

results_Class6$summary$shift
Class6_shifttimes

results_Class9$summary$shift
Class9_shifttimes

results_Class10$summary$shift
Class10_shifttimes



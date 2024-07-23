H=ArdecheRiverMeyrasGaugings$H
Q=ArdecheRiverMeyrasGaugings$Q
# obs=ArdecheRiverMeyrasGaugings$H
# obs=RhoneRiver$H
# obs=CongoRiverBrazzavilleMINAN$Q
# Q=CongoRiverBrazzavilleMINAN$Q
# H=CongoRiverBrazzavilleMINAN$H
time=ArdecheRiverMeyrasGaugings$Date
# time=RhoneRiver$Year
# time=CongoRiverBrazzavilleMINAN$time
uQ=ArdecheRiverMeyrasGaugings$uQ
# u=CongoRiverBrazzavilleMINAN$uQ
# u=ArdecheRiverMeyrasGaugings$uQ
# u=RhoneRiver$uH
# uQ=CongoRiverBrazzavilleMINAN$uQ
nS=2
nSmax=3
nMin= 2
nCycles=100
burn=0.5
nSlim=max(nCycles/10,1)
temp.folder=file.path(tempdir(),'BaM')
# funk=fitRC_SimplifiedBaRatin
funk=fitRC_BaRatin


results=recursive.ModelAndSegmentation(H=H,
                               Q=Q,
                               time=time,
                               uQ=uQ,
                               nSmax=nSmax,nMin=nMin,
                               funk=fitRC_SimplifiedBaRatin,
                               Hgrid=Hgrid,
                               temp.folder = file.path(tempdir(),'BaM'))



controlMatrix=matrix(c(1,0,0,0,1,1,0,0,1),ncol=3,nrow=3)

a1=parameter(name='a1',init=14.17,prior.dist='LogNormal',prior.par=c(2.66,1.54))
b1=parameter(name='b1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.58,1.49))
c1=parameter(name='c1',init=1.5,prior.dist='Gaussian',prior.par=c(1.5,0.025))
a2=parameter(name='a2',init=26.5165,prior.dist='LogNormal',prior.par=c(3.28,0.36))
b2=parameter(name='b2',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.58,1.49))
c2=parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
a3=parameter(name='a3',init=31.82,prior.dist='LogNormal',prior.par=c(3.46,0.397))
b3=parameter(name='b3',init=1.2,prior.dist='Gaussian',prior.par=c(1.2,0.2))
c3=parameter(name='c3',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))

a.object=list(a1,a2,a3)
b.object=list(b1,b2,b3)
c.object=list(c1,c2,c3)


# Recession

H=ArdecheRiverMeyrasStage$H
time=ArdecheRiverMeyrasStage$Date
chi=1.5
delta.t.min=0
delta.t.max=5
tgood=20
Nmin=3

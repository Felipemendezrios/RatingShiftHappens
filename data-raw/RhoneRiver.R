RhoneRiver=read.table('data-raw/RhoneRiver.txt',header=T,sep=';')
RhoneRiver$Date=as.Date(RhoneRiver$Date)
save(RhoneRiver,file='data/RhoneRiver.RData')

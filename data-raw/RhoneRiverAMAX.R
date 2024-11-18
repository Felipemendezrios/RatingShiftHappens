RhoneRiverAMAX=read.table('data-raw/RhoneRiverAMAX.txt',header=T,sep=';')
RhoneRiverAMAX$Date=as.Date(RhoneRiverAMAX$Date)
save(RhoneRiverAMAX,file='data/RhoneRiverAMAX.RData')

RhoneRiverGaugings=read.table('data-raw/RhoneRiverGaugings.txt',header=T,sep=';')
RhoneRiverGaugings$Date=as.Date(RhoneRiverGaugings$Date)
save(RhoneRiverGaugings,file='data/RhoneRiverGaugings.RData')

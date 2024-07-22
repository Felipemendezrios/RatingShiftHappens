
RhoneRiverStage=utils::read.csv('data-raw/PtBcrJournalier.csv',
                        header = T, sep = ';',  colClasses = c('Date','numeric'))

RhoneRiverStage=na.omit(RhoneRiverStage)

RhoneRiverStage$Date <- as.POSIXct(RhoneRiverStage$Date)

save(RhoneRiverStage,file='data/RhoneRiverStage.RData')

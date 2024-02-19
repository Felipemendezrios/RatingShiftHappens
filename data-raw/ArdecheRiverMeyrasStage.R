ArdecheRiverMeyrasStageRaw=read.table('data-raw/ArdecheRiverMeyrasStage.txt',header=T,sep='\t')
colnames(ArdecheRiverMeyrasStageRaw) <- c('NumericDate','H','Date','DateFromZero')

DateString<- paste(ArdecheRiverMeyrasStageRaw$Date,rep(0,length(ArdecheRiverMeyrasStageRaw$Date)),sep = ':')
ArdecheRiverMeyrasStage <- data.frame(Date = as.POSIXct(DateString, format = '%d/%m/%Y %H:%M:%S'),
                                      H = ArdecheRiverMeyrasStageRaw$H
                                      )
save(ArdecheRiverMeyrasStage,file='data/ArdecheRiverMeyrasStage.RData')

ArdecheRiverStage=read.table('data-raw/ArdecheRiverStage.txt',header=T,sep='\t')
colnames(ArdecheRiverStage) <- c('NumericDate','H','Date','DateFromZero')

DateString<- paste(ArdecheRiverStage$Date,rep(0,length(ArdecheRiverStage$Date)),sep = ':')
ArdecheRiverStage <- data.frame(Date = as.POSIXct(DateString, format = '%d/%m/%Y %H:%M:%S', tz='UTC'),
                                      H = ArdecheRiverStage$H/100 # Stage record in meter
                                      )
save(ArdecheRiverStage,file='data/ArdecheRiverStage.RData')

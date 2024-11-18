ArdecheRiverGaugingsRaw=read.table('data-raw/ArdecheRiverGaugings.txt',header=T,sep='\t')

ArdecheRiverGaugingsTemp=ArdecheRiverGaugingsRaw[,1:(length(ArdecheRiverGaugingsRaw)-2)]
colnames(ArdecheRiverGaugingsTemp)=c('Day','Month','Year','Hour','Minute','H','Q','Qsigma')
ArdecheRiverGaugingsDate=data.frame(ArdecheRiverGaugingsTemp[,1:5],
                          Second=rep(0,nrow(ArdecheRiverGaugingsTemp)))

DateString <- paste(ArdecheRiverGaugingsDate$Day, ArdecheRiverGaugingsDate$Month, ArdecheRiverGaugingsDate$Year, sep = "/")
TimeString <- paste(ArdecheRiverGaugingsDate$Hour, ArdecheRiverGaugingsDate$Minute, ArdecheRiverGaugingsDate$Second, sep = ":")
DateTimeString <- paste(DateString, TimeString, sep = " ")

ArdecheRiverGaugingsDate$Date <-  as.POSIXct(DateTimeString, format = '%d/%m/%Y %H:%M:%S' , tz='UTC')

ArdecheRiverGaugings=data.frame(ArdecheRiverGaugingsDate,
                                        H=ArdecheRiverGaugingsTemp$H,
                                        Q=ArdecheRiverGaugingsTemp$Q,
                                        uQ=round(ArdecheRiverGaugingsTemp$Q*ArdecheRiverGaugingsTemp$Qsigma/100/2,4) # to make sens with other file of gaugings at Ardecher river at Meyras
                                        )
save(ArdecheRiverGaugings,file='data/ArdecheRiverGaugings.RData')

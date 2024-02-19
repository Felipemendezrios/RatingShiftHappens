ArdecheRiverMeyrasGaugingsRaw=read.table('data-raw/ArdecheRiverMeyrasGaugings.txt',header=T,sep='\t')

ArdecheRiverMeyrasGaugingsTemp=ArdecheRiverMeyrasGaugingsRaw[,1:(length(ArdecheRiverMeyrasGaugingsRaw)-2)]
colnames(ArdecheRiverMeyrasGaugingsTemp)=c('Day','Month','Year','Hour','Minute','H','Q','Qsigma')
ArdecheRiverMeyrasGaugingsDate=data.frame(ArdecheRiverMeyrasGaugingsTemp[,1:5],
                          Second=rep(0,nrow(ArdecheRiverMeyrasGaugingsTemp)))

DateString <- paste(ArdecheRiverMeyrasGaugingsDate$Day, ArdecheRiverMeyrasGaugingsDate$Month, ArdecheRiverMeyrasGaugingsDate$Year, sep = "/")
TimeString <- paste(ArdecheRiverMeyrasGaugingsDate$Hour, ArdecheRiverMeyrasGaugingsDate$Minute, ArdecheRiverMeyrasGaugingsDate$Second, sep = ":")
DateTimeString <- paste(DateString, TimeString, sep = " ")

ArdecheRiverMeyrasGaugingsDate$Date <-  as.POSIXct(DateTimeString, format = '%d/%m/%Y %H:%M:%S')

ArdecheRiverMeyrasGaugingsDF=data.frame(ArdecheRiverMeyrasGaugingsDate,
                                        H=ArdecheRiverMeyrasGaugingsTemp$H,
                                        Q=ArdecheRiverMeyrasGaugingsTemp$Q,
                                        uQ=ArdecheRiverMeyrasGaugingsTemp$Q*ArdecheRiverMeyrasGaugingsTemp$Qsigma/100
                                        )
save(ArdecheRiverMeyrasGaugingsDF,file='data/ArdecheRiverMeyrasGaugingsDF.RData')

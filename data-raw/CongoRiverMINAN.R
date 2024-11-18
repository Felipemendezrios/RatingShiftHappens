CongoRiverMINANRaw=read.table('data-raw/Congo_hydro_1902-2024.csv',
                                 header=TRUE,
                                 quote='\"',
                                 sep=',',
                                 na.strings='NaN',
                                 colClasses = c('POSIXct',rep('numeric',5)))

colnames(CongoRiverMINANRaw) <- sub("\\..*","",colnames(CongoRiverMINANRaw))

CongoRiverMINANRaw$Time <- as.POSIXct(CongoRiverMINANRaw$Time,tz='CET')

uQ=(CongoRiverMINANRaw$Q_total_high-CongoRiverMINANRaw$Q_total_low)/(2*1.96)
date=CongoRiverMINANRaw$Time

year=as.numeric(substring(date,1,4))

dataset_Congo=data.frame(time=date,
                         year=year,
                         Q=CongoRiverMINANRaw$Q_maxpost,
                         uQ=uQ)

# The MINAN Congo river at Brazzaville
CongoRiverMINAN <- convert_list_to_dataframe(
  lapply(split(dataset_Congo,
               dataset_Congo$year),
         FUN = function(data){
           min_index <- which.min(data$Q)
           return(data[min_index,])
         }
         )
  )

save(CongoRiverMINAN,file='data/CongoRiverMINAN.RData')

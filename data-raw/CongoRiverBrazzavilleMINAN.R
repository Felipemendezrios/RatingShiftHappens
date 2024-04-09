CongoRiverBrazzavilleRaw=read.table('data-raw/Congo_hydro_1902-2024.csv',
                                 header=TRUE,
                                 quote='\"',
                                 sep=',',
                                 na.strings='NaN',
                                 colClasses = c('POSIXct',rep('numeric',5)))

colnames(CongoRiverBrazzavilleRaw) <- sub("\\..*","",colnames(CongoRiverBrazzavilleRaw))

CongoRiverBrazzavilleRaw$Time <- as.POSIXct(CongoRiverBrazzavilleRaw$Time,tz='CET')

uQ=(CongoRiverBrazzavilleRaw$Q_total_high-CongoRiverBrazzavilleRaw$Q_total_low)/(2*1.96)
date=CongoRiverBrazzavilleRaw$Time

year=as.numeric(substring(date,1,4))

dataset_Congo=data.frame(time=date,
                         year=year,
                         Q=CongoRiverBrazzavilleRaw$Q_maxpost,
                         uQ=uQ)

# The MINAN Congo river at Brazzaville
CongoRiverBrazzavilleMINAN <- convert_list_to_dataframe(
  lapply(split(dataset_Congo,
               dataset_Congo$year),
         FUN = function(data){
           min_index <- which.min(data$Q)
           return(data[min_index,])
         }
         )
  )

save(CongoRiverBrazzavilleMINAN,file='data/CongoRiverBrazzavilleMINAN.RData')

ArdecheRiverRaw=read.table('data-raw/ArdecheRiver.txt',header=T,sep=',')

ArdecheRiverDF=ArdecheRiverRaw[,c(1,2,3)]
colnames(ArdecheRiverDF)=c('Date','Date.h.min','stage')
ArdecheRiverDF[,1]=substr(ArdecheRiverDF$Date,start=1,stop=10)
ArdecheRiverDF[,2]=substr(ArdecheRiverDF$Date.h.min,start=1,stop=10)

ArdecheRiverDF$Date=as.Date(ArdecheRiverDF$Date)
ArdecheRiverDF$Date.h.min=as.Date(ArdecheRiverDF$Date.h.min)

data_by_year=split(ArdecheRiverDF,format(ArdecheRiverDF$Date,'%Y'))

# Subset by year
min_values <- lapply(data_by_year, function(subset) {
  min_value <- min(subset$stage)
  return(min_value)
})

# Find minimal stages by year and add uncertainty following a normal distribution
# with mean 0.08 meters and standard deviation of 0.04 meters
set.seed(1)
ArdecheRiver <- data.frame(
  Year = as.numeric(names(min_values)),
  H = unlist(min_values),
  uH = round(abs(rnorm(length(unlist(min_values)),mean = 0.08,sd=0.04)),3)
)
save(ArdecheRiver,file='data/ArdecheRiver.RData')

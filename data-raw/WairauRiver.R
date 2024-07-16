
alldata=xlsx::read.xlsx("data-raw/Gaugings_Wairau_BarnettsBank.xlsx",
                        sheetName = 'Export_Tideda')

WairauRiverGaugingsRaw=alldata[,c('time','stage.mm','flow.m3.s')]
WairauRiverGaugings = data.frame(Date=WairauRiverGaugingsRaw$time,
                                 H=WairauRiverGaugingsRaw$stage.mm/1000,
                                 Q=WairauRiverGaugingsRaw$flow.m3.s,
                                 uQ=WairauRiverGaugingsRaw$flow.m3.s*8/100)


save(WairauRiverGaugings,file='data/WairauRiverGaugings.RData')

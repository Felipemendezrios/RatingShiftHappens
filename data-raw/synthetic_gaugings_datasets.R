# File
file_path='data-raw/synthetic_gaugings_datasets.xlsx'

# get names of the sheets
wb=xlsx::loadWorkbook(file_path)
sheet_names <- names(xlsx::getSheets(wb))

# Get all data from all sheets
all_sheets <- lapply(sheet_names, function(sheet) {
  xlsx::read.xlsx(file_path, sheetName = sheet)
})

# names of shift time of each class
class_names <- stringr::str_extract(sheet_names, "Class \\d+")
class_names <- gsub(" ", "_", class_names)

class_names_f=paste0('Shift_time_',class_names)

SyntheticGaugingDataSets <- c()
j=1
# Separate synthetic data from shift time information and store in a list
for(i in 1:length(all_sheets)){

  # get synthetic data
  dataset.p <- all_sheets[[i]][, !names(all_sheets[[i]]) %in% "shift_time"]

  # get shift time
  df_shift_time.na <- all_sheets[[i]][, "shift_time", drop = FALSE]
  df_shift_time <- df_shift_time.na[!is.na(df_shift_time.na$shift_time),,drop=FALSE]

  SyntheticGaugingDataSets[[j]]=list(dataset.p,df_shift_time)
  names(SyntheticGaugingDataSets[[j]]) <- c(sheet_names[i],class_names_f[i])
  j=j+1
}

save(SyntheticGaugingDataSets,
     file='data/SyntheticGaugingDataSets.RData')

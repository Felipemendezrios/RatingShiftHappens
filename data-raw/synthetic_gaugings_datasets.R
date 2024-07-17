# File
file_path='data-raw/synthetic_gaugings_datasets.xlsx'

# get names of the sheets
wb=xlsx::loadWorkbook(file_path)
sheet_names <- names(xlsx::getSheets(wb))

# Get all data from all sheets
all_sheets <- lapply(sheet_names, function(sheet) {
  xlsx::read.xlsx(file_path, sheetName = sheet)
})

# Renames future list
class_names <- stringr::str_extract(sheet_names, "Class \\d+")
class_names <- gsub(" ", "_", class_names)

class_names_f=paste0('Shift_time_',class_names)


# Separate synthetic data from shift time information and store in a list
for(i in 1:length(all_sheets)){
  synthetic_gauging <- c()
  j=1
  # get sysnthetic data
  synthetic_gauging[[j]] <- all_sheets[[i]][, !names(all_sheets[[i]]) %in% "shift_time"]
  j=j+1
  # get shift time
  df_shift_time <- all_sheets[[i]][, "shift_time", drop = FALSE]
  synthetic_gauging[[j]] <- df_shift_time[!is.na(df_shift_time$shift_time),,drop=FALSE]

  names(synthetic_gauging) <- c(sheet_names[i],class_names_f[i])

  save(synthetic_gauging,
       file=paste0('data/synthetic_gauging_',class_names[i],'.RData'))
}

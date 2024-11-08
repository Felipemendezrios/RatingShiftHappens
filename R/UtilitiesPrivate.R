#' Check size vectors
#'
#' Function to check if the vectors have the same length
#' @param ... real vectors
#'
#' @return logical, return null if the vectors have not the same length
#' @keywords internal
#' @source \url{https://www.simonqueenborough.info/R/basic/lessons/lapply_and_sapply.html}
check_vector_lengths <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1){
    return(NULL)
  }else{
    return('ok')
  }
}

#' Transform numeric format
#'
#' Numeric date to input format of time
#'
#' @param numeric.date real vector, numeric date
#' @param origin.date value, origin date in character, POSIXct or Date format
#'
#' @return numeric date to input format of time
#' @keywords internal
NumericFormatTransform <- function(numeric.date,origin.date){
  if (lubridate::is.Date(origin.date)) {
    return(origin.date+numeric.date)
  }else if(lubridate::is.POSIXct(origin.date)){
    return(origin.date+numeric.date*86400)
  }else{ # character format
    # Attempt to parse the date using various formats
    parsed_date <- lubridate::parse_date_time(origin.date, tz = "UTC", orders  = c('ymd H:M:S', 'ymd', 'mdy',
                                                                                   'dmy', 'ymd HMS', 'y-m-d H:M:S',
                                                                                   'y/m/d H:M:S', 'y/m/d HMS' ))

    return(parsed_date+numeric.date*86400)
  }
}


#' Check square matrix
#'
#' @param x matrix, hydraulic controls
#'
#' @return logical, return null if matrix is not square
#' @keywords internal
check_square_matrix <- function(x){
  if(ncol(x)==nrow(x)){
    return('ok')
  }else{
    return(NULL)
  }
}

#' Check number of parameters depending of distribution
#'
#' @param distribution string, distribution
#' @param prior real vector, prior information about parameter
#'
#' @return logical, return null if number of parameter does not match with specified distribution
#' @keywords internal
check_param_distribution <- function(distribution, prior){
  if(distribution == 'Gaussian' && length(prior)==2){
    return('ok')
  }else if(distribution == 'LogNormal' && length(prior)==2){
    return('ok')
  }else if(distribution == 'Uniform' && length(prior)==3){
    return('ok')
  }else{
    return(NULL)
  }
}

#' Convert a list to a data frame if possible
#'
#' @param listx list
#'
#' @return data frame or list
#' @keywords internal
List_to_DF <- function(listx){
  col_nb=(unlist(lapply(listx, ncol)))
  all_equal <- all(col_nb==col_nb[1])

  if(all_equal){
    listxpost=do.call(rbind, listx)
    rownames(listxpost) <- NULL
  }else{
    listxpost=listx
  }
  return(listxpost)
}

#' Copy all files from a folder to a another folder
#'
#' @param dir.source character, source directory
#' @param dir.destination character, destination directory
#'
#' @keywords internal
copy_files_to_folder <- function(dir.source,
                                 dir.destination){
  # Ensure the destination directory exists
  if (!dir.exists(dir.destination)) {
    dir.create(dir.destination)
  }else{
    # erase old files
    unlink(file.path(dir.destination,"*"))
  }

  # List all files in the source directory recursively
  files <- list.files(dir.source, full.names = TRUE, recursive = FALSE)

  if(any(files==dir.destination)){
    files=files[-which(files==dir.destination)]
  }

  # Copy each file to the destination directory
  invisible(file.copy(files,
                      dir.destination,
                      recursive=TRUE))
}

#' Remove files excepting some specific files from the directory
#'
#' @param files_to_keep character, files to keep
#' @param dir.source character, source directory
#'
#' @keywords internal
remove_files <- function(dir.source,
                         files_to_keep=NULL){
  # Remove all files except Results_Cooking.txt and Residual folder
  all_files <- list.files(dir.source, full.names = TRUE)

  if(!is.null(files_to_keep)){
    files_to_keep_position <- which(!is.na(match(basename(all_files), files_to_keep)))

    files_to_remove <- all_files[-files_to_keep_position]

    # Remove the files
    lapply(files_to_remove, function(file) {
      file.remove(file)
    })
  }else{
    # Remove the files
    lapply(all_files, function(file) {
      file.remove(file)
    })
  }


}

#' Check if any column as a NA when convert to numeric format
#'
#' @param donnees data frame
#'
#' @return data frame, logical value. if `TRUE`, any NA was detected during transformation.
#' @keywords internal
character_check <- function(donnees){
  sapply(donnees, function(x) {
    # Attempt to convert to numeric and check for NA
    any(is.na(as.numeric(x)) & !is.na(x))
  })
}

#' Find recession models
#'
#' @return vector, return equation of all recession model available
#' @keywords internal
Recession_models_available <- function(){
  # List all functions exported from the package
  all_functions <- ls("package:RatingShiftHappens")

  # Filter functions with names containing both "Recession" and "Equation" (Be careful to use only this combination for recession models!!)
  recession_equation_functions <- all_functions[grep("Recession.*Equation", all_functions)]

  # Initialize an empty vector to store the return values as characters
  results <- character()

  # Execute each function and store its return as a character in the results vector
  for (func_name in recession_equation_functions) {
    # Retrieve the function using `get`
    func <- get(func_name, envir = asNamespace("RatingShiftHappens"))

    # Call the function and convert its return value to a character
    result <- as.character(func())

    # Append the result to the results vector
    results <- c(results, result)
  }
  return(results)
}

#' Identify model selected is in the catalog
#'
#' @param current_model name function, model selected
#' @param model_list list, models available in the file `RecessionCurveEquationsModels.R`
#'
#' @return equation of the model selected
#' @keywords internal
identify_model <- function(current_model, model_list) {
  # Check if exact or partial matches are present
  match_idx <- grep(current_model, model_list, fixed = TRUE)

  # If no exact match is found, try partial matching
  if (length(match_idx) == 0) {
    match_idx <- grep(substr(current_model, 1, 10), model_list)  # Adjust substring length as needed
  }

  # Return the matched model index (if any)
  if (length(match_idx) > 0) {
    return(model_list[match_idx[1]])
  } else {
    return(NA)  # No match found
  }
}
#' Convert list to data frame
#'
#' @param liste_df list, data frame
#'
#' @return data frame if possible, if not a list
#'
#' @keywords internal
#' @examples
#'
#' sample_data <- list(
#' A = data.frame(a=1, b=-2, c=0, d=4),
#' B = data.frame(a=-3, b=5, c=0, d=-1))
#'
#' # Convert to data frame
#' convert_list_to_dataframe(sample_data)
#'
#' sample_data <- list(
#' A = data.frame(a=1, c=0, d=4),
#' B = data.frame(a=-3, b=5, c=0, d=-1))
#'
#' # Convert to data frame
#' convert_list_to_dataframe(sample_data)
convert_list_to_dataframe <- function(liste_df){
  result=c()
  lengths <- sapply(liste_df, ncol)
  if (length(unique(lengths)) > 1) {
    warning("The number of the components in the list is not the same, replace by NA for missing values")
    maxncol <- max(lengths)
    list_maxncol=which(lengths==maxncol)[1]
    namescol=names(liste_df[[list_maxncol]])

    for (i in seq_along(liste_df)) {
      if(length(liste_df[[i]])==maxncol){
        result <- rbind(result, liste_df[[i]])
      }else{
        values=liste_df[[i]]
        values_with_NA = data.frame(rbind(rep(NA,maxncol)))
        colnames(values_with_NA)=namescol

        values_with_NA[,c(match(names(values),namescol))]=values

        result <- rbind(result,values_with_NA)
      }

    }
    return(result)
  }

  for (i in seq_along(liste_df)) {
    result <- rbind(result, liste_df[[i]])
  }
  rownames(result)<-NULL
  return(result)
}

#' Replace negative and/or zero values to NA or a specified value
#'
#' @param data_frame data frame, data to be checked and replaced
#' @param columns string, exact names of the columns from data frame
#' @param consider_zero logical, if `TRUE` zero values will be replace by NA. Otherwise, zero values are not be replaced.
#' @param replace NA or real value, data will be replace when condition are not accepted
#'
#' @return data frame transformed with NA values
#' @keywords internal
#'
#' @examples
#'
#' sample_data <- data.frame(
#' A = c(1, -2, 0, 4),
#' B = c(-3, 5, 0, -1),
#' C = c(0, 2, -7, 6) )
#'
#' # Display the original data frame
#' print(sample_data)
#'
#' # Function call with consider_zero = TRUE
#' replace_negatives_or_zero_values(sample_data, c("A", "B", "C"), consider_zero = TRUE)
#'
#' # Function call with consider_zero = FALSE
#' replace_negatives_or_zero_values(sample_data, c("A", "B", "C"), consider_zero = FALSE)
replace_negatives_or_zero_values <- function(data_frame, columns='all', consider_zero = TRUE, replace=NA) {

  if(all(columns=='all')){
    columns.match = seq(1,ncol(data_frame))
  }else{
    columns.match = match(columns,colnames(data_frame))
  }
  columns.df <- columns.match[!is.na(columns.match)]

  if(length(columns.df)==0)stop('Verify columns input. Names of the columns do not match with the columns names of the data frame')

  if(!is.na(replace)&&!is.numeric(replace))stop('replace must be NA or a numeric value')
  result <- apply(data_frame[, columns.df, drop = FALSE], 2, function(x) {
    if (consider_zero) {
      ifelse(x <= 0, replace, x)
    } else {
      ifelse(x < 0, replace, x)
    }
  })
  return(data.frame(result))
}

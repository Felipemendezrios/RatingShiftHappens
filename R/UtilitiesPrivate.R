#' Check size vectors
#'
#' Function to check if the vectors have the same length
#' @param ... real vectors
#'
#' @return logical, return null if the vectors have not the same length
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
NumericFormatTransform <- function(numeric.date,origin.date){
  date_time <- as.POSIXct(origin.date)
  date_only <- as.Date(date_time)
  contains_time_info <- date_time != as.POSIXct(date_only)
  if (any(contains_time_info)) {
    return(origin.date+numeric.date*86400)
  }else{
    return(origin.date+numeric.date)
  }
}


#' Check square matrix
#'
#' @param x matrix, hydraulic controls
#'
#' @return logical, return null if matrix is not square
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

#' Convert à list to a data frame if possible
#'
#' @param listx list
#'
#' @return data frame or list
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
character_check <- function(donnees){
  sapply(donnees, function(x) {
    # Attempt to convert to numeric and check for NA
    any(is.na(as.numeric(x)) & !is.na(x))
  })
}

#' Find recession models
#'
#' @return vector, return equation of all recession model available
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


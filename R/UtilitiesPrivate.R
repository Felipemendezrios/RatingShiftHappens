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
  return(origin.date+numeric.date*86400)
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
    dir.create(dir.destination, recursive = TRUE)
  }else{
    # erase old files
    unlink(file.path(dir.destination,"*"))
  }

  # List all files in the source directory recursively
  files <- list.files(dir.source, full.names = TRUE, recursive = TRUE)

  # Copy each file to the destination directory
  lapply(files, function(file) {
    # Construct the destination file path
    destination_file <- file.path(dir.destination)

    # Copy the file
    file.copy(file, destination_file, overwrite = TRUE)
  })
}

#' Remove files excepting some specific files from the directory
#'
#' @param files_to_keep character, files to keep
#' @param dir.source character, source directory
#'
remove_files <- function(dir.source,
                         files_to_keep){
  # Remove all files except Results_Cooking.txt and Residual folder
  all_files <- list.files(dir.source, full.names = TRUE)

  files_to_keep_position <- which(!is.na(match(basename(all_files), files_to_keep)))

  files_to_remove <- all_files[-files_to_keep_position]

  # Remove the files
  lapply(files_to_remove, function(file) {
    file.remove(file)
  })

}


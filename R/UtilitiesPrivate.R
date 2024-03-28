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
#' @param class string value, class of a vector (character, POSIXct, Date)
#' @param origin.date value, origin date in character, POSIXct or Date format
#'
#' @return numeric date to input format of time
NumericFormatTransform <- function(numeric.date,class,origin.date){
  original.format.date=switch(class,
                              'Date'=origin.date+numeric.date,
                              'POSIXct' = origin.date+numeric.date*86400,
                              'character' = origin.date+numeric.date*86400
  )
  return(original.format.date)
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

convert_list_to_dataframe <- function(liste_df){

  lengths <- sapply(liste_df, ncol)
  if (length(unique(lengths)) > 1) {
    warning("Parameters have not the same size for all estimation of the rating curves")
    return(liste_df)
  }

  result=c()
  for (i in seq_along(liste_df)) {
    result <- rbind(result, liste_df[[i]])
  }
  return(result)
}

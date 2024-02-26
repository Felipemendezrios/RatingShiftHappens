#' Check size vectors
#'
#' Function to check if the vectors have the same length
#' @param ... real vectors
#'
#' @return logical, return null if the vectors have not the same length
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

#' Check size vectors
#'
#' Function to check if the vectors have the same length
#' @param ... real vectors
#'
#' @return logical, return null if the vectors have not the same length
#' @export
#'
#' @examples
#'
#' check_vector_lengths(RhoneRiver$Time,RhoneRiver$H,RhoneRiver$Q,RhoneRiver$uQ)
#' check_vector_lengths(RhoneRiver$Time,RhoneRiver$H,RhoneRiver$Q,RhoneRiver$uQ[-1])
check_vector_lengths <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1){
    return(NULL)
  }else{
    return('ok')
  }
}

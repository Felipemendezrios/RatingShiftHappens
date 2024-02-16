# Function to check if the vectors have the same length

check_vector_lengths <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1){
    return(NULL)
  }else{
    return('ok')
  }
}

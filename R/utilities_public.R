#' Date to numeric format
#'
#' Get date format and transform date to numeric
#'
#' @param date vector, date in POSIXct, character or Date format
#'
#' @return List with the following components :
#' \enumerate{
#'   \item origin: value, denoting the oldest date from the input variable
#'   \item time: real vector, representing time in numeric format. This denotes the number of days from origin
#' }
#' @import lubridate
#' @export
#'
#' @examples
#' TransformDateFormat(
#'                    c(as.POSIXct('2024-04-19 12:30:00',
#'                      tz='UTC'),
#'                      as.POSIXct('2024-03-19 18:30:00',
#'                      tz='UTC')))
#'
#' TransformDateFormat(c(as.Date('2024/02/19'),
#'                       as.Date('2024/02/10')))
#'
#' TransformDateFormat(c(as.character('20240419'),
#'                       as.character('20240119')))
TransformDateFormat <- function(date){

  # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
  if(lubridate::is.Date(date)){
    origin=min(date)
    # Interval in days between origin and each date
    diff_days <- lubridate::time_length(lubridate::interval(origin, date), "day")
   return(list(origin=origin,
               time=diff_days))

  }else if(!lubridate::is.POSIXct(date)){

    date_string <- as.character(date)

    # Attempt to parse the date using various formats
    parsed_date <- lubridate::parse_date_time(date_string, tz = "UTC", orders  = c('ymd H:M:S', 'ymd', 'mdy',
                                                                                   'dmy', 'ymd HMS', 'y-m-d H:M:S',
                                                                                   'y/m/d H:M:S', 'y/m/d HMS' ))

    if(any(is.na(parsed_date)))stop('The format is not supported; please verify the input date or time format')

    # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
    origin=min(parsed_date)
    diff_days <- lubridate::time_length(lubridate::interval(origin, parsed_date), "day")
    return(list(origin=origin,
                time=diff_days))

  }else{
    # Read time zone
    date.transf <- date
  }
  origin=min(date.transf)
  diff_days <- lubridate::time_length(lubridate::interval(origin, date.transf), "day")

  return(list(origin=origin,
              time=diff_days))
}

#' Builder function to generate prior information on the parameter as a object
#'
#' @return object, prior information on the parameter defined
#' @export
#' @importFrom RBaM getCatalogue
Builder_Prior_Knowledge <- function() {

  name <- as.character(readline(prompt = 'Enter a name of the parameter: '))

  init <- as.numeric(readline(prompt = paste0('Enter initial guess of the parameter ',name,': ')))

  # Display Available distributions
  cat('Available distributions:\n ')
  print(RBaM::getCatalogue()$distributions)

  while(TRUE){
    prior.dist = as.character(readline(prompt = 'Enter the prior distribution:'))
    if(length(which(prior.dist==RBaM::getCatalogue()$distribution))==0){
      cat('Invalid input. Please enter a prior distribution available. Please respect capital letters')
      next
    }
    break
  }

  # Parameter of the distribution
  if(prior.dist=='FIX'){
    prior.par=NULL
  }else{
    while (TRUE) {
       if (prior.dist %in% c('Gaussian', 'Uniform', 'LogNormal')) {

        params_input <- readline(prompt = paste0("Enter the parameters for ", prior.dist, " distribution (separated by a space): "))
        params <- strsplit(params_input, " ")[[1]]

        # Check if two values were entered
        if (length(params) != 2) {
          cat("Invalid input. Please enter only two parameters separated by a space.\n")
          next
        }

        # Convert the inputs to numeric
        param1 <- as.numeric(params[1])
        param2 <- as.numeric(params[2])

        # Check if both inputs are valid numbers
        if (is.na(param1) || is.na(param2)) {
          cat("Invalid input. Please enter valid numbers for both parameters.\n")
          next
        }

        # If both inputs are valid, return the parameters as a vector
        prior.par = c(param1,param2)

      } else {
        params_input <- readline(prompt = "Enter the parameters of the distribution.\n If the number of the parameters is larger than one, separated by a space): \n")
        params <- strsplit(params_input, " ")[[1]]

        # Convert the inputs to numeric
        numeric_params <- as.numeric(params)

        # Check if any input is not a valid number
        if (any(is.na(numeric_params))) {
          cat("Invalid input. Please enter valid numbers for all parameters.\n")
          next  # Skip to the next iteration of the loop
        }
        prior.par=as.vector(numeric_params)
      }
      break
    }
  }

  return(RBaM::parameter(name=name,
                         init=init,
                         prior.dist = prior.dist,
                         prior.par = prior.par))
}


#' Interactive builder for the hydraulic control matrix
#'
#' @param ncontrols integer value, number of hydraulic controls
#'
#' @return matrix, hydraulic control. This matrix consists of 0s and 1s
#' @export
Builder_ControlMatrix <- function(ncontrols) {

  control_matrix_build <- c()
  for(i in 1:ncontrols){
    while (TRUE) {
      cat('
      For each stage range, specify active control(s) using:
        1 for active control
        0 for inactive control
        The numbers must be separated by a space and must be the same size as the number of controls')
      control_i_input <- c(readline(prompt = paste0('Stage range number ',i,' (from lowest to highest) \n')))
      control_i <- as.numeric(strsplit(control_i_input, " ")[[1]])

      if(ncontrols!=length(control_i)){
        cat('The specified control number must be identical to the number of the input data entered to respect a square matrix.\nPlease verify the separator used, it must be a space')
        next
      }

      if(any(control_i!=0 & control_i!=1)){
        cat('Hydraulic control must be filled by 1 (active) and 0 (inactive) for describing hydraulic controls')
        next
      }

      control_matrix_build = rbind(control_matrix_build,control_i)
      break
    }
  }
  rownames(control_matrix_build) <- NULL
  return(control_matrix_build)
}

#' Root Mean Squared Error weighted (RMSE weighted)
#'
#' @param observations vector, observed data
#' @param predicted vector, predicted data
#' @param weights vector, weights to each data to give more or less weights to some data
#'
#' @return real value, RMSE weighted value to qualify the fit
#' @export
RMSE_Weighted <- function(observations, predicted, weights) {
  sqrt(sum(weights * (observations - predicted)^2) / sum(weights))
}

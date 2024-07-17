#' Date to numeric format
#'
#' Get date format and transform date to numeric
#'
#' @param date vector, date in date, POSIXct or character format
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
#' DateFormatTransform(
#'                    c(as.POSIXct('2024-04-19 12:30:00',
#'                      tz='UTC'),
#'                      as.POSIXct('2024-03-19 18:30:00',
#'                      tz='UTC')))
#'
#' DateFormatTransform(c(as.Date('2024/02/19'),
#'                       as.Date('2024/02/10')))
#'
#' DateFormatTransform(c(as.character('20240419'),
#'                       as.character('20240119')))
DateFormatTransform <- function(date){

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

#' Replace negative and/or zero values to NA or a specified value
#'
#' @param data_frame data frame, data to be checked and replaced
#' @param columns string, exact names of the columns from data frame
#' @param consider_zero logical, if `TRUE` zero values will be replace by NA. Otherwise, zero values are not be replaced.
#' @param replace NA or real value, data will be replace when condition are not accepted
#'
#' @return data frame transformed with NA values
#' @export
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

#' Convert list to data frame
#'
#' @param liste_df list, data frame
#'
#' @return data frame if possible, if not a list
#' @export
#'
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

#' Builder function to generate prior information on the parameter as a object
#'
#' @return object, prior information on the parameter defined
#' @export
#' @importFrom RBaM getCatalogue
prior_infor_param_builder <- function() {

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


#' Builder hydraulic control matrix
#'
#' @param ncontrols integer value, number of hydraulic controls
#'
#' @return matrix, hydraulic control
#' @export
control_matrix_builder <- function(ncontrols) {

  control_matrix_build <- c()
  for(i in 1:ncontrols){
    while (TRUE) {
      cat('
      Describe the control-by-control matrix using :
        1 to active control
        0 to inactive control
        The number must be separed by a space and must be the same size as the number of controls')
      control_i_input <- c(readline(prompt = paste0('Hydraulic control N° ',i,' \n')))
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


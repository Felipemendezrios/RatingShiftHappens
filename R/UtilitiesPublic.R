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
#'
#' @export
#'
#' @examples
#' DateFormatTransform(c(as.POSIXct('2024-04-19 12:30:00',tz='UTC'),as.POSIXct('2024-03-19 18:30:00',tz='UTC')))
#'
#' DateFormatTransform(c(as.Date('2024/02/19'),as.Date('2024/02/10')))
#'
#' DateFormatTransform(c(as.character('20240419'),as.character('20240119')))
DateFormatTransform <- function(date){

  # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
  if(lubridate::is.Date(date)){

   return(list(origin=min(date),
               time=as.numeric(difftime(date, min(date), units = "days"))))

  }else if(!lubridate::is.POSIXct(date)){

    date_string <- as.character(date)

    # Attempt to parse the date using various formats
    parsed_date <- lubridate::parse_date_time(date_string, tz = "UTC", orders  = c('ymd H:M:S', 'ymd', 'mdy',
                                                                                   'dmy', 'ymd HMS', 'y-m-d H:M:S',
                                                                                   'y/m/d H:M:S', 'y/m/d HMS' ))

    if(any(is.na(parsed_date)))stop('The format is not supported; please verify the input date or time format')

    # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
    return(list(origin=min(parsed_date),
                time=as.numeric(difftime(parsed_date, min(parsed_date), units = "days"))))

  }else if(lubridate::tz(date)=='CET'){
    # Transform to UTC using force_tz
    warning('Time expressed in UTC')
    date_UTC <- lubridate::force_tz(date, tzone = "UTC")
  }else if(lubridate::tz(date)=="" || is.null(lubridate::tz(date))){
    warning('Timezone not explicityly set or not recognized. Assumption : UTC')
    date_UTC <- as.POSIXct(date, origin = '1970-01-01',tz='UTC')
  }else if(lubridate::tz(date)=='UTC'){
    date_UTC <- as.POSIXct(date, origin = '1970-01-01',tz='UTC')
  }

  return(list(origin=min(date_UTC),
              time=as.numeric(difftime(date_UTC, min(date_UTC), units = "days"))))
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
replace_negatives_or_zero_values <- function(data_frame, columns, consider_zero = TRUE, replace=NA) {
  if(!is.na(replace)&&!is.numeric(replace))stop('replace must be NA or a numeric value')
  result <- apply(data_frame[, columns, drop = FALSE], 2, function(x) {
    if (consider_zero) {
      ifelse(x <= 0, replace, x)
    } else {
      ifelse(x < 0, replace, x)
    }
  })
  return(result)
}


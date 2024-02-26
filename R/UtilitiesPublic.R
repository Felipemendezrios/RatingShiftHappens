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



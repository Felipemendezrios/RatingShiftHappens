#' Date to time format
#'
#' Get date format and transform date to numeric
#'
#' @param date vector, date
#'
#' @return List with the following components :
#' \enumerate{
#'   \item date: vector, date in POSIXct format
#'   \item time: real vector, time in numeric format where origin default is "1970-01-01 00:00:00 UTC"
#'   }
#' @export
#'
#' @examples
#' datetime=DateFormatTransform('2024-02-19 12:30:00')
#' datetime
#' as.POSIXct(datetime, origin = '1970-01-01',tz='UTC')
#'
#' datetime2=DateFormatTransform('2024/02/19')
#' datetime2
#' as.POSIXct(datetime2, origin = '1970-01-01',tz='UTC')
DateFormatTransform <- function(date){
  # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
  if(lubridate::is.Date(date)){

    # Attempt to parse the date using various formats
    parsed_date <- as.POSIXct(date, tz = "UTC")

    if(any(is.na(parsed_date)))stop('The format is not supported; please verify the input date or time format')

    return(as.numeric(parsed_date))

  }else if(!lubridate::is.POSIXct(date)){

    date_string <- as.character(date)

    # Attempt to parse the date using various formats
    parsed_date <- as.POSIXct(date_string, tz = "UTC", format = c('%Y-%m-%d %H:%M:%S',
                                             '%Y/%m/%d %H:%M:%S',
                                             '%Y%m%d %H:%M:%S',
                                             '%m/%d/%Y %H:%M:%S',
                                             '%Y%m%d%H%M%S',
                                             '%Y%m%dT%H%M%S',
                                             '%Y-%m-%d %H:%M',
                                             '%Y/%m/%d %H:%M',
                                             '%Y-%m-%d',
                                             '%Y/%m/%d',
                                             '%Y%m%d',
                                             '%m/%d/%Y'))

    parsed_date <- parsed_date[!is.na(parsed_date)][1]  # Select the first non-NA result

    if(any(is.na(parsed_date)))stop('The format is not supported; please verify the input date or time format')

    # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
    return(as.numeric(parsed_date))
  }else if(lubridate::tz(date)=='CET'){
    # Transform to UTC using force_tz
    warning('Time expressed in UTC')
    date <- lubridate::force_tz(date, tzone = "UTC")
  }else if(lubridate::tz(date)=="" || is.null(lubridate::tz(date))){
    warning('Timezone not explicityly set or not recognized. Assumption : UTC')
    date <- as.POSIXct(date, origin = '1970-01-01',tz='UTC')
  }
  return(as.numeric(date))
}



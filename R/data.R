#' Floods of the Rhone River at Beaucaire, France
#'
#' A data frame containing annual maximum stages (H, in meters) et discharges
#' (Q, in cubic meters per second) for the Rhone River at Beaucaire, France,
#' along with the associated uncertainties expressed as standard deviations (uH and uQ).
#' Details on the reconstruction of of these long series can be found in the article
#' by Lucas et al. (2023) referenced below.
#' Note that years 1968, 1969 and 1970 are missing and are not included in the data frame.
#'
#' @source \url{https://doi.org/10.1016/j.jhydrol.2023.129840}
"RhoneRiver"

#' Stage record of the Ardèche River at Meyras, France, provided by UHPC Grand Delta
#'
#' A data frame containing stages (H, in meters) for the Ardèche River
#' at Meyras, France, from 07/11/2001 until 29/10/2018
#'
#' @format
#' \describe{
#'   \item{Date}{Date time}
#'   \item{H}{Stage record}
#' }
#' More information about dataset available in :
#' @source\url{https://www.theses.fr/2021GRALU006}
"ArdecheRiverMeyrasStage"

#' Gauging of the Ardèche River at Meyras, France, provided by UHPC Grand Delta
#'
#' A data frame containing stages (H, in meters) and discharges ADCP measurements
#' (Q, in cubic meters per second) for the Ardèche River at Meyras, France,
#' along with the associated uncertainties expressed as standard deviations (uQ) from 2001 until 2018
#'
#' @format
#' \describe{
#'   \item{Day}{Day}
#'   \item{Month}{Month}
#'   \item{Year}{Year}
#'   \item{Hour}{Hour}
#'   \item{Minute}{Minute}
#'   \item{Second}{Second}
#'   \item{Date}{Date time}
#'   \item{H}{Stage record}
#'   \item{Q}{Discharge ADCP measurement}
#'   \item{uQ}{Uncertainty expressed as standard deviation}
#' }
#' More information about dataset available in :
#' @source\url{https://www.theses.fr/2021GRALU006}
"ArdecheRiverMeyrasGaugings"

#' MINAN of the discharge simulation obtained by BaRatinAGE for the Congo River at Brazzaville, Republic of the Congo.
#'
#' A data frame containing discharges simulations (Q, in cubic meters per second)
#' for the Congo River at Brazzaville, Republic of the Congo, along with
#' the associated uncertainties expressed as standard deviations (uQ) from 1902 until 2024
#'
#' @format
#' \describe{
#'   \item{time}{Date time}
#'   \item{year}{Year}
#'   \item{Q}{MINAN of the Discharge simulation }
#'   \item{uQ}{Uncertainty expressed as standard deviation}
#' }
"CongoRiverBrazzavilleMINAN"

#' Gauging of the Wairau River at Barnetts Bank, Tuamarina, New Zealand, Tideda software export.
#'
#' A data frame containing stages (H, in meters) and discharges ADCP measurements
#' (Q, in cubic meters per second) for the Wairau River at Barnetts Bank, Tuamarina, New Zealand,
#' from 1999 until 2015.
#' Any information about uncertainty in discharge measurement; hence a uncertainty of 8%
#' has been considered and it was expressed as standard deviations (uQ).
#'
#' @format
#' \describe{
#'   \item{Date}{Date time}
#'   \item{H}{Stage record}
#'   \item{Q}{Discharge ADCP measurement}
#'   \item{uQ}{Uncertainty expressed as standard deviation}
#' }
"WairauRiverGaugings"

#' Synthetic dataset (Darienzo, et al. 2021)
#'
#' A list containing sub-lists with six classes of synthetic datasets.
#' Each sub-list is composed of two parts. The first part is a data frame containing the stage (H, in meters)
#' and discharges ADCP measurements (Q, in cubic meters per second) for a hypothetical case, as used
#' by Darienzo et al. 2021, along with the associated uncertainties expressed as standard deviations (uQ).
#' In addition, the period is also defined, indicating the same gauging family.
#' The second part of the list contains the real shift time of each class.
#'
#' @format
#' All list are the classes. All classes have the same columns; hence, Only the two first data frames will be explained.
#' \describe{
#'   \item{\code{Class 1:}}{Few gaugings
#'      \itemize{
#'        \item time: Date time
#'        \item period: Period same family of gauging
#'        \item h: Stage record
#'        \item Q: Discharge ADCP measurement
#'        \item uQ: Uncertainty expressed as standard deviation
#'       }
#'   }
#'   \item{\code{Class 1 :}}{Shift time
#'      \itemize{
#'        \item shift_time: shift time detected in this class
#'       }
#'   }
#'   \item{\code{Class 3:}}{Many gaugins}
#'   \item{\code{Class 4:}}{Small shifts}
#'   \item{\code{Class 6:}}{Many shifts}
#'   \item{\code{Class 9:}}{Very uncertain gaugings}
#'   \item{\code{Class 10:}}{Three controls}
#' }
#' @source \url{https://doi.org/10.1029/2020WR028607}
"synthetic_gauging_datasets"


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
#' (Q, in cubic meters per second) for the Ardèche River at Meyras, France, along with
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


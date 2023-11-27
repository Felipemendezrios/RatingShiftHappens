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

#' Annual minimum stage record of the Ardèche River at Meyras, France
#'
#' A data frame containing annual minimum stages (H, in meters) for the Ardèche River
#' at Meyras, France, along with by associated random uncertainties from 1988 until 2022.
#' Assumption: uncertainties follow a normal distribution with a mean of 0.08 meters
#' and a standard deviation of 0.04 meters.
#'
#' @format
#' \describe{
#'   \item{Year}{Year}
#'   \item{H}{Stage record}
#'   \item{uH}{Uncertainty expressed as standard deviation (Random realisation)}
#' }
#' Details from stage record available in :
#' @source\url{https://hydro.eaufrance.fr/stationhydro/V500403001/series}
"ArdecheRiver"



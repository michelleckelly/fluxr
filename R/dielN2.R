#' Example diel dataset from 2024 stream mesocosm experiments
#'
#' This is a subset of data collected during a 2024 stream mesocosm experiment
#' with collaborators Jason Taylor and Amy Marcarelli at the University of
#' Mississippi Field Station. 
#' 
#' We collected dissolved gas samples at upstream and
#' downstream sampling locations (three replicates per location) every two hours
#' for a period of 30 hours. We then analyzed dissolved gas samples using a MIMS,
#' and converted raw signal data to N2, O2, and Ar concentrations using the
#' mimsy R package. 
#' 
#' We fit a GAM model to dissolved gas concentrations (2-hr
#' interval) to create a smoothed time series (15-min interval). light.calc was
#' calculated using the `calc_light` function from the `streamMetabolizer` R
#' package.
#'
#' @format ## `dielN2`
#' A data frame with 224 rows and 17 columns:
#' \describe{
#'    \item{solar.time}{Solar time}
#'    \item{location}{Sampling location}
#'    \item{N2.obs}{Dinitrogen gas concentration, mg L^{-1}}
#'    \item{N2.sat}{Dinitrogen gas concentration at saturation, mg L^{-1}}
#'    \item{DO.obs}{Dissolved oxygen concentration, mg L^{-1}}
#'    \item{DO.sat}{Dissolved oxygen concentration at saturation, mg L^{-1}}
#'    \item{Ar.obs}{Argon gas concentration, $mg L^{-1}$}
#'    \item{Ar.sat}{Argon gas concentration at saturation, mg L^{-1}}
#'    \item{temp.water}{Water temperature, deg C}
#'    \item{light.data}{PAR measured by light sensor next to sampling location, μmol m^{-2} s^{-1}}
#'    \item{light.smooth}{Measured PAR smoothed by GAM function, μmol m^{-2} s^{-1}}
#'    \item{light.calc}{PAR calculated with `streamMetabolizer::calc_light()`, μmol m^{-2} s^{-1}}
#'    \item{depth}{Stream depth, m}
#'    \item{discharge}{Stream discharge, m^3 s^{-1}}
#'    \item{tt}{Travel time between sampling locations, min}
#'    \item{pressure}{Barometric pressure, mbar}
#' }
"dielN2"

#' \code{Kcor_N2} Estimate kN2, given water temperature and K600
#'
#' Esimate k values (not temperature or gas normalized) from K600 (temperature and gas normalized). See Wanninkhof 1992 and Hall et al 2016 supplemental material.
#'
#' @param temp Water temperature in degrees C
#' @param K600 Normalized K600 reaeration rate
#' @param n Default is 0.5 for wavy water, but 0.667 refers to still water, see Jahne 1987
#'
#' @returns Numeric vector
#'
#' @references Raymond, P. A., Zappa, C. J., Butman, D., Bott, T. L., Potter, J., Mulholland, P., et al. (2012). Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers: Gas transfer velocity and hydraulic geometry. Limnology and Oceanography: Fluids and Environments, 2(1), 41–53. https://doi.org/10.1215/21573689-1597669
#'
#' @export
#'
Kcor_N2 <- function(temp, K600, n = 0.5) {
  # Values are from metaanalysis in Raymond et al 2012, Table 1
  A <- 1615
  B <- -92.15
  C <- 2.349
  D <- -0.0240

  # Calculate Schmidt number for gas
  sc_gas <- A + B*temp + C*(temp^2) + D*(temp^3)

  # Convert K600 values into k_gas values
  k_gas <- K600 / ((600/sc_gas)^(-n))
  return(k_gas)
}

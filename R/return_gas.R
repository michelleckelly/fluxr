#' \code{return_gas} Return a data.frame containing solar time, observed gas, and modeled gas
#' 
#' @keywords internal
return_gas <- function(modname, instresults, inputdata){

  if(modname == "o2_twostation"){
    gasout <-
      data.frame(
        solar.time = instresults$solar.time,
        DO.obs.down = 
          inputdata$DO_obs_down,
        DO.mod.down = 
          instresults$DO_mod_down_mean,
        DO.mod.down.lower = 
          instresults$DO_mod_down_mean - 
          instresults$DO_mod_down_sd,
        DO.mod.down.upper = 
          instresults$DO_mod_down_mean + 
          instresults$DO_mod_down_sd
      )
  }
  
  if(modname == "n2_twostation_nifong"){
    gasout <- 
      data.frame(
        solar.time = instresults$solar.time,
        N2.obs.down = 
          inputdata$N2_obs_down,
        N2.mod.down = 
          instresults$N2_mod_down_mean,
        N2.mod.down.lower = 
          instresults$N2_mod_down_mean - 
          instresults$N2_mod_down_sd,
        N2.mod.down.upper = 
          instresults$N2_mod_down_mean + 
          instresults$N2_mod_down_sd
      )
  }
  
  return(gasout)
}
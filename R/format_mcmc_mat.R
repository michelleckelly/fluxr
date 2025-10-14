#' Format MCMC output into a one-row data.frame
#'
#' For split_dates models. Formats output into a one-row data.frame for
#' row-binding with other such data.frames
#'
#' @param mcmc_mat matrix as extracted from Stan
#' @param names_params character vector of the names of the parameters
#' @param names_stats character vector of the names of the statistics
#' @import dplyr
#' @keywords internal
format_mcmc_mat <- function(mcmc_mat, #names_params, names_stats, 
                            keep_mcmc, runmcmc_out, data_list) {
  # Initialize a list to fill with modeled parameters, as in 
  # streamMetabolizer::format_mcmc_mat_nosplit()
  par_homes <-
    list(
      overall = c(
        "err_obs_iid_sigma", "err_obs_iid_sigma_scaled", "lp__"
      ),
      daily = c(
        "GPP", "ER", "GPP_daily",
        #"Pmax", "alpha",
        "ER_daily", "K600_daily", "DO_R2"
      ),
      inst = c(
        "DO_mod_down", "GPP_inst", "ER_inst", "KO2_inst", "err_obs_iid"
      )
    )
  
  # Get number of rows of each variable
  var_table <- table(gsub("\\[[[:digit:]|,]+\\]", "", rownames(mcmc_mat)))
  var_names <- names(var_table)
  
  # Determine which variables we have from par_homes,
  # Then set the par_homes list segment to just include the variables that we 
  # have
  par_homes$overall <- 
    par_homes$overall[which(par_homes$overall %in% var_names)]
  par_homes$daily <- par_homes$daily[which(par_homes$daily %in% var_names)]
  par_homes$inst <- par_homes$inst[which(par_homes$inst %in% var_names)]
  
  # Replace percentage symbol
  colnames(mcmc_mat) <- gsub("%", "pct", colnames(mcmc_mat))
  
  # Select the data from mcmc_mat that corresponds to par_homes names
  #data.frame(par_homes$overall)
  #overall_data <- 
  #matrix(data = overall_data, 
   #      ncol = length(par_homes$overall),
    #     dimnames = list(c(NA), c(par_homes$overall)))
  
  # Overall estimates
  par_overall <- mcmc_mat[rownames(mcmc_mat) %in% par_homes$overall,]
  par_overall <- pivotfunc(par_overall, varList = rownames(par_overall))
  
  # Daily estimates
  par_daily <- mcmc_mat[rownames(mcmc_mat) %in% par_homes$daily,]
  par_daily <- pivotfunc(par_daily, varList = rownames(par_daily))
  
  # Instantaneous estimates ----------------------------------------------------
  # Format dataframe
  par_inst <- as.data.frame(mcmc_mat[str_extract(rownames(mcmc_mat), 
                                   pattern = "(\\w+)\\[(\\d+)\\]") %in% 
                         rownames(mcmc_mat),])
  # Add indexing for date times
  par_inst$datetime_index <- 
    as.numeric(str_extract(rownames(par_inst), 
                           pattern = "(\\w+)\\[(\\d+)\\]", group = 2))
  par_inst <- pivotfunc(par_inst, 
                        varList = str_extract(rownames(par_inst), 
                                              pattern = "(\\w+)\\[(\\d+)\\]", 
                                              group = 1))
  # Connect real dates to indexing numbers
  par_inst$solar.time <- NA
  # Set first entry to start time
  par_inst$solar.time[1] <- data_list$startTime
  # Loop through time index to calc solar time at each timestep
  for (i in 2:nrow(par_inst)){
    par_inst$solar.time[i] <- 
      par_inst$solar.time[i-1] + 
      data_list$timestep * 60
  } # EDITING NOTE: left off here, change below to reference par_inst 
  # rather than stan_out$inst
  
  # Convert from UNIX to UTC
  stan_out$inst$solar.time <- as.POSIXct(stan_out$inst$solar.time, 
                                         origin = "1970-01-01",
                                         tz = "UTC")
  stan_out$inst <- relocate(stan_out$inst, solar.time)
  
  mcmc_out <- 
    list(
      overall = par_overall,
      daily = par_daily,
      inst = par_inst
    )
  
  return(mcmc_out)
}

# Internal function for pivoting and widening model output df
pivotfunc <- function(data, varList){
  data <- as.data.frame(data)
  data$vars <- varList
  data <- pivot_longer(data, cols = mean:Rhat)
  data <- pivot_wider(data, names_from = c(vars, name), values_from = value)
  return(data)
}

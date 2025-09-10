#' \code{prepdata_bayes_twostation} Prepare two-station data for Bayesian modeling
#'
#' @param data Dataset containing columns: `solar.time`, `location`, `DO.obs`, `DO.sat`, `temp.water`, `light`, `depth`, `tt`, `pressure`
#' @param ply_date Date to model, YYYY-MM-DD
#' @param specs List returned by `set_specs()`
#' @param up.name Character string denoting the name of the upstream sampling station (must match name in `location` column)
#' @param down.name Character string denoting the name of the downstream sampling station (must match name in `location` column)
#'
#' @returns Numeric vector
#'
#' @references Raymond, P. A., Zappa, C. J., Butman, D., Bott, T. L., Potter, J., Mulholland, P., et al. (2012). Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers: Gas transfer velocity and hydraulic geometry. Limnology and Oceanography: Fluids and Environments, 2(1), 41–53. https://doi.org/10.1215/21573689-1597669
#'
#' @export
prepdata_bayes_twostation <- function(data, ply_date, specs, 
                                      up.name = "upstream", 
                                      down.name = "downstream"){
  # Add date column if not there already
  data$date <- lubridate::date(data$solar.time) 
  # Glue location names together
  location.names <- c(up.name, down.name)
  
  # NOTE: fixing number of days here for now, but in future this should be 
  # modified to reflect the number of days of data passed to the prep data 
  # function
  d <- 1
  
  # Get timesteps --------------------------------------------------------------
  obs_times <- unique(data$solar.time)
  # Grab time difference between measurements, min
  timestep_min <- round(as.numeric(obs_times[2] - obs_times[1]))
  # Calc n24, number of measurements in full 24 hours
  n24 <- round(1 / (timestep_min / 1440))
  
  # Get dates ------------------------------------------------------------------
  obs_dates <- unique(data$date)
  num_dates <- length(obs_dates)
  
  # Get average parameters per date --------------------------------------------
  data_daily <- 
    as.data.frame(
      data %>%
        group_by(date) %>%
        summarise(# mean travel time 
          tt_min = mean(tt),
          # mean lag
          lag = round(tt_min / timestep_min),
          # mean depth
          depth = mean(depth),
          # mean discharge
          discharge = mean(discharge)))
  
  # Pivot dataset wider --------------------------------------------------------
  # Rename location names to s1 and s2
  data$location <- 
    case_match(data$location, 
               location.names[1] ~ "s1", 
               location.names[2] ~ "s2")
  
  # Select only columns of interest, wide-transform dataset
  data_wide <- 
    as.data.frame(data %>%
                    select(date, solar.time, location, DO.obs, DO.sat, 
                           temp.water, light, depth) %>%
                    pivot_wider(names_from = location, 
                                values_from = DO.obs:depth))
  
  # Get number of unique observations per date ---------------------------------
  date_table <- table(data_wide$date)
  num_daily_obs <- unname(date_table)
  
  # For future editing: loop over total number of dates here, separating data 
  # from each date into its own data_ply ---------------------------------------
  data_ply <- data_wide[data_wide$date == ply_date,]
  daily_ply <- data_daily[data_daily$date == ply_date,]
  obs_ply <- date_table[names(date_table) == ply_date]
  
  # Pull day's parameters out of data_daily
  tt_min <- daily_ply$tt_min
  lag <- daily_ply$lag
  depth <- daily_ply$depth
  discharge <- daily_ply$discharge
  moddate <- names(obs_ply)
  nobs <- as.numeric(obs_ply)
  
  # Initialize empty vectors for mean light fraction and water temperature
  lightfrac <- vector(mode = "numeric", length = nobs)
  temp <- vector(mode = "numeric", length = nobs)
  depthvec <- vector(mode = "numeric", length = nobs)
  
  # Loop through observations in day
  for (i in 1:nobs){
    # Mean light -------------------------------------------------------------
    # Calculate fraction of daily light seen at upstream station
    lightstep_s1 <- sum(data_ply$light_s1[i:(i+lag)])
    lighttotal_s1 <- sum(data_ply$light_s1, na.rm = TRUE)
    lightfrac_s1 <- lightstep_s1 / lighttotal_s1
    
    # Calculate fraction of daily light seen at downstream station
    lightstep_s2 <- sum(data_ply$light_s2[i:(i+lag)])
    lighttotal_s2 <- sum(data_ply$light_s2, na.rm = TRUE)
    lightfrac_s2 <- lightstep_s2 / lighttotal_s2
    
    # Calculate mean light fraction, add to vector
    lightfrac[i] <- mean(c(lightfrac_s1, lightfrac_s2), na.rm = TRUE)
    
    # Mean water temperature -------------------------------------------------
    temp[i] <- mean(c(data_ply$temp.water_s1[i:(i+lag)], 
                      data_ply$temp.water_s2[i:(i+lag)]), na.rm = TRUE)
    
    # Mean depth -------------------------------------------------------------
    depthvec[i] <- mean(c(data_ply$depth_s1[i:(i+lag)],
                          data_ply$depth_s2[i:(i+lag)]), na.rm = TRUE)
  }
  
  # Format list of data
  data_list = c(
    list(
      d = d, # see note above
      date = moddate, # date
      timestep = timestep_min, # Length of each timestep in minutes
      n = nobs, # Number of observations on date
      n24 = n24 ,# Number of observations in 24 hours, given timestep
      tt = tt_min / 1440, # Avg travel time on day
      lag = lag, # Time lag between s1 and s2
      #depth = depth, #average depth on day
      discharge = discharge # average discharge on day
    ),
    
    list(
      light_mult_GPP = lightfrac[1:(nobs-lag)],
      lightfrac = lightfrac[1:(nobs-lag)],
      const_mult_ER = rep(1, length = nobs-lag),
      KO2_conv = Kcor_O2(temp = temp, K600 = 1)[1:(nobs-lag)],
      depth = depthvec[1:(nobs-lag)],
      DO_obs_up = data_ply$DO.obs_s1[1:(nobs-lag)],
      DO_obs_down = data_ply$DO.obs_s2,
      DO_sat_up = data_ply$DO.sat_s1[1:(nobs-lag)],
      DO_sat_down = data_ply$DO.sat_s2
    ),
    
    specs
  )
  
  return(data_list)
}
#' \code{prepdata_bayes_twostation} Prepare two-station data for Bayesian modeling
#'
#' @param data Dataset containing columns: `solar.time`, `location`, `DO.obs`, `DO.sat`, `temp.water`, `light`, `depth`, `tt`, `pressure`
#' @param ply_date Date to model, YYYY-MM-DD
#' @param specs List returned by `set_specs()`

prepdata_bayes_twostation <- function(data, specs, up.name, down.name){
  # Grab model name from specs
  modname <- specs$model_name
  
  # Verify that modname is an accepted model name
  if(!(modname %in% c("o2_twostation", "n2_twostation"))){
    stop("model name must be `o2_twostation` or `n2_twostation`")
  }
  
  # Glue location names together
  location.names <- c(up.name, down.name)
  
  # Seperate upstream and downstream data for checking procedures
  updata <- filter(data, location == up.name)
  downdata <- filter(data, location == down.name)
  
  date_table_up <- table(updata$date)
  date_table_down <- table(downdata$date)
  
  # Verify that number of dates is the same
  num_dates <- unique(c(length(date_table_up), length(date_table_down)))
  if(length(num_dates) > 1){
    stop("upstream and downstream sensor have differing number of observation dates")
  }
  
  # Check if dates have same number of observations at upstream and downstream sensor
  if(!identical(date_table_up, date_table_down)){
    stop("upstream and downstream sensor have differing number of observations")
  }
  # Check if each date has same number of observations
  tot_daily_obs <- unique(c(unname(date_table_up), unname(date_table_down)))
  # Return error if there's unequal observations across days
  if(length(tot_daily_obs) > 1){
    stop("dates have differing numbers of rows")
  }
  
  num_daily_obs <- tot_daily_obs
  
  time_by_date_matrix <- function(vec){
    matrix(data = vec, nrow = num_daily_obs, ncol = num_dates, byrow = FALSE)
  }
  
  # Confirm that number of daily observations is the same at both sensor stations
  obs_dates_up <- time_by_date_matrix(as.numeric(updata$date, format = "%Y-%m-%d"))
  obs_dates_down <- time_by_date_matrix(as.numeric(downdata$date, format = "%Y-%m-%d"))
  if(!all.equal(obs_dates_up, obs_dates_down)){
    stop("number of daily observations differs at upstream and downstream locations")
  }
  
  # Confirm that data collected at both sensor stations has the same timestep
  obs_times_up <- time_by_date_matrix(as.numeric(updata$solar.time - updata$solar.time[1], units = "days"))
  obs_times_down <- time_by_date_matrix(as.numeric(downdata$solar.time - downdata$solar.time[1], units = "days"))
  if(!all.equal(obs_times_up, obs_times_down)){
    stop("time of measurement differs at upstream and downstream locations")
  }
  
  # Get timesteps --------------------------------------------------------------
  timestep_up <- as.numeric(diff(obs_times_up), units = "days")
  timestep_down <- as.numeric(diff(obs_times_down), units = "days")
  if(!all.equal(timestep_up, timestep_down)){
    stop("time between successive measurments differs at upstream and downstream locations")
  }
  timestep_eachday <- c(timestep_up, timestep_down)
  if(length(unique(round(timestep_eachday, digits = 10))) != 1){
    stop("could not determine a single timestep for all observations")
  }
  timestep_days <- mean(timestep_eachday)
  timestep_min <- timestep_days * 1440
  n24 <- round(1/timestep_days)
  
  # Get dates ------------------------------------------------------------------
  obs_dates <- unique(data$date)
  num_dates <- length(obs_dates)
  
  # Get average parameters per date --------------------------------------------
  data_daily <- 
    as.data.frame(
      data %>%
        group_by(date) %>%
        summarise(
          # mean travel time 
          tt_min = mean(tt),
          # mean lag
          lag = round(tt_min / timestep_min),
          # mean depth
          depth = mean(depth),
          # mean discharge
          discharge = mean(discharge),
          # number of observations each day, less lag
          n = length(unique(solar.time)) - lag))
  
  # Pivot dataset wider --------------------------------------------------------
  # Rename location names to s1 and s2
  data$location <- 
    case_match(data$location, 
               location.names[1] ~ "s1", 
               location.names[2] ~ "s2")
  
  # Function for adding lag to s2 data points
  lagged <- function(widedata){
    # initialize i
    dayi <- 1
    # Intialize empty list
    laglist <- vector(mode = "list", length = length(unique(widedata$date)))
    
    # Loop
    for(dayi in seq_along(unique(widedata$date))){
      # Get date
      day <- unique(widedata$date)[dayi]
      
      # Segment data_wide to just data from day
      day_data <- filter(widedata, date == day)
      # Grab lag
      lagn <- filter(data_daily, date == day)$lag
      # Grab tot number of observations less lag
      n <- filter(data_daily, date == day)$n
      
      # Mutate 
      # Note that light_s2 does not need lagging
      day_lagged <- 
        day_data %>%
        mutate_at(vars(!starts_with("light") & ends_with("_s2")), 
                  ~lag(., n = lagn))
      
      # Rename columns to denote that they've been lagged
      names(day_lagged) <- gsub("_s2", "_s2.lag", names(day_lagged))
      # Except for light
      names(day_lagged) <- gsub("light_s2.lag", "light_s2", names(day_lagged))
      
      # Trim size of dataframe so that we don't have NA lag rows at the start
      # index should begin at lagn + 1
      day_lagged <- day_lagged[(lagn+1):(n+lagn),]

      # Grab the average light
      day_lagged$light_avg <- rowMeans(day_lagged[c("light_s1", "light_s2")])
      
      # Calculate light fraction
      lighttot <- sum(day_lagged$light_avg, na.rm = TRUE)
      day_lagged$lightfrac <- day_lagged$light_avg / lighttot
      
      # Add to output list
      laglist[[dayi]] <- day_lagged
    }
    
    # Glue output together
    lagframe <- bind_rows(laglist)
    return(lagframe)
  }
  
  # Select only columns of interest and wide-transform dataset
  # & lag s2 columns
  if(modname == "o2_twostation"){
    data_wide <- 
      as.data.frame(
        data %>%
          select(date, solar.time, location, DO.obs, DO.sat, 
                 temp.water, light, depth) %>%
          pivot_wider(names_from = location, 
                      values_from = DO.obs:depth))
    
    # Add lag for s2 columns
    data_wide <- lagged(data_wide)
    
  }
  if(str_detect(modname, "n2_twostation")){
    data_wide <-
      as.data.frame(
        data %>%
          select(date, solar.time, location, N2.obs, N2.sat, 
                 temp.water, light, depth)%>%
          pivot_wider(names_from = location, 
                      values_from = N2.obs:depth))
    
    # Add lag for s2 columns
    data_wide <- lagged(data_wide)
  }
  
  time_by_date_matrix_lagged <- function(vec){
    matrix(data = vec, nrow = nrow(data_wide), ncol = num_dates, byrow = FALSE)
  }
  
  # Grab avg depth
  # Since we've lagged s2, we can calc this with a simple row average
  data_wide$depth_avg <- rowMeans(data_wide[c("depth_s1", "depth_s2.lag")])
  
  # Grab avg water temp
  data_wide$temp.water_avg <- rowMeans(data_wide[c("temp.water_s1", 
                                                   "temp.water_s2.lag")])
  
  # Format list of data
  data_list = c(
    list(
      
      # Overall
      d = num_dates, 
      timestep = timestep_min, # Length of each timestep in minutes
      n24 = n24 # Number of total potential observations in 24 hours
    ),
    
    list(
      # Data that has one value per day
      tt = data_daily$tt_min / 1400,
      lag = data_daily$lag,
      n = data_daily$n
    ),
    
    list(
      
      depth = time_by_date_matrix_lagged(data_wide$depth_avg),
      temp = time_by_date_matrix_lagged(data_wide$temp.water_avg)
      
    ),
    
    specs
  )
  
  # Add O2 data to output
  if(modname == "o2_twostation"){
    data_list$light_mult_GPP <- time_by_date_matrix_lagged(data_wide$lightfrac)
    data_list$const_mult_ER <- time_by_date_matrix_lagged(1)
    
    data_list$KO2_conv_vec <- Kcor_O2(K600 = 1, temp = data_wide$temp.water_avg)
    data_list$KO2_conv <- time_by_date_matrix_lagged(data_list$KO2_conv_vec)
    
    data_list$DO_obs_up <- time_by_date_matrix_lagged(data_wide$DO.obs_s1)
    data_list$DO_obs_down <- time_by_date_matrix_lagged(data_wide$DO.obs_s2.lag)
    data_list$DO_sat_up <- time_by_date_matrix_lagged(data_wide$DO.sat_s1)
    data_list$DO_sat_down <- time_by_date_matrix_lagged(data_wide$DO.sat_s2.lag)
  }
  
  # Add N2 data to output
  if(str_detect(modname, "n2_twostation")){
    data_list$light_mult_N2consume <- time_by_date_matrix_lagged(data_wide$lightfrac)
    data_list$const_mult_DN <- time_by_date_matrix_lagged(1)
    
    data_list$KN2_conv_vec <- Kcor_N2(K600 = 1, temp = data_wide$temp.water_avg)
    data_list$KN2_conv <- time_by_date_matrix_lagged(data_list$KN2_conv_vec)
    
    data_list$N2_obs_up <- time_by_date_matrix_lagged(data_wide$N2.obs_s1)
    data_list$N2_obs_down <- time_by_date_matrix_lagged(data_wide$N2.obs_s2.lag)
    data_list$N2_sat_up <- time_by_date_matrix_lagged(data_wide$N2.sat_s1)
    data_list$N2_sat_down <- time_by_date_matrix_lagged(data_wide$N2.sat_s2.lag)
  }
  
  return(data_list)
}

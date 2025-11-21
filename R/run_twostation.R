#' \code{run_twostation} Fit a two-station Bayesian metabolism model for O2 or N2 flux
#'
#' @param data dataframe containing columns `solar.time`, `location`, `DO.obs`, `DO.sat`, `temp.water`, `light`, `depth`, `tt`, `pressure`
#' @param specs list of model specifications as returned by `set_twostation_specs()`
#' @param upname a character string specifying the name of the upstream sampling station in the `location` column. Defaults to "upstream"
#' @param downname a character string specifying the name of the downstream sampling station in the `location` column. Defaults to "downstream"
#' @param return_fit if TRUE, include dataframe of daily Rhat and R2 and dataframe of modeled gas concentrations in output
#' @import dplyr
#' 
#' @returns List of model outputs
#'
#' @references Raymond, P. A., Zappa, C. J., Butman, D., Bott, T. L., Potter, J., Mulholland, P., et al. (2012). Scaling the gas transfer velocity and hydraulic geometry in streams and small rivers: Gas transfer velocity and hydraulic geometry. Limnology and Oceanography: Fluids and Environments, 2(1), 41–53. https://doi.org/10.1215/21573689-1597669
#'
#' @export
run_twostation <- function(data, specs, upname = "upstream", 
                           downname = "downstream", return_fit = TRUE){
  
  # Add date column to data in case it's not there already
  data$date <- lubridate::date(data$solar.time) 
  
  # To hold warning outputs
  warn_strs <- character(0)
  stop_strs <- character(0)
  
  # Execute the bayes function
  bayes_allday <- withCallingHandlers(
    tryCatch({
      if(is.null(data) || nrow(data) == 0) stop("no valid days of data")
      # Try to fit bayes
      data_list <- prepdata_bayes_twostation(
        data = data, specs = specs, up.name = upname, down.name = downname)
      #do.call(runstan_bayes_twostation, c(list(data_list = data_list), specs))
      bayes_out <- runstan_bayes_twostation(data_list = data_list, specs = specs)
    }, error = function(err){
      # Currently, returning error & warnings are normal - how I have the 
      # function indexing over the lags creates NaNs in the output rows 1:lag,
      # which stan hates. have to fix it but being very lazy for now
      stop_strs <<- c(stop_strs, err$message)
      return(bayes_out)
    }), message = function(war){
      warn_strs <<- c(warn_strs, war$message)
      #invokeRestart("muffleWarning")
    })
  
  # Matching date and time info to indicies
  date_df <- tibble::tibble(
    date = as.Date(unique(data$date)),
    date_index = seq_len(data_list$d)
  )
  
  datetime_df <- tibble::tibble(
    solar.time = unique(data$solar.time)[(1+data_list$lag):(data_list$lag + data_list$n)],
    date_index = rep(seq_len(data_list$d), each = data_list$n),
    time_index = rep(seq_len(data_list$n), times = data_list$d)) %>%
    left_join(date_df, by = "date_index")
  
  # Skipping some safety checks/error handling here. Being dangerous
  
  # Match indicies to estimates
  date_index <- time_index <- index <- "dplyr.var"
  bayes_allday$daily <- 
    bayes_allday$daily %>%
    left_join(date_df, by = "date_index") %>%
    select(-date_index, -time_index, -index) %>%
    select(date, everything())
  bayes_allday$inst <- 
    bayes_allday$inst %>%
    left_join(datetime_df, by = c("date_index", "time_index")) %>%
    select(-date_index, -time_index, -index) %>%
    select(date, solar.time, everything()) 
  
  if(return_fit){
    # Pull rhat and R2 values
    dayfits <- 
      data.frame(date = bayes_allday$daily$date,
                 bayes_allday$daily[str_detect(names(bayes_allday$daily), 
                                               pattern = "daily_Rhat")],
                 bayes_allday$daily[str_detect(names(bayes_allday$daily), 
                                               pattern = "R2_mean")])
    
    # Reformat the observed gas data
    obsgas <- data_list[str_detect(names(data_list), pattern = "obs_down")]
    gasvec <- as.data.frame(unlist(obsgas, use.names = FALSE))
    names(gasvec) <- names(obsgas)
    
    # Pull modeled and measured gas
    gasfits <-
      data.frame(solar.time = bayes_allday$inst$solar.time,
               gasvec,
               bayes_allday$inst[str_detect(names(bayes_allday$inst), 
                                             pattern = "mod_down_mean")],
               bayes_allday$inst[str_detect(names(bayes_allday$inst), 
                                            pattern = "mod_down_sd")]) %>%
      na.omit()
    
    returned_fits <-
      list(fit_stats = dayfits,
           modeled_gas = gasfits)
  }
  
  # Return output to user
  c(bayes_allday,
    list(fit_info = if(return_fit) returned_fits else NULL,
         mcmc_data = if(specs$keep_mcmc_data) data_list else NULL,
         messages = trimws(unique(warn_strs)),
         errors = trimws(unique(stop_strs))))
}

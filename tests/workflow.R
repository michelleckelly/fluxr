library(rstan)
library(devtools)
library(tidyverse)
library(roxygen2)

# Build package
roxygen2::roxygenize()

# Compile stan model
modl <- stanc(file = "inst/stan/o2_twostation.stan")
str(modl)

## Note: for cleaning new datasets: will need to import convert_UTC_to_solartime function from streammetabolizer

# Load internal dataset ########################################################
# Already has equal time breaks, travel time calculated, no NAs
data(diel)

# Pull just the O2 data 
data_all <- 
  diel %>%
  select(solar.time, location, DO.obs, DO.sat, temp.water, light.calc,
         depth, discharge, tt, pressure) %>%
  rename(light = light.calc)

# Add date column
data_all$date <- date(data_all$solar.time)

# Set specifications for modeling ##############################################
mod_specs <- set_specs(model_name = "o2_twostation",
                       verbose = TRUE)

# Pass data and specs to prepping function #####################################
data_prepped <- prepdata_bayes_twostation(data = data_all, 
                                          ply_date = "2024-09-19", 
                                          specs = mod_specs)

# Run model ####################################################################
#runstan_bayes(data_list = data_prepped, specs = data_prepped$specs)

runstan_bayes <- function(data_list){
  # Pull out model name
  model_name <- data_list$model_name
  
  # Pull settings from specs list
  ncores <- data_list$n_cores
  n_chains <- data_list$n_chains
  verbose <- data_list$verbose
  warmup <- data_list$burnin_steps
  iter <- data_list$saved_steps + warmup
  thin_steps <- data_list$thin_steps
  params_out <- data_list$params_out
  
  # Determine number of cores to use
  tot_cores <- parallel::detectCores()
  # Compare to number of cores specified in set_specs()
  n_cores <- min(tot_cores, ncores)
  
  # Tell user how many cores Stan is going to use so they can adjust settings
  # if they want
  if(verbose){message(paste0("\nStan requesting ", n_chains, " chains on ", 
                             n_cores, " of ", tot_cores, " available cores."))}
  
  # Get filepath to compiled & uncompiled stan model
  model_path <- paste0("inst/stan/", model_name, ".stan")
  mobj_path <- paste0("inst/stan/", model_name, ".rds")
  
  # Check for existence of compiled model. If model doesn't exist, compile it.
  if(file.exists(mobj_path)){
    if(verbose){message(paste0("\nUsing previously compiled ", model_name, 
                               " model ", compiled_path))}
  } else{
    # Start clock for compilation time
    start.time <- Sys.time()
    # If file doesn't exist, compile model
    if(verbose){message(paste0("\nCompiling ", model_name, ".stan model."))}
    stan_model(file = model_path, auto_write = TRUE)
    # Stop clock
    end.time <- Sys.time()
    
    # Verify that compilation was successful & present done message
    if(file.exists(mobj_path)){
      elapsed <- round(as.numeric(end.time - start.time, 
                                  units = "hours") * 60 * 60,
                       digits = 2)
      if(verbose){message(paste0("\nDone! Compile time was ", elapsed, 
                                 " seconds."))}
    } else{
      message("ERROR Stan model was not compiled successfully")
    }
  }
  
  # streamMetabolizer has the following garbage collection line with the 
  # comment "this humble line saves us from many horrible R crashes" which
  # sounds good so replicating that here
  gc()
  
  # Read compiled Stan file into memory
  stan_mobj <- readRDS(mobj_path)
  
  # Run Stan
  if(verbose){message("\nSampling Stan model")}
  
  start.time <- Sys.time()
  
  rstan::sampling(
    object = stan_mobj,
    data = data_list, # think i have to do something with linlight$light_mult_GPP here
    pars = params_out, #NOTE haven't set up params_out
    include = TRUE,
    chains = n_chains,
    warmup = warmup,
    iter = iter, 
    thin = thin_steps,
    init = "random",
    verbose = verbose,
    open_progress = FALSE,
    cores = n_cores
  )
  
  end.time <- Sys.time()
  elapsed <- round(as.numeric(end.time - start.time, 
                              units = "hours") * 60, digits = 2)
  
  
  
  
  # run stan
}